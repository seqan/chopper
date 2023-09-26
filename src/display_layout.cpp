// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <iostream>
#include <set>

#include <robin_hood.h>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <chopper/adjust_seed.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/sketch/hyperloglog.hpp>

struct config
{
    std::filesystem::path input{};
    std::filesystem::path output_prefix{};
};

struct stats
{
    std::vector<size_t> ibf_sizes{};
    std::vector<size_t> ibf_levels{};
    std::vector<double> ibf_load_factor{};
};

static void print_progress(size_t const percentage)
{
    assert(percentage <= 100u);
    std::cerr << '[';
    for (size_t i{}; i < percentage; ++i)
        std::cerr << '=';
    if (percentage < 100u)
    {
        std::cerr << '>';
        for (size_t i{1u}; i < 100u - percentage; ++i)
            std::cerr << ' ';
    }
    std::cerr << "] " << percentage << " %\r" << std::flush;
}

struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<dna4_traits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

// Using two sets and erasing from shared is slower
void keep_duplicates(robin_hood::unordered_set<uint64_t> & shared, std::vector<uint64_t> const & current)
{
    // static + calling clear is slower
    robin_hood::unordered_set<uint64_t> result{};

    for (uint64_t value : current)
    {
        if (shared.contains(value))
            result.emplace(value);
    }

    shared = std::move(result);
}

void process_file(std::string const & filename,
                  std::vector<uint64_t> & current_kmers,
                  seqan::hibf::sketch::hyperloglog & sketch,
                  bool const fill_current_kmers,
                  uint8_t const kmer_size)
{
    if (filename.ends_with(".minimiser"))
    {
        uint64_t hash{};
        char * const hash_data{reinterpret_cast<char *>(&hash)};
        std::streamsize const hash_bytes{sizeof(hash)};

        std::ifstream infile{filename, std::ios::binary};

        if (fill_current_kmers)
        {
            while (infile.read(hash_data, hash_bytes))
            {
                current_kmers.push_back(hash);
                sketch.add(hash_data, hash_bytes);
            }
        }
        else
        {
            while (infile.read(hash_data, hash_bytes))
            {
                sketch.add(hash_data, hash_bytes);
            }
        }
    }
    else
    {
        sequence_file_type fin{filename};

        seqan3::shape shape{seqan3::ungapped{kmer_size}};
        auto minimizer_view = seqan3::views::minimiser_hash(shape,
                                                            seqan3::window_size{kmer_size},
                                                            seqan3::seed{chopper::adjust_seed(shape.count())});
        if (fill_current_kmers)
        {
            for (auto && [seq] : fin)
            {
                for (uint64_t hash_value : seq | minimizer_view)
                {
                    current_kmers.push_back(hash_value);
                    sketch.add(reinterpret_cast<char *>(&hash_value), sizeof(hash_value));
                }
            }
        }
        else
        {
            for (auto && [seq] : fin)
            {
                for (uint64_t hash_value : seq | minimizer_view)
                {
                    sketch.add(reinterpret_cast<char *>(&hash_value), sizeof(hash_value));
                }
            }
        }
    }
}

int execute(config const & cfg)
{
    std::ifstream layout_file{cfg.input};

    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"}; // GCOVR_EXCL_LINE

    auto [filenames, chopper_config, hibf_layout] = chopper::layout::read_layout_file(layout_file);
    auto const & hibf_config = chopper_config.hibf_config;

    std::filesystem::path const output_path = [&cfg]()
    {
        std::filesystem::path output_path{cfg.output_prefix};
        output_path += ".stats";
        return output_path;
    }();
    std::ofstream output_stream{output_path};

    if (!output_stream.good() || !output_stream.is_open())
        throw std::logic_error{"Could not open file " + output_path.string() + " for reading"}; // GCOVR_EXCL_LINE

    // Fetch all file sizes such that sorting by file size doesn't have to access the filesystem too often.
    // n = filenames.size()
    // Constructing this vector has `n` filesystem accesses.
    // Sorting without pre-fetching has `O(n * log(n))` accesses.
    std::vector<std::uintmax_t> const filesizes{[&filenames]()
                                                {
                                                    std::vector<std::uintmax_t> result{};
                                                    result.reserve(filenames.size());
                                                    for (auto const & filename : filenames)
                                                        result.push_back(std::filesystem::file_size(filename.front()));
                                                    return result;
                                                }()};

    // Sorts by the technical bin indices in the top-level IBF:
    // split bins: storage_TB_id, previous_TB_indices is empty
    // merged bins: previous_TB_indices[0]
    // If the index is the same, sort by file sizes (happens for merged bins).
    // Using the smallest file to initialise the shared k-mers later will be less work.
    std::ranges::sort(
        hibf_layout.user_bins,
        [&filesizes](seqan::hibf::layout::layout::user_bin const & lhs,
                     seqan::hibf::layout::layout::user_bin const & rhs)
        {
            size_t const first_idx = lhs.previous_TB_indices.empty() ? lhs.storage_TB_id : lhs.previous_TB_indices[0];
            size_t const second_idx = rhs.previous_TB_indices.empty() ? rhs.storage_TB_id : rhs.previous_TB_indices[0];
            return first_idx < second_idx || (first_idx == second_idx && filesizes[lhs.idx] < filesizes[rhs.idx]);
        });

    seqan::hibf::sketch::hyperloglog sketch{hibf_config.sketch_bits};
    robin_hood::unordered_set<uint64_t> shared_kmers{};
    // We can't use `shared_kmers.size() == 0` instead of `shared_kmers_initialised`, because keep_duplicates
    // will result in a size of 0 when there are no shared k-mers.
    bool shared_kmers_initialised{false};
    std::vector<uint64_t> current_kmers{};
    size_t ub_count{};    // How many user bins are stored in the current technical bin? Always 1 for split bins.
    size_t split_count{}; // Into how many techincal bins is the user bin split? Always 1 for merged bins.

    std::vector<chopper::layout::hibf_statistics::bin_kind> bin_kinds(
        hibf_config.tmax,
        chopper::layout::hibf_statistics::bin_kind::split);

    for (auto const & max_bin : hibf_layout.max_bins)
    {
        // max_bin.previous_TB_indices.size() == 1: true for merged bins, false for split bins
        // max_bin.previous_TB_indices[0]: technical bin index on the top-level
        if (max_bin.previous_TB_indices.size() == 1)
        {
            bin_kinds[max_bin.previous_TB_indices[0]] = chopper::layout::hibf_statistics::bin_kind::merged;
        }
    }

    size_t const total_ub{hibf_layout.user_bins.size()};               // For progress bar
    size_t const ub_percentage{std::max<size_t>(1u, total_ub / 100u)}; // For progress bar, 1 % of all user bins

    size_t current_idx{}; // The current top-level technical bin index

    // Stats file header
    output_stream << "# Layout: " << cfg.input.c_str() << '\n' //
                  << "tb_index\t"
                  << "size\t"
                  << "shared_size\t"
                  << "ub_count\t"
                  << "kind\t"
                  << "splits" << '\n';

    auto print_result_line = [&]()
    {
        bool const is_merged{bin_kinds[current_idx] == chopper::layout::hibf_statistics::bin_kind::merged};
        size_t const avg_kmer_count = (sketch.estimate() + split_count - 1u) / split_count;

        for (size_t i{}, total{split_count}; i < total; ++i)
        {
            output_stream << current_idx + i << '\t'                  //
                          << avg_kmer_count << '\t'                   //
                          << shared_kmers.size() << '\t'              //
                          << ub_count << '\t'                         //
                          << (is_merged ? "merged" : "split") << '\t' //
                          << split_count << '\n';
            split_count = 0u; // Subsequent split bins display 0, the first split bin displays the actual split count.
        }
    };

    // Iterate over all user bins. They are sorted by their technical bin index in the top-level.
    // Hence, we can process one top-level technical bin, print the results, clear our stored data and continue with
    // the next.
    for (size_t ub_index{}; ub_index < hibf_layout.user_bins.size(); ++ub_index)
    {
        if (ub_index % ub_percentage == 0)
            print_progress(ub_index / ub_percentage);

        auto const & user_bin = hibf_layout.user_bins[ub_index];
        current_kmers.clear();

        // The top-level technical bin index for the current user bin.
        // user_bin.previous_TB_indices.size() == 0: true for split bins, false for merged bins
        // user_bin.storage_TB_id: technical bin index on the lowest level
        // user_bin.previous_TB_indices[0]: technical bin index on the top-level
        size_t const idx =
            (user_bin.previous_TB_indices.size() == 0) ? user_bin.storage_TB_id : user_bin.previous_TB_indices[0];

        // We processed all user bins that belong to the `current_idx`th top-level technical bin.
        // Print results and update data.
        if (idx != current_idx)
        {
            print_result_line();
            sketch.clear();
            shared_kmers.clear();
            shared_kmers_initialised = false;
            ub_count = 0u;
            split_count = 0u;
            current_idx = idx;
        }

        bool const is_merged = bin_kinds[idx] == chopper::layout::hibf_statistics::bin_kind::merged;
        // For user bins in a merged bin, `user_bin.number_of_technical_bins` is the number of technical bins on
        // the lowest level. A user bin could be part of a merged bin on the top-level, but still be split into
        // `user_bin.number_of_technical_bins` many bins on the lowest level.
        split_count = is_merged ? 1u : user_bin.number_of_technical_bins;

        // We don't need to keep the current_kmers if there are no shared k-mers to merge them with.
        bool const fill_current_kmers = is_merged && !(shared_kmers_initialised && shared_kmers.empty());

        for (auto const & filename : filenames[user_bin.idx])
        {
            ++ub_count; // This assumes that each user bin has exactly one associated file. Currently the case.

            process_file(filename, current_kmers, sketch, fill_current_kmers, chopper_config.k);
        }

        // Compute set intersection: shared_kmers = shared_kmers ∩ current_kmers
        // This happens for each user bin that belongs to a merged bin.
        if (fill_current_kmers)
        {
            if (!shared_kmers_initialised)
            {
                shared_kmers_initialised = true;
                shared_kmers.insert(current_kmers.begin(), current_kmers.end());
            }
            else
            {
                keep_duplicates(shared_kmers, current_kmers);
            }
        }
    }

    // Print results for the last top-level technical bin.
    print_result_line();

    // The progress bar uses a carriage return '\r' to only use a single line.
    std::cerr << '\n';

    return 0;
}

void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   seqan::hibf::build::build_data const & data,
                   seqan::hibf::layout::layout::user_bin const & record)
{
    data.config.input_fn(record.idx, std::inserter(kmers, kmers.begin()));
}

void update_parent_kmers(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                         robin_hood::unordered_flat_set<uint64_t> const & kmers)
{
    parent_kmers.insert(kmers.begin(), kmers.end());
}

size_t compute_ibf_size(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                        robin_hood::unordered_flat_set<uint64_t> & kmers,
                        size_t const number_of_bins,
                        seqan::hibf::layout::graph::node const & ibf_node,
                        seqan::hibf::build::build_data & data,
                        size_t current_hibf_level)
{
    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(kmers.size()) / number_of_bins))};
    double const bin_bits{
        static_cast<double>(seqan::hibf::build::bin_size_in_bits({.fpr = data.config.maximum_false_positive_rate,
                                                                  .hash_count = data.config.number_of_hash_functions,
                                                                  .elements = kmers_per_bin}))};
    seqan::hibf::bin_size const bin_size{
        static_cast<size_t>(std::ceil(bin_bits * data.fpr_correction[number_of_bins]))};
    seqan::hibf::bin_count const bin_count{ibf_node.number_of_technical_bins};

    size_t const bin_words = (bin_count.value + 63) >> 6; // = ceil(bins/64)
    size_t const technical_bins = bin_words << 6;         // = bin_words * 64

    size_t const ibf_size = technical_bins * bin_size.value;

    if (current_hibf_level > 0 /* not top level */)
        update_parent_kmers(parent_kmers, kmers);

    return ibf_size;
}

size_t hierarchical_stats(stats & stats_data,
                          robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          seqan::hibf::layout::graph::node const & current_node,
                          seqan::hibf::build::build_data & data,
                          size_t current_hibf_level)
{
    size_t const ibf_pos{data.request_ibf_idx()};
    std::vector<size_t> size_waste_per_tb(current_node.number_of_technical_bins);

    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [&]() -> size_t
    {
        if (current_node.favourite_child_idx.has_value()) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            hierarchical_stats(stats_data,
                               kmers,
                               current_node.children[current_node.favourite_child_idx.value()],
                               data,
                               current_hibf_level + 1);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = current_node.remaining_records[0];
            compute_kmers(kmers, data, record);

            return record.number_of_technical_bins;
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers();
    size_t const ibf_size = compute_ibf_size(parent_kmers, kmers, max_bin_tbs, current_node, data, current_hibf_level);
    size_t const amount_of_max_bin_kmers = kmers.size();
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    auto loop_over_children = [&]()
    {
        if (current_node.children.empty())
            return;

        std::vector<seqan::hibf::layout::graph::node> children = current_node.children; // copy for threads

        size_t const number_of_mutex = (current_node.number_of_technical_bins + 63) / 64;
        std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

        size_t number_of_threads{};
        std::vector<size_t> indices(children.size());
        std::iota(indices.begin(), indices.end(), size_t{});

        // We do not want to process the favourite child. It has already been processed prior.
        // https://godbolt.org/z/6Yav7hrG1
        if (current_node.favourite_child_idx.has_value())
            std::erase(indices, current_node.favourite_child_idx.value());

        if (current_hibf_level == 0)
        {
            // Shuffle indices: More likely to not block each other. Optimal: Interleave
            std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
            number_of_threads = data.config.threads;
        }
        else
        {
            number_of_threads = 1u;
        }

#pragma omp parallel for schedule(dynamic, 1) num_threads(number_of_threads)
        for (auto && index : indices)
        {
            auto & child = children[index];

            robin_hood::unordered_flat_set<uint64_t> kmers{};
            hierarchical_stats(stats_data, kmers, child, data, current_hibf_level + 1);
            auto parent_bin_index = child.parent_bin_index;

            {
                size_t const mutex_id{parent_bin_index / 64};
                std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};

                if (current_hibf_level > 0 /* not top level */)
                    update_parent_kmers(parent_kmers, kmers);

                assert(amount_of_max_bin_kmers >= kmers.size());
                size_waste_per_tb[parent_bin_index] =
                    amount_of_max_bin_kmers - std::min(amount_of_max_bin_kmers, kmers.size());
            }
        }
    };

    loop_over_children();

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node.favourite_child_idx.has_value()) ? 0u : 1u};
    for (size_t i = start; i < current_node.remaining_records.size(); ++i)
    {
        auto const & record = current_node.remaining_records[i];

        compute_kmers(kmers, data, record);

        size_t const kmers_per_tb = kmers.size() / record.number_of_technical_bins;
        size_t const diff = amount_of_max_bin_kmers - std::min(amount_of_max_bin_kmers, kmers_per_tb);

        if (diff > amount_of_max_bin_kmers)
            throw std::runtime_error{"Wrong? diff:" + std::to_string(diff)
                                     + " amount_of_max_bin_kmers:" + std::to_string(amount_of_max_bin_kmers)};

        // size_t const kmers_per_tb_corrected = kmers_per_tb - (0.03 * kmers_per_tb);
        // if (amount_of_max_bin_kmers < kmers_per_tb_corrected)
        // {
        //     throw std::runtime_error{"amount_of_max_bin_kmers:" + std::to_string(amount_of_max_bin_kmers)
        //                              + " was smaller than current size: " + std::to_string(kmers.size()) + "/"
        //                              + std::to_string(record.number_of_technical_bins) + "="
        //                              + std::to_string(kmers.size() / record.number_of_technical_bins)};
        // }

        assert(size_waste_per_tb.size() >= record.storage_TB_id + record.number_of_technical_bins);
        for (size_t i = record.storage_TB_id; i < record.storage_TB_id + record.number_of_technical_bins; ++i)
            size_waste_per_tb[i] = diff;

        if (current_hibf_level > 0 /* not top level */)
            update_parent_kmers(parent_kmers, kmers);

        kmers.clear();
    }

    stats_data.ibf_sizes[ibf_pos] = ibf_size;
    stats_data.ibf_levels[ibf_pos] = current_hibf_level;
    // compute load factor
    size_t const waste = std::accumulate(size_waste_per_tb.begin(), size_waste_per_tb.end(), size_t{});
    size_t const total = amount_of_max_bin_kmers * size_waste_per_tb.size();

    if (waste > total)
        throw std::runtime_error{"Wrong? waste:" + std::to_string(waste) + " total:" + std::to_string(total)};

    stats_data.ibf_load_factor[ibf_pos] = (total - waste) / static_cast<double>(total);

    return ibf_pos;
}

size_t hierarchical_stats(stats & stats_data,
                          seqan::hibf::layout::graph::node const & root_node,
                          seqan::hibf::build::build_data & data)
{
    robin_hood::unordered_flat_set<uint64_t> root_kmers{};
    return hierarchical_stats(stats_data, root_kmers, root_node, data, 0);
}

void execute_general_stats(config const & cfg)
{
    std::ifstream layout_file{cfg.input};

    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"}; // GCOVR_EXCL_LINE

    auto [filenames, chopper_config, hibf_layout] = chopper::layout::read_layout_file(layout_file);

    auto input_lambda = [&filenames, &chopper_config](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        std::vector<uint64_t> current_kmers;
        seqan::hibf::sketch::hyperloglog sketch;

        if (filenames[user_bin_id].size() > 1)
            throw std::runtime_error{"No multi files accepted yet."};

        process_file(filenames[user_bin_id][0], current_kmers, sketch, true, chopper_config.k);

        for (auto const kmer : current_kmers)
            it = kmer;
    };
    chopper_config.hibf_config.input_fn = input_lambda;

    auto const & hibf_config = chopper_config.hibf_config;

    std::filesystem::path const output_path = [&cfg]()
    {
        std::filesystem::path output_path{cfg.output_prefix};
        output_path += ".general_stats";
        return output_path;
    }();
    std::ofstream output_stream{output_path};

    if (!output_stream.good() || !output_stream.is_open())
        throw std::logic_error{"Could not open file " + output_path.string() + " for reading"}; // GCOVR_EXCL_LINE

    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1u;

    stats stats_data;
    stats_data.ibf_sizes.resize(number_of_ibfs);
    stats_data.ibf_levels.resize(number_of_ibfs);
    stats_data.ibf_load_factor.resize(number_of_ibfs, 1.0);

    seqan::hibf::build::build_data data{.config = hibf_config, .ibf_graph = {hibf_layout}};

    seqan::hibf::layout::graph::node const & root_node = data.ibf_graph.root;

    size_t const t_max{root_node.number_of_technical_bins};
    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = hibf_config.maximum_false_positive_rate,
                                                     .hash_count = hibf_config.number_of_hash_functions,
                                                     .t_max = t_max});

    hierarchical_stats(stats_data, root_node, data);

    // accumulate statistics per level
    size_t const number_of_levels = std::ranges::max(stats_data.ibf_levels) + 1u;
    std::vector<size_t> size_per_level(number_of_levels, 0u);
    std::vector<double> load_factor_per_level(number_of_levels, 0.0);
    std::vector<size_t> num_ibfs_per_level(number_of_levels, 0u);

    for (size_t i = 0; i < stats_data.ibf_sizes.size(); ++i)
    {
        size_t const ibf_level{stats_data.ibf_levels[i]};
        size_per_level[ibf_level] += stats_data.ibf_sizes[i];
        load_factor_per_level[ibf_level] += stats_data.ibf_load_factor[i];
        ++num_ibfs_per_level[ibf_level];
    }

    for (size_t i = 0; i < number_of_levels; ++i)
        load_factor_per_level[i] /= num_ibfs_per_level[i];

    output_stream << "#Levels: " << number_of_levels << '\n';
    output_stream << "LEVEL-IDX\tSIZE-IN-BITS\tNUMBER-OF-IBFS\tLOAD-FACTOR-IN-%\n";
    for (size_t i = 0; i < number_of_levels; ++i)
        output_stream << i << '\t' << size_per_level[i] << '\t' << num_ibfs_per_level[i] << '\t'
                      << load_factor_per_level[i] << '\n';
}

inline void set_up_parser(sharg::parser & parser, config & cfg)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Compute a top-level HIBF layout figure file";

    parser.info.description.emplace_back("Computes an table to display the top-level layout.");

    parser.add_subsection("Main options:");
    parser.add_option(
        cfg.input,
        sharg::config{.short_id = '\0',
                      .long_id = "input",
                      .description = "The input must be a layout file computed via chopper layout or raptor layout. ",
                      .required = true,
                      .validator = sharg::input_file_validator{}});
    parser.add_option(cfg.output_prefix,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output-prefix",
                                    .description = "The output-prefix. Create .stats and .general_stats file. ",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"layout_stats", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    config cfg{};
    set_up_parser(parser, cfg);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    // execute(cfg);
    execute_general_stats(cfg);
}
