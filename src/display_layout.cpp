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
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/detail/sketch/hyperloglog.hpp>

struct config
{
    std::filesystem::path input{};
    std::filesystem::path output{};
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

        if (fill_current_kmers)
        {
            for (auto && [seq] : fin)
            {
                for (uint64_t hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{kmer_size}))
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
                for (uint64_t hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{kmer_size}))
                {
                    sketch.add(reinterpret_cast<char *>(&hash_value), sizeof(hash_value));
                }
            }
        }
    }
}

int execute(config const & cfg)
{
    chopper::configuration chopper_config;
    std::ifstream layout_file{cfg.input};

    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"}; // GCOVR_EXCL_LINE

    std::vector<std::vector<std::string>> filenames = chopper::layout::read_filenames_from(layout_file);
    chopper::layout::read_config_from(chopper_config, layout_file);
    chopper_config.hibf_config.read_from(layout_file);
    seqan::hibf::layout::layout hibf_layout{};
    hibf_layout.read_from(layout_file);
    auto const & hibf_config = chopper_config.hibf_config;

    std::ofstream output_stream{cfg.output};

    if (!output_stream.good() || !output_stream.is_open())
        throw std::logic_error{"Could not open file " + cfg.output.string() + " for reading"}; // GCOVR_EXCL_LINE

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
    parser.add_option(cfg.output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
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

    execute(cfg);
}
