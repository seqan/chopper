// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <chopper/layout/input.hpp>

#include <hibf/build/build_data.hpp>
#include <hibf/config.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include "compute_ibf_size.hpp"
#include "shared.hpp"

struct ibf_stats
{
    size_t size_in_bits{};
    size_t level{};
    size_t max_elements{};
    size_t tbs_too_big{};
    size_t tbs_too_many_elements{};
    double load_factor{};
};

struct per_level_stats
{
    size_t number_of_levels{};
    std::vector<size_t> size_in_bits{};
    std::vector<size_t> max_elements{};
    std::vector<size_t> tbs_too_big{};
    std::vector<size_t> tbs_too_many_elements{};
    std::vector<size_t> num_ibfs{};
    std::vector<double> load_factor{};

    explicit per_level_stats(std::vector<ibf_stats> const & stats)
    {
        number_of_levels = std::ranges::max(stats,
                                            [](ibf_stats const & lhs, ibf_stats const & rhs)
                                            {
                                                return lhs.level < rhs.level;
                                            })
                               .level
                         + 1u;

        size_in_bits.resize(number_of_levels);
        max_elements.resize(number_of_levels);
        tbs_too_big.resize(number_of_levels);
        tbs_too_many_elements.resize(number_of_levels);
        num_ibfs.resize(number_of_levels);
        load_factor.resize(number_of_levels);

        for (size_t i = 0; i < stats.size(); ++i)
        {
            auto const & stat = stats[i];
            size_in_bits[stat.level] += stat.size_in_bits;
            max_elements[stat.level] += stat.max_elements;
            tbs_too_big[stat.level] += stat.tbs_too_big;
            tbs_too_many_elements[stat.level] += stat.tbs_too_many_elements;
            load_factor[stat.level] += stat.load_factor;
            ++num_ibfs[stat.level];
        }

        for (size_t i = 0; i < number_of_levels; ++i)
        {
            assert(num_ibfs[i] > 0u);
            load_factor[i] /= num_ibfs[i];
            max_elements[i] /= num_ibfs[i];
            tbs_too_many_elements[i] /= tbs_too_big[i] ? tbs_too_big[i] : 1u;
        }
    }

    void print(std::ostream & stream, size_t const number_of_user_bins) const
    {
        stream << std::fixed << std::setprecision(2);
        stream << "# Levels: " << number_of_levels << '\n';
        stream << "# User bins: " << number_of_user_bins << '\n';
        stream << "LEVEL\t"                    //
               << "BIT_SIZE\t"                 //
               << "IBFS\t"                     //
               << "AVG_LOAD_FACTOR\t"          //
               << "TBS_TOO_BIG\t"              //
               << "AVG_TBS_TOO_BIG_ELEMENTS\t" //
               << "AVG_MAX_ELEMENTS\n";

        for (size_t idx = 0; idx < number_of_levels; ++idx)
        {
            stream << idx << '\t'                        //
                   << size_in_bits[idx] << '\t'          //
                   << num_ibfs[idx] << '\t'              //
                   << load_factor[idx] * 100 << '\t'     //
                   << tbs_too_big[idx] << '\t'           //
                   << tbs_too_many_elements[idx] << '\t' //
                   << max_elements[idx] << '\n';
        }
    }
};

void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   seqan::hibf::build::build_data const & data,
                   seqan::hibf::layout::layout::user_bin const & record)
{
    data.config.input_fn(record.idx, seqan::hibf::insert_iterator{kmers});
}

size_t hierarchical_stats(std::vector<ibf_stats> & stats,
                          robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          seqan::hibf::layout::graph::node const & current_node,
                          seqan::hibf::build::build_data & data,
                          size_t const current_hibf_level)
{
    size_t const ibf_pos{data.request_ibf_idx()};
    std::vector<size_t> size_waste_per_tb(current_node.number_of_technical_bins);

    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [&]() -> size_t
    {
        if (current_node.favourite_child_idx.has_value()) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            hierarchical_stats(stats,
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
    std::atomic<uint64_t> tbs_too_big{};
    std::atomic<uint64_t> tbs_too_many_elements{};
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    if (!current_node.children.empty())
    {
        std::mutex merge_kmer_mutex{};

        size_t number_of_threads{};
        std::vector<size_t> indices(current_node.children.size());
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
            auto & child = current_node.children[index];

            robin_hood::unordered_flat_set<uint64_t> local_kmers{};
            hierarchical_stats(stats, local_kmers, child, data, current_hibf_level + 1u);

            if (current_hibf_level > 0 /* not top level */)
            {
                std::lock_guard<std::mutex> guard{merge_kmer_mutex};
                update_parent_kmers(parent_kmers, local_kmers);
            }

            if (amount_of_max_bin_kmers < local_kmers.size())
            {
                ++tbs_too_big;
                tbs_too_many_elements += local_kmers.size() - amount_of_max_bin_kmers;
            }
            else
            {
                size_waste_per_tb[child.parent_bin_index] = amount_of_max_bin_kmers - local_kmers.size();
            }
        }
    }

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node.favourite_child_idx.has_value()) ? 0u : 1u};
    for (size_t i = start; i < current_node.remaining_records.size(); ++i)
    {
        auto const & record = current_node.remaining_records[i];

        compute_kmers(kmers, data, record);

        size_t const kmers_per_tb = kmers.size() / record.number_of_technical_bins + 1u;

        if (amount_of_max_bin_kmers < kmers_per_tb)
        {
            tbs_too_big += record.number_of_technical_bins;
            tbs_too_many_elements += (kmers_per_tb - amount_of_max_bin_kmers) * record.number_of_technical_bins;
        }
        else
        {
            assert(size_waste_per_tb.size() >= record.storage_TB_id + record.number_of_technical_bins);
            for (size_t i = record.storage_TB_id; i < record.storage_TB_id + record.number_of_technical_bins; ++i)
                size_waste_per_tb[i] = amount_of_max_bin_kmers - kmers_per_tb;
        }

        if (current_hibf_level > 0 /* not top level */)
            update_parent_kmers(parent_kmers, kmers);

        kmers.clear();
    }

    double const load_factor = [&]()
    {
        size_t const waste = std::accumulate(size_waste_per_tb.begin(), size_waste_per_tb.end(), size_t{});
        size_t const total = amount_of_max_bin_kmers * current_node.number_of_technical_bins;
        assert(waste <= total);
        assert(total > 0u);
        return (total - waste) / static_cast<double>(total);
    }();

    auto & current_stats = stats[ibf_pos];
    current_stats.size_in_bits = ibf_size;
    current_stats.level = current_hibf_level;
    current_stats.max_elements = amount_of_max_bin_kmers;
    current_stats.tbs_too_big = tbs_too_big.load();
    current_stats.tbs_too_many_elements = tbs_too_many_elements.load();
    current_stats.load_factor = load_factor;

    return ibf_pos;
}

size_t hierarchical_stats(std::vector<ibf_stats> & stats,
                          seqan::hibf::layout::graph::node const & root_node,
                          seqan::hibf::build::build_data & data)
{
    robin_hood::unordered_flat_set<uint64_t> root_kmers{};
    return hierarchical_stats(stats, root_kmers, root_node, data, 0);
}

void execute_general_stats(config const & cfg)
{
    // Read config and layout
    std::ifstream layout_file{cfg.input};
    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"};

// https://godbolt.org/z/PeKnxzjn1
#if defined(__clang__)
    auto tuple = chopper::layout::read_layouts_file(layout_file);
    // https://godbolt.org/z/WoWf55KPb
    auto filenames = std::move(std::get<0>(tuple));
    auto chopper_config = std::move(std::get<1>(tuple));
    auto hibf_layouts = std::move(std::get<2>(tuple));
#else
    auto [filenames, chopper_config, hibf_layouts] = chopper::layout::read_layouts_file(layout_file);
#endif

    // Prepare configs
    chopper_config.hibf_config.threads = cfg.threads;
    auto input_lambda = [&filenames, &chopper_config](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        std::vector<uint64_t> current_kmers;

        for (std::string const & filename : filenames[user_bin_id])
            process_file(filename, current_kmers, chopper_config.k, chopper_config.window_size);

        for (auto const kmer : current_kmers)
            it = kmer;
    };
    chopper_config.hibf_config.input_fn = input_lambda;
    auto const & hibf_config = chopper_config.hibf_config;

    // Prepare stats
    assert(hibf_layouts.size() > 0);
    // size_t part = (hibf_layouts.size() == 1) ? 0 : 1;
    for (auto const & hibf_layout : hibf_layouts)
    {
        size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1u;
        std::vector<ibf_stats> stats(number_of_ibfs);

        // Prepare data
        seqan::hibf::build::build_data data{.config = hibf_config, .ibf_graph = {hibf_layout}};
        seqan::hibf::layout::graph::node const & root_node = data.ibf_graph.root;
        size_t const t_max{root_node.number_of_technical_bins};
        data.fpr_correction = seqan::hibf::layout::compute_fpr_correction(
            {.fpr = hibf_config.maximum_fpr, .hash_count = hibf_config.number_of_hash_functions, .t_max = t_max});

        // Get stats
        hierarchical_stats(stats, root_node, data);

        // Get stats per level
        per_level_stats const level_stats{stats};

        // Output
        std::ofstream output_stream(cfg.output.string(), std::ios::app);
        if (!output_stream.good() || !output_stream.is_open())
            throw std::logic_error{"Could not open file " + cfg.output.string() + " for reading (appending)."};

        level_stats.print(output_stream, hibf_layout.user_bins.size());
        // ++part;
    }
}

void execute_sizes(config const & cfg)
{
    execute_general_stats(cfg);
}
