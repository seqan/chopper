// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <iostream>

#include <robin_hood.h>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include "shared.hpp"

struct stats
{
    std::vector<size_t> ibf_sizes{};
    std::vector<size_t> ibf_levels{};
    std::vector<double> ibf_load_factor{};
};

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
                        size_t const current_hibf_level)
{
    size_t const kmers_per_bin = std::ceil(static_cast<double>(kmers.size()) / number_of_bins);
    size_t const bin_size =
        std::ceil(seqan::hibf::build::bin_size_in_bits({.fpr = data.config.maximum_false_positive_rate,
                                                        .hash_count = data.config.number_of_hash_functions,
                                                        .elements = kmers_per_bin})
                  * data.fpr_correction[number_of_bins]);

    size_t const ibf_size = ibf_node.number_of_technical_bins * bin_size;

    if (current_hibf_level > 0 /* not top level */)
        update_parent_kmers(parent_kmers, kmers);

    return ibf_size;
}

size_t hierarchical_stats(stats & stats_data,
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

            robin_hood::unordered_flat_set<uint64_t> kmers2{};
            hierarchical_stats(stats_data, kmers2, child, data, current_hibf_level + 1u);

            if (current_hibf_level > 0 /* not top level */)
            {
                std::lock_guard<std::mutex> guard{merge_kmer_mutex};
                update_parent_kmers(parent_kmers, kmers2);
            }

            size_waste_per_tb[child.parent_bin_index] =
                amount_of_max_bin_kmers - std::min(amount_of_max_bin_kmers, kmers2.size());
        }
    }

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node.favourite_child_idx.has_value()) ? 0u : 1u};
    for (size_t i = start; i < current_node.remaining_records.size(); ++i)
    {
        auto const & record = current_node.remaining_records[i];

        compute_kmers(kmers, data, record);

        size_t const kmers_per_tb = kmers.size() / record.number_of_technical_bins + 1u;
        size_t const diff = amount_of_max_bin_kmers - std::min(amount_of_max_bin_kmers, kmers_per_tb);

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

    assert(waste <= total);

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
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"};

    auto [filenames, chopper_config, hibf_layout] = chopper::layout::read_layout_file(layout_file);
    chopper_config.hibf_config.threads = cfg.threads;

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

    std::ofstream output_stream{cfg.output};

    if (!output_stream.good() || !output_stream.is_open())
        throw std::logic_error{"Could not open file " + cfg.output.string() + " for reading"};

    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1u;

    stats stats_data;
    stats_data.ibf_sizes.resize(number_of_ibfs);
    stats_data.ibf_levels.resize(number_of_ibfs);
    stats_data.ibf_load_factor.resize(number_of_ibfs);

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
    std::vector<size_t> size_per_level(number_of_levels);
    std::vector<double> load_factor_per_level(number_of_levels);
    std::vector<size_t> num_ibfs_per_level(number_of_levels);

    for (size_t i = 0; i < number_of_ibfs; ++i)
    {
        size_t const ibf_level{stats_data.ibf_levels[i]};
        size_per_level[ibf_level] += stats_data.ibf_sizes[i];
        load_factor_per_level[ibf_level] += stats_data.ibf_load_factor[i];
        ++num_ibfs_per_level[ibf_level];
    }

    for (size_t i = 0; i < number_of_levels; ++i)
        load_factor_per_level[i] /= num_ibfs_per_level[i];

    output_stream << std::fixed << std::setprecision(2);
    output_stream << "# Levels: " << number_of_levels << '\n';
    output_stream << "LEVEL\tBIT_SIZE\tIBFS\tLOAD_FACTOR\n";
    for (size_t i = 0; i < number_of_levels; ++i)
        output_stream << i << '\t' << size_per_level[i] << '\t' << num_ibfs_per_level[i] << '\t'
                      << load_factor_per_level[i] * 100 << '\n';
}

void execute_sizes(config const & cfg)
{
    execute_general_stats(cfg);
}
