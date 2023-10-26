// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cstdint>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/layout/graph.hpp>

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
