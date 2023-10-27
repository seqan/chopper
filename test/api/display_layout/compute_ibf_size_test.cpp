// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <compute_ibf_size.hpp>

#include <hibf/build/construct_ibf.hpp> // for construct_ibf
#include <hibf/layout/compute_fpr_correction.hpp>

TEST(compute_ibf_size_test, merged_bin_is_max_bin)
{
    size_t const number_of_bins = 1;
    size_t const current_hibf_level = 0;

    robin_hood::unordered_flat_set<uint64_t> parent_kmers;
    robin_hood::unordered_flat_set<uint64_t> kmers;

    for (size_t i = 0; i < 1000; ++i)
        kmers.insert(i);

    seqan::hibf::layout::graph::node ibf_node{
        .children = {seqan::hibf::layout::graph::node{}}, // one merged bin
        .parent_bin_index = 0,
        .max_bin_index = 0,
        .number_of_technical_bins = 64,
        .favourite_child_idx = 0, // max bin is a merged bin and it is the first in `children`
        .remaining_records = {}   // not needed for compute_ibf_size
    };

    double fpr = 0.0001;
    size_t hash = 2;
    size_t tmax = 64;

    seqan::hibf::build::build_data data{
        .config = {.number_of_user_bins = 123,
                   .number_of_hash_functions = hash,
                   .maximum_fpr = fpr,
                   .threads = 1,
                   .tmax = tmax},
        .fpr_correction = seqan::hibf::layout::compute_fpr_correction({.fpr = fpr, .hash_count = hash, .t_max = tmax})};

    auto ibf = seqan::hibf::build::construct_ibf(parent_kmers, kmers, number_of_bins, ibf_node, data, true);

    auto ibf_size = compute_ibf_size(parent_kmers, kmers, number_of_bins, ibf_node, data, current_hibf_level);

    EXPECT_EQ(ibf_size, ibf.bit_size());
    EXPECT_TRUE(parent_kmers.empty());
}

TEST(compute_ibf_size_test, split_bin_is_max_bin)
{
    size_t const number_of_bins = 4;
    size_t const current_hibf_level = 1;

    robin_hood::unordered_flat_set<uint64_t> parent_kmers;
    robin_hood::unordered_flat_set<uint64_t> kmers;

    for (size_t i = 0; i < 1000; ++i)
        kmers.insert(i);

    seqan::hibf::layout::graph::node ibf_node{
        .children = {seqan::hibf::layout::graph::node{}}, // one merged bin
        .parent_bin_index = 0,
        .max_bin_index = 0,
        .number_of_technical_bins = 64,
        .favourite_child_idx = std::nullopt, // max bin is a split bin
        .remaining_records = {}              // not needed for compute_ibf_size
    };

    double fpr = 0.0001;
    size_t hash = 2;
    size_t tmax = 64;

    seqan::hibf::build::build_data data{
        .config = {.number_of_user_bins = 123,
                   .number_of_hash_functions = hash,
                   .maximum_fpr = fpr,
                   .threads = 1,
                   .tmax = tmax},
        .fpr_correction = seqan::hibf::layout::compute_fpr_correction({.fpr = fpr, .hash_count = hash, .t_max = tmax})};

    auto ibf = seqan::hibf::build::construct_ibf(parent_kmers, kmers, number_of_bins, ibf_node, data, true);

    auto ibf_size = compute_ibf_size(parent_kmers, kmers, number_of_bins, ibf_node, data, current_hibf_level);

    EXPECT_EQ(ibf_size, ibf.bit_size());
    EXPECT_FALSE(parent_kmers.empty());
}
