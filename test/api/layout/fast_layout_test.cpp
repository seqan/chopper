// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <format>    // for format
#include <stdexcept> // for invalid_argument
#include <vector>    // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/fast_layout.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace
{

// splitmix64 finaliser, see lsh_distributed_ibf_layout_test.cpp.
uint64_t scramble(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

using diagnostic = seqan::hibf::layout::layout::diagnostic;

auto ignore_unexpected_empty_bins_and_throw = [](diagnostic const & finding)
{
    if (finding.what != diagnostic::code::unexpected_empty_bins)
        throw std::invalid_argument{std::format("{}", finding)};
};

} // namespace

TEST(fast_layout_test, recursion)
{
    // With tmax = 4, a merged bin is laid out recursively with the fast layout if it contains at least 64 * 4 = 256
    // user bins of similar size. 2000 user bins in 4 top-level technical bins give about 500 per merged bin.
    size_t const tmax{4};
    size_t const number_of_user_bins{2000};

    chopper::configuration config{};
    config.fast_layout = true;
    config.hibf_config.tmax = tmax;
    config.hibf_config.threads = 2;
    config.hibf_config.number_of_user_bins = number_of_user_bins;
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement
    config.hibf_config.input_fn = [](size_t const ub, seqan::hibf::insert_iterator it)
    {
        uint64_t const offset = static_cast<uint64_t>(ub) << 32;
        for (uint64_t i = 0; i < 3'000 + ub % 7 * 100; ++i)
            it = scramble(offset + i);
    };

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

    std::vector<size_t> cardinalities(number_of_user_bins);
    for (size_t i = 0; i < number_of_user_bins; ++i)
        cardinalities[i] = sketches[i].estimate();

    std::vector<size_t> positions = seqan::hibf::iota_vector(number_of_user_bins);
    seqan::hibf::layout::layout hibf_layout{};
    chopper::layout::fast_layout(config, positions, cardinalities, sketches, minHash_sketches, hibf_layout);

    EXPECT_NO_THROW(hibf_layout.validate(config.hibf_config, ignore_unexpected_empty_bins_and_throw));
    EXPECT_GE(hibf_layout.number_of_levels(), 2u);
}
