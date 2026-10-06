// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <algorithm> // for max, ranges::max
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <map>       // for map
#include <numeric>   // for iota
#include <set>       // for set
#include <vector>    // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/fast_layout.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/misc/next_multiple_of_64.hpp>
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

// Checks that `hibf_layout` describes a consistent HIBF over `number_of_user_bins` user bins:
// * user_bins[i].idx == i
// * In every IBF (identified by its path of merged bin indices), each technical bin is either a merged bin (it is on the
//   path of some user bin) or holds exactly one user bin.
// * max_bins contains exactly one entry per lower-level IBF.
// * All technical bin indices are below next_multiple_of_64(tmax): The lowest levels of the DP layout use
//   next_multiple_of_64(#user bins) technical bins, which is only bounded by tmax if tmax is a multiple of 64.
void check_layout(seqan::hibf::layout::layout const & hibf_layout, size_t const number_of_user_bins, size_t const tmax)
{
    ASSERT_EQ(hibf_layout.user_bins.size(), number_of_user_bins);

    size_t const max_technical_bins{seqan::hibf::next_multiple_of_64(tmax)};

    constexpr size_t merged_bin{static_cast<size_t>(-1)};
    constexpr size_t empty_bin{static_cast<size_t>(-2)};
    std::map<std::vector<size_t>, std::vector<size_t>> ibfs{}; // path -> content of each technical bin

    auto get_ibf = [&](std::vector<size_t> const & path) -> std::vector<size_t> &
    {
        return ibfs.try_emplace(path, max_technical_bins, empty_bin).first->second;
    };

    for (size_t ub = 0; ub < number_of_user_bins; ++ub)
    {
        auto const & user_bin = hibf_layout.user_bins[ub];
        ASSERT_EQ(user_bin.idx, ub);

        std::vector<size_t> path{};
        for (size_t const tb : user_bin.previous_TB_indices)
        {
            ASSERT_LT(tb, max_technical_bins);
            auto & bin = get_ibf(path)[tb];
            ASSERT_TRUE(bin == empty_bin || bin == merged_bin) << "user bin " << ub << " passes through a filled bin";
            bin = merged_bin;
            path.push_back(tb);
        }

        ASSERT_GE(user_bin.number_of_technical_bins, 1u);
        ASSERT_LE(user_bin.storage_TB_id + user_bin.number_of_technical_bins, max_technical_bins);
        auto & ibf = get_ibf(path);
        for (size_t tb = user_bin.storage_TB_id; tb < user_bin.storage_TB_id + user_bin.number_of_technical_bins; ++tb)
        {
            ASSERT_EQ(ibf[tb], empty_bin) << "user bin " << ub << " is stored in an occupied bin";
            ibf[tb] = ub;
        }
    }

    EXPECT_LT(hibf_layout.top_level_max_bin_id, tmax); // the top level is laid out by the fast layout

    std::set<std::vector<size_t>> lower_level_ibfs{};
    for (auto const & [path, bins] : ibfs)
        if (!path.empty())
            lower_level_ibfs.insert(path);

    std::set<std::vector<size_t>> max_bin_ibfs{};
    for (auto const & max_bin : hibf_layout.max_bins)
    {
        EXPECT_LT(max_bin.id, max_technical_bins);
        EXPECT_TRUE(max_bin_ibfs.insert(max_bin.previous_TB_indices).second) << "duplicate max bin entry";
    }

    EXPECT_EQ(max_bin_ibfs, lower_level_ibfs);
}

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

    std::vector<size_t> positions(number_of_user_bins);
    std::iota(positions.begin(), positions.end(), 0u);

    seqan::hibf::layout::layout hibf_layout{};
    chopper::layout::fast_layout(config, positions, cardinalities, sketches, minHash_sketches, hibf_layout);

    check_layout(hibf_layout, number_of_user_bins, tmax);

    // Recursive layouts produce at least three levels.
    size_t max_depth{0};
    for (auto const & user_bin : hibf_layout.user_bins)
        max_depth = std::max(max_depth, user_bin.previous_TB_indices.size());
    EXPECT_GE(max_depth, 2u);
}
