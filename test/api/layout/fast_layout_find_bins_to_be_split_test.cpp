// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EXIT, TEST

#include <cstdlib>  // for exit
#include <numeric>  // for iota
#include <unistd.h> // for alarm
#include <vector>   // for vector

#include <chopper/layout/fast_layout_find_bins_to_be_split.hpp>

TEST(find_bins_to_be_split_test, dissertation_mehringer_tiny_example)
{
    std::vector<size_t> const cardinalities{120, 60, 30, 30};
    std::vector<size_t> const sorted_positions{0, 1, 2, 3};

    // t_s = number of split techincal bins
    auto const [idx, t_s] = chopper::layout::find_bins_to_be_split(sorted_positions, cardinalities, /*t=*/45, 3);
    EXPECT_EQ(t_s, 2); // two split bins
    EXPECT_EQ(idx, 1); // points to the end of the range of indices to split -> only user bin 0 (cardinality 120)
}

TEST(find_bins_to_be_split_test, all_bins_split)
{
    std::vector<size_t> const cardinalities(60, 1100); // 60 times a value of 1100
    std::vector<size_t> sorted_positions(cardinalities.size());
    std::iota(sorted_positions.begin(), sorted_positions.end(), 0);

    // t_s = number of split techincal bins
    auto const [idx, t_s] = chopper::layout::find_bins_to_be_split(sorted_positions, cardinalities, /*t=*/1000, 63);
    EXPECT_EQ(t_s, 63); // 63 split bins
    EXPECT_EQ(idx, 60); // all user bins shall be split
}

TEST(find_bins_to_be_split_test, small_number_edge_case)
{
    // small number edge case: the `std::max<size_t>(threshold + 1, ... )` catches here
    // but increasing by 1 results in none of the bins split. They will be "merged" then
    // which is not optimal but fine. The merging algorithm will distribute each user bin
    // into one technical bin and 4 technical bins will be empty. With these small numbers
    // that is fine.
    std::vector<size_t> const cardinalities(60, 11); // 60 times a value of 11
    std::vector<size_t> sorted_positions(cardinalities.size());
    std::iota(sorted_positions.begin(), sorted_positions.end(), 0);

    // t_s = number of split techincal bins
    auto const [idx, t_s] = chopper::layout::find_bins_to_be_split(sorted_positions, cardinalities, /*t=*/10, 63);
    EXPECT_EQ(t_s, 0); // 63 split bins
    EXPECT_EQ(idx, 0); // no user bins shall be split
}

TEST(find_bins_to_be_split_test, threshold_zero_death_in_debug)
{
#ifdef NDEBUG
    GTEST_SKIP() << "Debug-only test";
#endif
    std::vector<size_t> const cardinalities{10, 20, 30};
    std::vector<size_t> const sorted_positions{0, 1, 2};

    EXPECT_DEATH(chopper::layout::find_bins_to_be_split(sorted_positions, cardinalities, /*t=*/0, 2), "threshold > 0");
}
