#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <robin_hood.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/layout/compute_fp_correction.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/hierarchical_binning.hpp>

TEST(hierarchical_binning_test, small_example)
{
    chopper::configuration config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{500, 1000, 500, 500, 500, 500, 500, 500};

    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:3

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{1}, 22}, {{2}, 22}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{7, {}, 1, 0},
                                                                      {4, {1}, 22, 0},
                                                                      {5, {1}, 21, 22},
                                                                      {6, {1}, 21, 43},
                                                                      {0, {2}, 22, 0},
                                                                      {2, {2}, 21, 22},
                                                                      {3, {2}, 21, 43},
                                                                      {1, {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, another_example)
{
    chopper::configuration config;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{50, 1000, 1000, 50, 5, 10, 10, 5};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0}, 56}, {{0}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{6, {0, 0}, 42, 0},
                                                                      {5, {0, 0}, 14, 42},
                                                                      {7, {0, 0}, 4, 56},
                                                                      {4, {0, 0}, 4, 60},
                                                                      {0, {0}, 2, 1},
                                                                      {3, {0}, 2, 3},
                                                                      {2, {}, 2, 1},
                                                                      {1, {}, 2, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, high_level_max_bin_id_is_0)
{
    chopper::configuration config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{500, 500, 500, 500};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{3, {}, 1, 0},
                                                                      {2, {}, 1, 1},
                                                                      {1, {}, 1, 2},
                                                                      {0, {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, knuts_example)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{60, 600, 1000, 800, 800};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0}, 63}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{1, {0}, 63, 0},
                                                                      {0, {0}, 1, 63},
                                                                      {4, {}, 1, 1},
                                                                      {3, {}, 1, 2},
                                                                      {2, {}, 2, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, four_level_hibf)
{
    chopper::configuration config;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{11090, 5080, 3040, 1020, 510, 500};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0, 0, 0}, 33},
                                                                    {{0, 0, 0}, 1},
                                                                    {{0, 0}, 1},
                                                                    {{0}, 1}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{4, {0, 0, 0, 0}, 33, 0},
                                                                      {5, {0, 0, 0, 0}, 31, 33},
                                                                      {3, {0, 0, 0}, 1, 1},
                                                                      {2, {0, 0}, 1, 1},
                                                                      {1, {0}, 1, 1},
                                                                      {0, {}, 1, 1}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{500, 500, 500, 500};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0}, 0}, {{1}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{2, {0}, 32, 0},
                                                                      {3, {0}, 32, 32},
                                                                      {0, {1}, 32, 0},
                                                                      {1, {1}, 32, 32}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin_and_leads_to_recursive_call)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{{}, {}, {}, {}};
    std::vector<size_t> kmer_counts{500, 500, 500, 500, 500, 500, 500, 500};
    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = kmer_counts};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0}, 0},
                                                                    {{0, 1}, 0},
                                                                    {{0}, 0},
                                                                    {{1, 0}, 0},
                                                                    {{1, 1}, 0},
                                                                    {{1}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{5, {0, 0}, 32, 0},
                                                                      {4, {0, 0}, 32, 32},
                                                                      {7, {0, 1}, 32, 0},
                                                                      {6, {0, 1}, 32, 32},
                                                                      {1, {1, 0}, 32, 0},
                                                                      {0, {1, 0}, 32, 32},
                                                                      {3, {1, 1}, 32, 0},
                                                                      {2, {1, 1}, 32, 32}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}
