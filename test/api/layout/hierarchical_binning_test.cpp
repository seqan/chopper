#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <robin_hood.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/layout/compute_fp_correction.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/hierarchical_binning.hpp>

TEST(hierarchical_binning_test, filenames_and_kmer_counts_size_differs)
{
    chopper::configuration config;
    config.tmax = 4;

    chopper::layout::layout hibf_layout{};
    chopper::layout::hibf_statistics global_stats_dummy{};

    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1"},    // 2 filenames
                             .kmer_counts = {500, 1000, 500}}; // 3 kmer_counts
    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);

    EXPECT_THROW((chopper::layout::hierarchical_binning{data, config}), std::runtime_error);
}

TEST(hierarchical_binning_test, small_example)
{
    chopper::configuration config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                             .kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;
    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:3

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{1}, 22}, {{2}, 22}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq7", {}, 1, 0},
                                                                      {"seq4", {1}, 22, 0},
                                                                      {"seq5", {1}, 21, 22},
                                                                      {"seq6", {1}, 21, 43},
                                                                      {"seq0", {2}, 22, 0},
                                                                      {"seq2", {2}, 21, 22},
                                                                      {"seq3", {2}, 21, 43},
                                                                      {"seq1", {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, another_example)
{
    chopper::configuration config;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                             .kmer_counts = {50, 1000, 1000, 50, 5, 10, 10, 5}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0}, 56}, {{0}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq6", {0, 0}, 42, 0},
                                                                      {"seq5", {0, 0}, 14, 42},
                                                                      {"seq7", {0, 0}, 4, 56},
                                                                      {"seq4", {0, 0}, 4, 60},
                                                                      {"seq0", {0}, 2, 1},
                                                                      {"seq3", {0}, 2, 3},
                                                                      {"seq2", {}, 2, 1},
                                                                      {"seq1", {}, 2, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, high_level_max_bin_id_is_0)
{
    chopper::configuration config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3"},
                             .kmer_counts = {500, 500, 500, 500}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq3", {}, 1, 0},
                                                                      {"seq2", {}, 1, 1},
                                                                      {"seq1", {}, 1, 2},
                                                                      {"seq0", {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, knuts_example)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::data_store data{.hibf_layout = &hibf_layout};
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4"};
    data.kmer_counts = {60, 600, 1000, 800, 800};
    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0}, 63}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq1", {0}, 63, 0},
                                                                      {"seq0", {0}, 1, 63},
                                                                      {"seq4", {}, 1, 1},
                                                                      {"seq3", {}, 1, 2},
                                                                      {"seq2", {}, 2, 3}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, four_level_hibf)
{
    chopper::configuration config;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::layout hibf_layout{};
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5"},
                             .kmer_counts = {11090, 5080, 3040, 1020, 510, 500}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0, 0, 0}, 33},
                                                                    {{0, 0, 0}, 1},
                                                                    {{0, 0}, 1},
                                                                    {{0}, 1}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq4", {0, 0, 0, 0}, 33, 0},
                                                                      {"seq5", {0, 0, 0, 0}, 31, 33},
                                                                      {"seq3", {0, 0, 0}, 1, 1},
                                                                      {"seq2", {0, 0}, 1, 1},
                                                                      {"seq1", {0}, 1, 1},
                                                                      {"seq0", {}, 1, 1}};

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
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3"},
                             .kmer_counts = {500, 500, 500, 500}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0}, 0}, {{1}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq2", {0}, 32, 0},
                                                                      {"seq3", {0}, 32, 32},
                                                                      {"seq0", {1}, 32, 0},
                                                                      {"seq1", {1}, 32, 32}};

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
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                             .kmer_counts = {500, 500, 500, 500, 500, 500, 500, 500}};

    data.fp_correction = chopper::layout::compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<chopper::layout::layout::max_bin> expected_max_bins{{{0, 0}, 0},
                                                                    {{0, 1}, 0},
                                                                    {{0}, 0},
                                                                    {{1, 0}, 0},
                                                                    {{1, 1}, 0},
                                                                    {{1}, 0}};

    std::vector<chopper::layout::layout::user_bin> expected_user_bins{{"seq5", {0, 0}, 32, 0},
                                                                      {"seq4", {0, 0}, 32, 32},
                                                                      {"seq7", {0, 1}, 32, 0},
                                                                      {"seq6", {0, 1}, 32, 32},
                                                                      {"seq1", {1, 0}, 32, 0},
                                                                      {"seq0", {1, 0}, 32, 32},
                                                                      {"seq3", {1, 1}, 32, 0},
                                                                      {"seq2", {1, 1}, 32, 32}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}
