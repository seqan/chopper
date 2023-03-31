#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/simple_binning.hpp>

TEST(simple_binning_test, small_example)
{
    chopper::layout::layout hibf_layout;
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq1", "seq2", "seq3", "seq4"},
                             .kmer_counts = {100, 40, 20, 20}};

    data.fp_correction = std::vector<double>(65, 1.0);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::simple_binning algo{data, 9};
    size_t max_bin = algo.execute();

    std::vector<chopper::layout::layout::user_bin> expected{{"seq4", {}, 1, 0},
                                                            {"seq3", {}, 1, 1},
                                                            {"seq2", {}, 2, 2},
                                                            {"seq1", {}, 5, 4}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, uniform_distribution)
{
    chopper::layout::layout hibf_layout;
    chopper::data_store data{.hibf_layout = &hibf_layout,
                             .filenames = {"seq1", "seq2", "seq3", "seq4"},
                             .kmer_counts = {20, 20, 20, 20}};

    data.fp_correction = std::vector<double>(65, 1.0);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::simple_binning algo{data, 4u};
    size_t max_bin = algo.execute();

    std::vector<chopper::layout::layout::user_bin> expected{{"seq4", {}, 1, 0},
                                                            {"seq3", {}, 1, 1},
                                                            {"seq2", {}, 1, 2},
                                                            {"seq1", {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    chopper::layout::layout hibf_layout;
    chopper::layout::hibf_statistics global_stats_dummy{};

    chopper::data_store data{.stats = &global_stats_dummy.top_level_ibf,
                             .hibf_layout = &hibf_layout,
                             .filenames = {"seq1", "seq2", "seq3", "seq4"},
                             .kmer_counts = {100, 40, 20, 20},
                             .fp_correction = std::vector<double>(65, 1.0)};

    EXPECT_THROW((chopper::layout::simple_binning{data, 2}), std::runtime_error);
}
