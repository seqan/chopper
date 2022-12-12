#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/simple_binning.hpp>

TEST(simple_binning_test, small_example)
{
    std::stringstream output_buffer;
    chopper::data_store data{.output_buffer = &output_buffer, .header_buffer = &output_buffer};
    data.kmer_counts = {100, 40, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};
    data.fp_correction = std::vector<double>(65, 1.0);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::simple_binning algo{data, 9};
    size_t max_bin = algo.execute();

    std::string expected{"seq4\t0\t1\n"
                         "seq3\t1\t1\n"
                         "seq2\t2\t2\n"
                         "seq1\t4\t5\n"};

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, uniform_distribution)
{
    std::stringstream output_buffer;
    chopper::data_store data{.output_buffer = &output_buffer, .header_buffer = &output_buffer};
    data.kmer_counts = {20, 20, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};
    data.fp_correction = std::vector<double>(65, 1.0);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::simple_binning algo{data, 4u};
    size_t max_bin = algo.execute();

    std::string expected{"seq4\t0\t1\n"
                         "seq3\t1\t1\n"
                         "seq2\t2\t1\n"
                         "seq1\t3\t1\n"};

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    std::stringstream output_buffer;
    chopper::layout::hibf_statistics global_stats_dummy{};

    chopper::data_store data{.output_buffer = &output_buffer,
                             .header_buffer = &output_buffer,
                             .stats = &global_stats_dummy.top_level_ibf,
                             .filenames = {"seq1", "seq2", "seq3", "seq4"},
                             .kmer_counts = {100, 40, 20, 20},
                             .fp_correction = std::vector<double>(65, 1.0)};

    EXPECT_THROW((chopper::layout::simple_binning{data, 2}), std::runtime_error);
}
