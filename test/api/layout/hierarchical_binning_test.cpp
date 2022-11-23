#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <robin_hood.h>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/hierarchical_binning.hpp>

#include "../api_test.hpp"

TEST(hierarchical_binning_test, filenames_and_kmer_counts_size_differs)
{
    chopper::configuration config;
    config.tmax = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    data.filenames = {"seq0", "seq1"};   // 2 filenames
    data.kmer_counts = {500, 1000, 500}; // 3 kmer_counts :(

    EXPECT_THROW((chopper::layout::hierarchical_binning{data, config}), std::runtime_error);
}

TEST(hierarchical_binning_test, small_example)
{
    chopper::configuration config;
    config.tmax = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;
    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:3

    std::string expected_file
    {
        "#MERGED_BIN_1 max_bin_id:22\n"
        "#MERGED_BIN_2 max_bin_id:22\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq7\t0\t1\n"
        "seq4\t1;0\t1;22\n"
        "seq5\t1;22\t1;21\n"
        "seq6\t1;43\t1;21\n"
        "seq0\t2;0\t1;22\n"
        "seq2\t2;22\t1;21\n"
        "seq3\t2;43\t1;21\n"
        "seq1\t3\t1\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file) << output_buffer.str() << std::endl << expected_file << std::endl;
}

TEST(hierarchical_binning_test, another_example)
{
    chopper::configuration config;
    config.tmax = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {50, 1000, 1000, 50, 5, 10, 10, 5};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::string expected_file
    {
        "#MERGED_BIN_0;0 max_bin_id:56\n"
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq6\t0;0;0\t1;1;31\n"
        "seq5\t0;0;31\t1;1;25\n"
        "seq7\t0;0;56\t1;1;4\n"
        "seq4\t0;0;60\t1;1;4\n"
        "seq0\t0;1\t1;2\n"
        "seq3\t0;3\t1;2\n"
        "seq2\t1\t2\n"
        "seq1\t3\t2\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, high_level_max_bin_id_is_0)
{
    chopper::configuration config;
    config.tmax = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3"};
    data.kmer_counts = {500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::string expected_file
    {
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq3\t0\t1\n"
        "seq2\t1\t1\n"
        "seq1\t2\t1\n"
        "seq0\t3\t1\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, knuts_example)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4"};
    data.kmer_counts = {60, 600, 1000, 800, 800};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u);

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:63\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq1\t0;0\t1;63\n"
        "seq0\t0;63\t1;1\n"
        "seq4\t1\t1\n"
        "seq3\t2\t1\n"
        "seq2\t3\t2\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, four_level_hibf)
{
    chopper::configuration config;
    config.tmax = 2;
    // config.debug = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5"};
    data.kmer_counts = {11090, 5080, 3040, 1020, 510, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::string expected_file
    {
        "#MERGED_BIN_0;0;0;0 max_bin_id:33\n"
        "#MERGED_BIN_0;0;0 max_bin_id:1\n"
        "#MERGED_BIN_0;0 max_bin_id:1\n"
        "#MERGED_BIN_0 max_bin_id:1\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq4\t0;0;0;0;0\t1;1;1;1;33\n"
        "seq5\t0;0;0;0;33\t1;1;1;1;31\n"
        "seq3\t0;0;0;1\t1;1;1;1\n"
        "seq2\t0;0;1\t1;1;1\n"
        "seq1\t0;1\t1;1\n"
        "seq0\t1\t1\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 2;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3"};
    data.kmer_counts = {500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq2\t0;0\t1;32\n"
        "seq3\t0;32\t1;32\n"
        "seq0\t1;0\t1;32\n"
        "seq1\t1;32\t1;32\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file) << output_buffer.str();
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin_with_debug)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 2;
    config.debug = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3"};
    data.kmer_counts = {500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX\n"
        "seq2\t0;0\t1;32\t1000;16\t1000;140\t1.00;9.02\t2;64\n"
        "seq3\t0;32\t1;32\t1000;16\t1000;140\t1.00;9.02\t2;64\n"
        "seq0\t1;0\t1;32\t1000;16\t1000;140\t1.00;9.02\t2;64\n"
        "seq1\t1;32\t1;32\t1000;16\t1000;140\t1.00;9.02\t2;64\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file) << output_buffer.str();
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin_and_leads_to_recursive_call)
{
    chopper::configuration config;
    config.alpha = 1;
    config.tmax = 2;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    chopper::layout::data_store data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"};
    data.kmer_counts = {500, 500, 500, 500, 500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2, config.tmax);
    chopper::layout::hibf_statistics global_stats_dummy{};
    data.stats = &global_stats_dummy.top_level_ibf;

    chopper::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::string expected_file
    {
        "#MERGED_BIN_0;0 max_bin_id:0\n"
        "#MERGED_BIN_0;1 max_bin_id:0\n"
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#MERGED_BIN_1;0 max_bin_id:0\n"
        "#MERGED_BIN_1;1 max_bin_id:0\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq5\t0;0;0\t1;1;32\n"
        "seq4\t0;0;32\t1;1;32\n"
        "seq7\t0;1;0\t1;1;32\n"
        "seq6\t0;1;32\t1;1;32\n"
        "seq1\t1;0;0\t1;1;32\n"
        "seq0\t1;0;32\t1;1;32\n"
        "seq3\t1;1;0\t1;1;32\n"
        "seq2\t1;1;32\t1;1;32\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file) << output_buffer.str();
}
