#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/hierarchical_binning.hpp>

#include "../api_test.hpp"

TEST(hierarchical_binning_test, small_example)
{
    pack_config config;
    config.bins = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 2);

    std::string expected_file
    {
        "#MERGED_BIN_2;3 max_bin_id:0\n"
        "#MERGED_BIN_2 max_bin_id:3\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq7\t0\t1\t500\n"
        "seq6\t1\t1\t500\n"
        "seq0\t2;0\t1;1\t2500;500\n"
        "seq2\t2;1\t1;1\t2500;500\n"
        "seq3\t2;2\t1;1\t2500;500\n"
        "seq5\t2;3;0\t1;1;32\t2500;1000;16\n"
        "seq4\t2;3;32\t1;1;32\t2500;1000;16\n"
        "seq1\t3\t1\t1000\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, another_example)
{
    pack_config config;
    config.bins = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {50, 1000, 1000, 50, 5, 10, 10, 5};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1);

    std::string expected_file
    {
        "#MERGED_BIN_0;0 max_bin_id:0\n"
        "#MERGED_BIN_0 max_bin_id:3\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq7\t0;0;0\t1;1;58\t130;10;1\n"
        "seq4\t0;0;58\t1;1;6\t130;10;1\n"
        "seq5\t0;1\t1;1\t130;10\n"
        "seq6\t0;2\t1;1\t130;10\n"
        "seq0\t0;3\t1;1\t130;50\n"
        "seq3\t0;4\t1;1\t130;50\n"
        "seq2\t1\t2\t500\n"
        "seq1\t3\t2\t500\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, knuts_example)
{
    pack_config config;
    config.alpha = 1;
    config.bins = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4"};
    data.kmer_counts = {60, 600, 1000, 800, 800};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1);

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq1\t0;0\t1;58\t660;11\n"
        "seq0\t0;58\t1;6\t660;10\n"
        "seq4\t1\t1\t800\n"
        "seq3\t2\t1\t800\n"
        "seq2\t3\t2\t500\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}
