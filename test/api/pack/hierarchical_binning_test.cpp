#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/hierarchical_binning.hpp>
#include <robin_hood.h>

#include "../api_test.hpp"

TEST(hierarchical_binning_test, small_example)
{
    pack_config config;
    config.t_max = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500};
    data.compute_fp_correction(0.05, 2);
    hierarchical_binning algo{data, config};
    EXPECT_EQ(std::get<0>(algo.execute()), 3); // #HIGH_LEVEL_IBF max_bin_id:3

    std::string expected_file
    {
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#MERGED_BIN_2 max_bin_id:0\n"
        "#MERGED_BIN_3 max_bin_id:52\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq7\t0\t1\n"
        "seq5\t1;0\t1;32\n"
        "seq6\t1;32\t1;32\n"
        "seq3\t2;0\t1;32\n"
        "seq4\t2;32\t1;32\n"
        "seq1\t3;0\t1;52\n"
        "seq0\t3;52\t1;6\n"
        "seq2\t3;58\t1;6\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, another_example)
{
    pack_config config;
    config.t_max = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {50, 1000, 1000, 50, 5, 10, 10, 5};
    data.compute_fp_correction(0.05, 2);

    hierarchical_binning algo{data, config};
    EXPECT_EQ(std::get<0>(algo.execute()), 1); // #HIGH_LEVEL_IBF max_bin_id:1

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

TEST(hierarchical_binning_test, knuts_example)
{
    pack_config config;
    config.alpha = 1;
    config.t_max = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4"};
    data.kmer_counts = {60, 600, 1000, 800, 800};
    data.compute_fp_correction(0.05, 2);

    hierarchical_binning algo{data, config};
    EXPECT_EQ(std::get<0>(algo.execute()), 1);

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
