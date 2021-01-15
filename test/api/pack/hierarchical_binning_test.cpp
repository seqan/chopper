#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/pack/hierarchical_binning.hpp>

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
    algo.execute();

    std::string expected_file
    {
        "#MERGED_BIN_2 max_bin_id:16\n"
        "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2\n"
        "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "SPLIT_BIN_0\tseq7\t1\t500\n"
        "SPLIT_BIN_1\tseq6\t1\t500\n"
        "MERGED_BIN_2_0\tseq0\t16\t32\n"
        "MERGED_BIN_2_16\tseq2\t12\t42\n"
        "MERGED_BIN_2_28\tseq3\t12\t42\n"
        "MERGED_BIN_2_40\tseq4\t12\t42\n"
        "MERGED_BIN_2_52\tseq5\t12\t42\n"
        "SPLIT_BIN_3\tseq1\t1\t1000\n"
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
    algo.execute();

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:35\n"
        "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_1\n"
        "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "MERGED_BIN_0_0\tseq0\t35\t2\n"
        "MERGED_BIN_0_35\tseq3\t17\t3\n"
        "MERGED_BIN_0_52\tseq5\t4\t3\n"
        "MERGED_BIN_0_56\tseq6\t4\t3\n"
        "MERGED_BIN_0_60\tseq4\t2\t3\n"
        "MERGED_BIN_0_62\tseq7\t2\t3\n"
        "SPLIT_BIN_1\tseq2\t2\t500\n"
        "SPLIT_BIN_3\tseq1\t2\t500\n"
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
    algo.execute();

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_1\n"
        "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "MERGED_BIN_0_0\tseq1\t58\t11\n"
        "MERGED_BIN_0_58\tseq0\t6\t10\n"
        "SPLIT_BIN_1\tseq4\t1\t800\n"
        "SPLIT_BIN_2\tseq3\t1\t800\n"
        "SPLIT_BIN_3\tseq2\t2\t500\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}
