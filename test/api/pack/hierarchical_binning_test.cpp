#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/pack/hierarchical_binning.hpp>

TEST(hierarchical_binning_test, small_example)
{
    std::vector<std::string> seq_ids{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    std::vector<size_t> counts{500, 1000, 500, 500, 500, 500, 500, 500};
    seqan3::test::tmp_filename output_filename{"output.binning"};

    pack_config config;
    config.bins = 4;
    config.output_filename = output_filename.get_path();

    hierarchical_binning algo{seq_ids, counts, config};
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

    // high level ibf file:
    std::ifstream output_file{output_filename.get_path()};
    std::string output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file);
}
