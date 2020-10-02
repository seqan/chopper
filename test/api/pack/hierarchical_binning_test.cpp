#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/hierarchical_binning.hpp>

TEST(hierarchical_binning_test, small_example)
{
    std::vector<std::string> seq_ids{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    std::vector<size_t> counts{500, 1000, 500, 500, 500, 500, 500, 500};

    pack_config config;
    config.bins = 4;

    hierarchical_binning algo{seq_ids, counts, config};
    algo.dp_algorithm();

    std::istringstream expected_file{std::string
    {
        "BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "SPLIT_BIN_0\tseq7\t1\t500\n"
        "SPLIT_BIN_1\tseq6\t1\t500\n"
        "COLORFUL_MERGED_BIN_2_0\tseq0\t16\t32\n"
        "COLORFUL_MERGED_BIN_2_1\tseq2\t12\t42\n"
        "COLORFUL_MERGED_BIN_2_2\tseq3\t12\t42\n"
        "COLORFUL_MERGED_BIN_2_3\tseq4\t12\t42\n"
        "COLORFUL_MERGED_BIN_2_4\tseq5\t12\t42\n"
        "SPLIT_BIN_3\tseq1\t1\t1000\n"
    }};

    std::string expected_line;
    std::string output_line;

    // high level ibf file:
    {
        std::ifstream output_file{"output.binning"};
        while (std::getline(expected_file, expected_line) && std::getline(output_file, output_line))
            EXPECT_EQ(expected_line, output_line);

        EXPECT_FALSE(std::getline(expected_file, expected_line)); // both files are exhausted
        EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted
    }
}
