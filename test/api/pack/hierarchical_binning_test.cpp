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

    hierarchical_binning algo{seq_ids, counts, 4};
    algo.dp_algorithm();

    std::istringstream expected_high_level{std::string
    {
        "FILE_OR_COLOR_ID\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "seq7\t1\t500\n"
        "seq6\t1\t500\n"
        "COLORFUL_MERGED_BIN_0\t1\t2500\n"
        "seq1\t1\t1000\n"
    }};

    std::istringstream expected_low_level{std::string
    {
        "COLOR_ID\tFILE_ID\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "COLORFUL_MERGED_BIN_0\tseq0\t16\t32\n"
        "COLORFUL_MERGED_BIN_0\tseq2\t12\t42\n"
        "COLORFUL_MERGED_BIN_0\tseq3\t12\t42\n"
        "COLORFUL_MERGED_BIN_0\tseq4\t12\t42\n"
        "COLORFUL_MERGED_BIN_0\tseq5\t12\t42\n"
    }};

    std::string expected_line;
    std::string output_line;

    // high level ibf file:
    {
        std::ifstream output_file{"high_level_ibf.binning"};
        while (std::getline(expected_high_level, expected_line) && std::getline(output_file, output_line))
            EXPECT_EQ(expected_line, output_line);

        EXPECT_FALSE(std::getline(expected_high_level, expected_line)); // both files are exhausted
        EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted
    }

    // high level ibf file:
    {
        std::ifstream output_file{"low_level_ibfs.binning"};
        while (std::getline(expected_low_level, expected_line) && std::getline(output_file, output_line))
            EXPECT_EQ(expected_line, output_line);

        EXPECT_FALSE(std::getline(expected_low_level, expected_line)); // both files are exhausted
        EXPECT_FALSE(std::getline(output_file, output_line)) << output_line << std::endl; // both files are exhausted
    }
}
