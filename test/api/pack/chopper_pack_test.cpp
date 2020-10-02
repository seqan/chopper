#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/pack/chopper_pack.hpp>

TEST(chopper_pack_test, small_example)
{
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    {
        std::ofstream fout{data_filename.get_path()};
        fout << "seq0\t500\n"
             << "seq1\t1000\n"
             << "seq2\t500\n"
             << "seq3\t500\n"
             << "seq4\t500\n"
             << "seq5\t500\n"
             << "seq6\t500\n"
             << "seq7\t500\n";
    }

    const char * argv[] = {"./chopper-pack", "-b", "4", "-f", data_filename.get_path().c_str()};
    int argc = 5;
    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, false};

    chopper_pack(pack_parser);

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

// std::string input_filename = DATADIR"filenames_counts_and_extra_information.tsv";
