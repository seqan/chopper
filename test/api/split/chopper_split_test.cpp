#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/split/split_data.hpp>
#include <chopper/split/sequence_input.hpp>
#include <chopper/split/chopper_split.hpp>

TEST(chopper_split_test, simple_example)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename output_filename{"small_traverse.out"};
    const char * argv[] = {"./chopper-split", "-k", "15", "-w", "25", "-b", "3",
                           "-s", input_filename.c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 11;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, false};

    chopper_split(split_parser);

    // compare results
    std::ifstream expected_file{DATADIR"small_traverse.out"};
    std::ifstream output_file{output_filename.get_path()};

    std::string expected_line;
    std::string output_line;

    while (std::getline(expected_file, expected_line) && std::getline(output_file, output_line))
        EXPECT_EQ(expected_line, output_line);

    EXPECT_FALSE(std::getline(expected_file, expected_line)); // both files are exhausted
    EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted
}

TEST(chopper_split_test, no_s_or_f_option)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename output_filename{"small_traverse.out"};
    const char * argv[] = {"./chopper-split", "-k", "15"};
    int argc = 3;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, false};

    EXPECT_THROW(chopper_split(split_parser), std::runtime_error);
}

TEST(chopper_split_test, high_level_ibf)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    {
        std::ofstream fout{data_filename.get_path()};
        fout << "FILE_OR_COLOR_ID\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
             << input_filename + "\t2\t500\n"
             << input_filename + "\t2\t500\n"
             << "COLORFUL_MERGED_BIN_0\t1\t2500\n"
             << input_filename + "\t3\t1000\n";
    }

    seqan3::test::tmp_filename output_filename{"traverse"};

    const char * argv[] = {"./chopper-split", "-k", "15", "-w", "25",
                           "-f", data_filename.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 9;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, false};

    EXPECT_EQ(chopper_split(split_parser), 0);

    std::vector<std::string> const expected_output
    {
        {
            "[(0,209),(0,289),(0,209)]\n"
            "[(209,400),(289,480),(209,481)]\n"
        },
        {
            "[(0,209),(0,289),(0,209)]\n"
            "[(209,400),(289,480),(209,481)]\n"
        },
        {
            "[(0,163),(0,186),(0,163)]\n"
            "[(163,247),(186,327),(163,284)]\n"
            "[(247,400),(327,480),(284,481)]\n"
        }
    };

    // compare results
    for (size_t batch_number = 0; batch_number < 3; ++batch_number)
    {
        std::ifstream output_file{output_filename.get_path().string() + "_" + std::to_string(batch_number) + ".out"};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        EXPECT_EQ(output_file_str, expected_output[batch_number]);
    }
}
