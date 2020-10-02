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
