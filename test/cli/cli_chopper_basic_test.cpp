// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string> // strings

#include <seqan3/test/tmp_directory.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("chopper");
    std::string expected{"chopper - Compute an HIBF layout\n"
                         "================================\n"
                         "    chopper --input <file> [--output <file>] [--threads <number>] [--kmer\n    <number>] "
                         "[--fpr <number>] [--hash <number>] [--disable-estimate-union]\n    "
                         "[--disable-rearrangement]\n    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, chopper_cmd_error_unknown_option)
{
    cli_test_result result = execute_app("chopper", "--unkown-option");
    std::string expected{"[ERROR] Option --input is required but not set.\n"};
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_empty_file)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path empty_file{tmp_dir.path() / "empty.count"};

    {
        std::ofstream ofs{empty_file.string()}; // opens file, s.t. it exists but is empty
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax",
                                         "64", /* required option */
                                         "--input",
                                         empty_file.c_str());

    std::string expected{"[ERROR] The file " + empty_file.string() + " appears to be empty.\n"};
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_non_existing_path_in_input)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path non_existing_content{tmp_dir.path() / "file_with_one_path_that_does_not.exist"};

    {
        std::ofstream ofs{non_existing_content.string()}; // opens file, s.t. it exists but is empty
        ofs << "/I/do/no/exist.fa\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax",
                                         "64", /* required option */
                                         "--input",
                                         non_existing_content.c_str());

    std::string expected{"[ERROR] File /I/do/no/exist.fa does not exist!\n"};
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
