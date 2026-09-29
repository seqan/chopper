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

TEST_F(cli_test, chopper_cmd_kmer_bigger_than_window)
{
    std::filesystem::path const input_filename{"bins.filenames"};
    {
        std::ofstream fout{input_filename};
    }

    cli_test_result result =
        execute_app("chopper", "--tmax 64", "--input", input_filename.c_str(), "--kmer 20", "--window 10");

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[ERROR] The k-mer size cannot be bigger than the window size.\n"});
}

TEST_F(cli_test, chopper_user_bin_with_few_kmers)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{tmp_dir.path() / "data.filenames"};
    std::filesystem::path const short_filename{tmp_dir.path() / "short.fa"};
    std::filesystem::path const layout_filename{tmp_dir.path() / "output.binning"};

    // 200 bases give 182 k-mers, too few to fill the MinHash sketches that only the fast layout needs.
    {
        std::ofstream fout{short_filename};
        fout << ">short\n";
        for (size_t i = 0; i < 200; ++i)
            fout << "ACGT"[(i * 7 + i / 3) % 4];
        fout << '\n';
    }

    {
        std::ofstream fout{input_filename};
        fout << data("seq1.fa").string() << '\n' << short_filename.string() << '\n';
    }

    cli_test_result result =
        execute_app("chopper", "--input", input_filename.c_str(), "--tmax", "64", "--output", layout_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
    EXPECT_TRUE(std::filesystem::exists(layout_filename));

    // The fast layout needs MinHash sketches and fails.
    result = execute_app("chopper",
                         "--fast-layout",
                         "--input",
                         input_filename.c_str(),
                         "--tmax",
                         "64",
                         "--output",
                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_TRUE(result.err.starts_with("[ERROR] Not enough kmers")) << result.err;
}
