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

TEST_F(cli_test, timing_output)
{
    std::string const seq_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{tmp_dir.path() / "data.tsv"};
    std::filesystem::path const layout_filename{tmp_dir.path() / "output.layout"};
    std::filesystem::path const timing_filename{tmp_dir.path() / "output.timings"};

    {
        std::ofstream fout{input_filename};
        fout << seq_filename << '\n' << seq_filename << '\n' << seq_filename << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--output",
                                         layout_filename.c_str(),
                                         "--timing-output",
                                         timing_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    EXPECT_TRUE(std::filesystem::exists(timing_filename)); // file should have been written
    // not not check output since it is not relevant how exectly the timings look like
}
