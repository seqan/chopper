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

TEST_F(cli_test, chopper_layout)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{tmp_dir.path() / "data.filenames"};
    std::filesystem::path const layout_filename{tmp_dir.path() / "output.binning"};
    std::filesystem::path const sketches_dir{tmp_dir.path() / "sketches"};

    {
        std::ofstream fout{input_filename};
        fout << data("seq1.fa").string() << '\n'
             << data("seq2.fa").string() << '\n'
             << data("seq3.fa").string() << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--kmer",
                                         "15",
                                         "--input",
                                         input_filename.c_str(),
                                         "--tmax",
                                         "64",
                                         "--output-sketches-to",
                                         sketches_dir.c_str(),
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    // Layout file should exist. Content is checked in other test.
    EXPECT_TRUE(std::filesystem::exists(layout_filename));

    EXPECT_TRUE(std::filesystem::exists(sketches_dir)); // directory exists
    EXPECT_TRUE(std::filesystem::exists(sketches_dir / "seq1.hll"));
    EXPECT_TRUE(std::filesystem::exists(sketches_dir / "seq2.hll"));
    EXPECT_TRUE(std::filesystem::exists(sketches_dir / "seq3.hll"));
}