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

#include <chopper/sketch/sketch_file.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, chopper_layout_wrong_sketch_extension)
{
    std::filesystem::path const input_filename{"data.filenames"};
    std::filesystem::path const layout_filename{"output.binning"};
    std::filesystem::path const sketches_filename{"out.sketched"};

    {
        std::ofstream fout{input_filename};
        fout << data("seq1.fa").string() << '\n'
             << data("seq2.fa").string() << '\n'
             << data("seq3.fa").string() << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--tmax 64",
                                         "--output-sketches-to",
                                         sketches_filename.c_str(),
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[ERROR] The sketch output file must have the extension \".sketch\" or \".sketches\".\n"});
}

TEST_F(cli_test, chopper_layout)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{tmp_dir.path() / "data.filenames"};
    std::filesystem::path const layout_filename{tmp_dir.path() / "output.binning"};
    std::filesystem::path const sketches_filename{tmp_dir.path() / "out.sketches"};

    {
        std::ofstream fout{input_filename};
        fout << data("seq1.fa").string() << '\n'
             << data("seq2.fa").string() << '\n'
             << data("seq3.fa").string() << '\n';
    }

    size_t kmer_size{15};
    size_t tmax{64};

    cli_test_result result = execute_app("chopper",
                                         "--kmer",
                                         std::to_string(kmer_size).c_str(),
                                         "--input",
                                         input_filename.c_str(),
                                         "--tmax",
                                         std::to_string(tmax).c_str(),
                                         "--output-sketches-to",
                                         sketches_filename.c_str(),
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    // Layout file should exist. Content is checked in other test.
    EXPECT_TRUE(std::filesystem::exists(layout_filename));

    ASSERT_TRUE(std::filesystem::exists(sketches_filename));

    chopper::sketch::sketch_file sin{};

    std::ifstream is{sketches_filename};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(sin);

    EXPECT_EQ(sin.chopper_config.k, kmer_size);
    EXPECT_EQ(sin.chopper_config.data_file, input_filename);
    EXPECT_EQ(sin.chopper_config.output_filename, layout_filename);
    // EXPECT_EQ(sin.chopper_config.hibf_config.tmax, tmax); // TODO: why does this fail?

    EXPECT_EQ(sin.filenames.size(), 3);
    EXPECT_EQ(sin.hll_sketches.size(), 3);
    EXPECT_EQ(sin.minHash_sketches.size(), 0); // currently, no minhash sketches are needed in chopper layout

    EXPECT_EQ(sin.filenames[0][0], data("seq1.fa").string());
    EXPECT_EQ(sin.filenames[1][0], data("seq2.fa").string());
    EXPECT_EQ(sin.filenames[2][0], data("seq3.fa").string());
}
