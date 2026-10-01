// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <seqan3/test/tmp_directory.hpp>

#include <chopper/layout/input.hpp>

#include "../api/api_test.hpp"
#include "cli_test.hpp"

namespace
{

std::filesystem::path write_input_file(std::filesystem::path const & directory)
{
    std::filesystem::path const input_filename{directory / "data.tsv"};
    std::ofstream fout{input_filename};
    fout << data("seq1.fa").string() << '\n'
         << data("seq2.fa").string() << '\n'
         << data("seq3.fa").string() << '\n'
         << data("small.fa").string() << '\n';
    return input_filename;
}

} // namespace

TEST_F(cli_test, chopper_layout_phibf)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const layout_filename{tmp_dir.path() / "phibf.layout"};

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--number-of-partitions",
                                         "2",
                                         "--partitioning-approach",
                                         "1",
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::ifstream layout_stream{layout_filename};
    auto [filenames, config, layouts] = chopper::layout::read_layouts_file(layout_stream);

    EXPECT_EQ(filenames.size(), 4u);
    EXPECT_EQ(config.number_of_partitions, 2u);
    EXPECT_EQ(config.partitioning_approach, 1);
    ASSERT_EQ(layouts.size(), 2u);

    std::vector<size_t> user_bins{};
    for (auto const & layout : layouts)
        for (auto const & user_bin : layout.user_bins)
            user_bins.push_back(user_bin.idx);
    std::ranges::sort(user_bins);
    EXPECT_EQ(user_bins, (std::vector<size_t>{0, 1, 2, 3}));
}

TEST_F(cli_test, chopper_layout_phibf_with_fast_layout)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const layout_filename{tmp_dir.path() / "phibf.layout"};

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--number-of-partitions",
                                         "2",
                                         "--fast-layout",
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{"[ERROR] You cannot combine --fast-layout with --number-of-partitions.\n"});
}

TEST_F(cli_test, chopper_layout_phibf_with_determine_best_tmax)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const layout_filename{tmp_dir.path() / "phibf.layout"};

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--number-of-partitions",
                                         "2",
                                         "--determine-best-tmax",
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[ERROR] You cannot combine --determine-best-tmax with --number-of-partitions.\n"});
}

TEST_F(cli_test, chopper_layout_phibf_from_sketch_file_without_minhashes)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const sketch_filename{tmp_dir.path() / "default.sketch"};

    // The default layout does not compute MinHash sketches.
    {
        std::filesystem::path const layout_filename{tmp_dir.path() / "default.layout"};
        cli_test_result result = execute_app("chopper",
                                             "--input",
                                             input_filename.c_str(),
                                             "--output-sketches-to",
                                             sketch_filename.c_str(),
                                             "--output",
                                             layout_filename.c_str());
        ASSERT_EQ(result.exit_code, 0);
    }

    // lsh_sim needs MinHash sketches.
    {
        std::filesystem::path const layout_filename{tmp_dir.path() / "lsh_sim.layout"};
        cli_test_result result = execute_app("chopper",
                                             "--input",
                                             sketch_filename.c_str(),
                                             "--number-of-partitions",
                                             "2",
                                             "--partitioning-approach",
                                             "6",
                                             "--output",
                                             layout_filename.c_str());

        EXPECT_NE(result.exit_code, 0);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err,
                  std::string{"[ERROR] The sketch file does not contain MinHash sketches, which the chosen "
                              "--partitioning-approach needs. Create the sketch file with --number-of-partitions.\n"});
    }

    // sorted does not need MinHash sketches.
    {
        std::filesystem::path const layout_filename{tmp_dir.path() / "sorted.layout"};
        cli_test_result result = execute_app("chopper",
                                             "--input",
                                             sketch_filename.c_str(),
                                             "--number-of-partitions",
                                             "2",
                                             "--partitioning-approach",
                                             "1",
                                             "--output",
                                             layout_filename.c_str());

        EXPECT_EQ(result.exit_code, 0);
        EXPECT_EQ(result.out, std::string{});
        EXPECT_EQ(result.err, std::string{});

        std::ifstream layout_stream{layout_filename};
        auto [filenames, config, layouts] = chopper::layout::read_layouts_file(layout_stream);
        EXPECT_EQ(layouts.size(), 2u);
    }
}

TEST_F(cli_test, chopper_layout_phibf_unknown_partitioning_approach)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const layout_filename{tmp_dir.path() / "phibf.layout"};

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--number-of-partitions",
                                         "2",
                                         "--partitioning-approach",
                                         "7",
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[ERROR] Validation failed for option --partitioning-approach: Value 7 is not in range "
                          "[0,6].\n"});
}

TEST_F(cli_test, chopper_layout_phibf_more_partitions_than_user_bins)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};
    std::filesystem::path const layout_filename{tmp_dir.path() / "phibf.layout"};

    cli_test_result result = execute_app("chopper",
                                         "--input",
                                         input_filename.c_str(),
                                         "--number-of-partitions",
                                         "5",
                                         "--output",
                                         layout_filename.c_str());

    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err,
              std::string{"[ERROR] The number of partitions (5) must not exceed the number of user bins (4).\n"});
}

TEST_F(cli_test, chopper_layout_phibf_tmax)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{write_input_file(tmp_dir.path())};

    // The number of technical bins of the top-level IBF of each layout.
    auto top_level_technical_bins = [](std::filesystem::path const & layout_filename)
    {
        std::ifstream layout_stream{layout_filename};
        auto [filenames, config, layouts] = chopper::layout::read_layouts_file(layout_stream);
        std::vector<size_t> result;
        for (auto const & layout : layouts)
        {
            size_t tbs{};
            for (auto const & user_bin : layout.user_bins)
                tbs = std::max(tbs,
                               user_bin.previous_TB_indices.empty()
                                   ? user_bin.storage_TB_id + user_bin.number_of_technical_bins
                                   : user_bin.previous_TB_indices[0] + 1);
            result.push_back(tbs);
        }
        return result;
    };

    for (bool const set_tmax : {false, true})
    {
        std::filesystem::path const layout_filename{tmp_dir.path() / (set_tmax ? "tmax.layout" : "default.layout")};
        std::string const tmax_option{set_tmax ? "--tmax 128" : ""};

        cli_test_result result = execute_app("chopper",
                                             "--input",
                                             input_filename.c_str(),
                                             "--number-of-partitions",
                                             "2",
                                             "--partitioning-approach",
                                             "1",
                                             tmax_option,
                                             "--output",
                                             layout_filename.c_str());
        ASSERT_EQ(result.exit_code, 0) << result.err;

        // Without --tmax, each partition (2 user bins) gets next_multiple_of_64(ceil(sqrt(2))) = 64.
        size_t const expected = set_tmax ? 128u : 64u;
        EXPECT_EQ(top_level_technical_bins(layout_filename), (std::vector<size_t>{expected, expected}));
    }
}
