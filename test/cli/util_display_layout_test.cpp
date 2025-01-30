// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "../api/api_test.hpp"
#include "cli_test.hpp"

std::string get_layout_with_correct_filenames(std::string_view const seq1_filename,
                                              std::string_view const seq2_filename,
                                              std::string_view const seq3_filename,
                                              std::string_view const small_filename,
                                              std::string_view const output_filename)
{
    return {"@CHOPPER_USER_BINS\n"
            "@0 "
            + std::string{seq1_filename}
            + "\n" //
              "@1 "
            + seq2_filename.data() + " " + seq2_filename.data() + // ensure that multi filename works
            +"\n"                                                 //
             "@2 "
            + seq3_filename.data()
            + "\n" //
              "@3 "
            + small_filename.data()
            + "\n" //
              "@CHOPPER_USER_BINS_END\n"
              "@CHOPPER_CONFIG\n"
              "@{\n"
              "@    \"chopper_config\": {\n"
              "@        \"version\": 2,\n"
              "@        \"data_file\": {\n"
              "@            \"value0\": \"foo.txt\"\n"
              "@        },\n"
              "@        \"debug\": false,\n"
              "@        \"sketch_directory\": {\n"
              "@            \"value0\": \"\"\n"
              "@        },\n"
              "@        \"k\": 15,\n"
              "@        \"window_size\": 15,\n"
              "@        \"disable_sketch_output\": true,\n"
              "@        \"precomputed_files\": false,\n"
              "@        \"output_filename\": {\n"
              "@            \"value0\": \""
            + output_filename.data()
            + "\"\n"
              "@        },\n"
              "@        \"determine_best_tmax\": false,\n"
              "@        \"force_all_binnings\": false\n"
              "@    }\n"
              "@}\n"
              "@CHOPPER_CONFIG_END\n"
              "@HIBF_CONFIG\n"
              "@{\n"
              "@    \"hibf_config\": {\n"
              "@        \"version\": 2,\n"
              "@        \"number_of_user_bins\": 4,\n"
              "@        \"number_of_hash_functions\": 2,\n"
              "@        \"maximum_fpr\": 0.05,\n"
              "@        \"relaxed_fpr\": 0.3,\n"
              "@        \"threads\": 2,\n"
              "@        \"sketch_bits\": 12,\n"
              "@        \"tmax\": 4,\n"
              "@        \"empty_bin_fraction\": 0.0,\n"
              "@        \"alpha\": 1.2,\n"
              "@        \"max_rearrangement_ratio\": 0.5,\n"
              "@        \"disable_estimate_union\": false,\n"
              "@        \"disable_rearrangement\": false\n"
              "@    }\n"
              "@}\n"
              "@HIBF_CONFIG_END\n"
              "#TOP_LEVEL_IBF fullest_technical_bin_idx:0\n"
              "#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0\n"
              "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
              "0\t0;0\t2\n"
              "1\t0;2\t2\n"
              "2\t1\t1\n"
              "3\t2\t2\n"};
}

TEST_F(cli_test, display_layout_general)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_filename{tmp_dir.path() / "small.layout"};
    std::filesystem::path const general_filename{tmp_dir.path() / "small.layout.general"};

    {
        std::ofstream fout{layout_filename};
        fout << get_layout_with_correct_filenames(seq1_filename,
                                                  seq2_filename,
                                                  seq3_filename,
                                                  small_filename,
                                                  layout_filename.string());
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    cli_test_result result = execute_app("display_layout",
                                         "general",
                                         "--input",
                                         layout_filename.c_str(),
                                         "--output",
                                         general_filename.c_str());

    ASSERT_EQ(result.exit_code, 0) << "PWD: " << result.pwd << "\nCMD: " << result.command;
    EXPECT_EQ(result.out, std::string{});
    // std err will have a progress bar

    ASSERT_TRUE(std::filesystem::exists(general_filename));

    std::string expected_general_file{
        "# Layout: " + layout_filename.string() + "\n" +
        R"(tb_index	exact_size	estimated_size	fpr_corrected_size	shared_size	ub_count	kind	splits
0	479	483	153	0	2	merged	1
1	466	466	466	0	1	split	1
2	287	289	420	0	1	split	2
3	287	289	420	0	1	split	0
)"};

    std::string const actual_file{string_from_file(general_filename)};
    EXPECT_EQ(expected_general_file, actual_file);
}

TEST_F(cli_test, display_layout_general_with_shared_kmers)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    std::filesystem::path const layout_filename{"small.layout"};
    std::filesystem::path const general_filename{"small.layout.general"};

    {
        std::ofstream fout{layout_filename};
        fout << get_layout_with_correct_filenames(seq1_filename,
                                                  seq2_filename,
                                                  seq3_filename,
                                                  small_filename,
                                                  layout_filename.string());
    }

    cli_test_result result = execute_app("display_layout",
                                         "general",
                                         "--output-shared-kmers",
                                         "--input",
                                         layout_filename.c_str(),
                                         "--output",
                                         general_filename.c_str());

    ASSERT_EQ(result.exit_code, 0) << "PWD: " << result.pwd << "\nCMD: " << result.command;
    EXPECT_EQ(result.out, std::string{});
    // std err will have a progress bar

    ASSERT_TRUE(std::filesystem::exists(general_filename));

    std::string expected_general_file{
        "# Layout: " + layout_filename.string() + "\n" +
        R"(tb_index	exact_size	estimated_size	fpr_corrected_size	shared_size	ub_count	kind	splits
0	479	483	153	371	2	merged	1
1	466	466	466	0	1	split	1
2	287	289	420	0	1	split	2
3	287	289	420	0	1	split	0
)"};

    std::string const actual_file{string_from_file(general_filename)};
    EXPECT_EQ(expected_general_file, actual_file);
}

TEST_F(cli_test, display_layout_sizes)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_filename{tmp_dir.path() / "small.layout"};
    std::filesystem::path const sizes_filename{tmp_dir.path() / "small.layout.sizes"};

    {
        std::ofstream fout{layout_filename};
        fout << get_layout_with_correct_filenames(seq1_filename,
                                                  seq2_filename,
                                                  seq3_filename,
                                                  small_filename,
                                                  layout_filename.string());
    }

    ASSERT_TRUE(std::filesystem::exists(layout_filename));

    cli_test_result result =
        execute_app("display_layout", "sizes", "--input", layout_filename.c_str(), "--output", sizes_filename.c_str());

    ASSERT_EQ(result.exit_code, 0) << "PWD: " << result.pwd << "\nCMD: " << result.command;
    EXPECT_EQ(result.out, std::string{});
    // std err will have a progress bar

    ASSERT_TRUE(std::filesystem::exists(sizes_filename));

    std::string expected_general_file{R"(# Levels: 2
# User bins: 4
LEVEL	BIT_SIZE	IBFS	AVG_LOAD_FACTOR	TBS_TOO_BIG	AVG_TBS_TOO_BIG_ELEMENTS	AVG_MAX_ELEMENTS
0	4832	1	79.44	0	0	479
1	8916	1	80.26	0	0	385
)"};

    std::string const actual_file{string_from_file(sizes_filename)};
    EXPECT_EQ(expected_general_file, actual_file);
}
