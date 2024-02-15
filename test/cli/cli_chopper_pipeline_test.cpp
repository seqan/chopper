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

#include "../api/api_test.hpp"
#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_layout)
{
    std::string const seq_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const taxa_filename{tmp_dir.path() / "data.tsv"};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename};
        fout << seq_filename << '\t' << "TAX1\n"
             << seq_filename << '\t'
             << "TAX2\n"
             /* << seq_filename << '\t' << "TAX2\n" */
             << seq_filename << '\t' << "TAX3\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--kmer",
                                         "15",
                                         "--threads",
                                         "2",
                                         "--input",
                                         taxa_filename.c_str(),
                                         "--tmax",
                                         "64",
                                         "--output",
                                         binning_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected_file{"@CHOPPER_USER_BINS\n"
                                    "@0 "
                                    + seq_filename
                                    + "\n"
                                      "@1 "
                                    + seq_filename
                                    + "\n"
                                      "@2 "
                                    + seq_filename
                                    + "\n"
                                      "@CHOPPER_USER_BINS_END\n"
                                      "@CHOPPER_CONFIG\n"
                                      "@{\n"
                                      "@    \"chopper_config\": {\n"
                                      "@        \"version\": 2,\n"
                                      "@        \"data_file\": {\n"
                                      "@            \"value0\": \""
                                    + taxa_filename.string()
                                    + "\"\n"
                                      "@        },\n"
                                      "@        \"debug\": false,\n"
                                      "@        \"sketch_directory\": {\n"
                                      "@            \"value0\": \"\"\n"
                                      "@        },\n"
                                      "@        \"k\": 15,\n"
                                      "@        \"window_size\": 15,\n"
                                      "@        \"disable_sketch_output\": true,\n"
                                      "@        \"precomputed_files\": false,\n"
                                      "@        \"maximum_index_size\": 0,\n"
                                      "@        \"number_of_partitions\": 0,\n"
                                      "@        \"output_filename\": {\n"
                                      "@            \"value0\": \""
                                    + binning_filename.string()
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
                                      "@        \"version\": 1,\n"
                                      "@        \"number_of_user_bins\": 3,\n"
                                      "@        \"number_of_hash_functions\": 2,\n"
                                      "@        \"maximum_fpr\": 0.05,\n"
                                      "@        \"relaxed_fpr\": 0.3,\n"
                                      "@        \"threads\": 2,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 64,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": false,\n"
                                      "@        \"disable_rearrangement\": false\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:22\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"};
    // Order of IDs is not deterministic
    // "2\t0\t22\n"
    // "1\t22\t21\n"
    // "0\t43\t21\n"};

    std::string const actual_file{string_from_file(binning_filename)};
    EXPECT_TRUE(actual_file.starts_with(expected_file));
}

TEST_F(cli_test, chopper_layout2)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const taxa_filename{tmp_dir.path() / "data.tsv"};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};

    // we need to have filenames from the user
    {
        std::ofstream fout{taxa_filename};
        fout << seq1_filename << '\n' << seq2_filename << '\n' << seq3_filename << '\n' << seq4_filename << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--threads",
                                         "2",
                                         "--sketch-bits",
                                         "12",
                                         "--input",
                                         taxa_filename.c_str(),
                                         "--tmax",
                                         "64",
                                         "--output",
                                         binning_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string const expected_file{"@CHOPPER_USER_BINS\n"
                                    "@0 "
                                    + seq1_filename
                                    + "\n"
                                      "@1 "
                                    + seq2_filename
                                    + "\n"
                                      "@2 "
                                    + seq3_filename
                                    + "\n"
                                      "@3 "
                                    + seq4_filename
                                    + "\n"
                                      "@CHOPPER_USER_BINS_END\n"
                                      "@CHOPPER_CONFIG\n"
                                      "@{\n"
                                      "@    \"chopper_config\": {\n"
                                      "@        \"version\": 2,\n"
                                      "@        \"data_file\": {\n"
                                      "@            \"value0\": \""
                                    + taxa_filename.string()
                                    + "\"\n"
                                      "@        },\n"
                                      "@        \"debug\": false,\n"
                                      "@        \"sketch_directory\": {\n"
                                      "@            \"value0\": \"\"\n"
                                      "@        },\n"
                                      "@        \"k\": 19,\n"
                                      "@        \"window_size\": 19,\n"
                                      "@        \"disable_sketch_output\": true,\n"
                                      "@        \"precomputed_files\": false,\n"
                                      "@        \"maximum_index_size\": 0,\n"
                                      "@        \"number_of_partitions\": 0,\n"
                                      "@        \"output_filename\": {\n"
                                      "@            \"value0\": \""
                                    + binning_filename.string()
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
                                      "@        \"version\": 1,\n"
                                      "@        \"number_of_user_bins\": 4,\n"
                                      "@        \"number_of_hash_functions\": 2,\n"
                                      "@        \"maximum_fpr\": 0.05,\n"
                                      "@        \"relaxed_fpr\": 0.3,\n"
                                      "@        \"threads\": 2,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 64,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": false,\n"
                                      "@        \"disable_rearrangement\": false\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:16\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
                                      "1\t0\t16\n"
                                      "3\t16\t23\n"
                                      "2\t39\t15\n"
                                      "0\t54\t10\n"};

    std::string const actual_file{string_from_file(binning_filename)};
    EXPECT_EQ(actual_file, expected_file);
}
