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

#include <chopper/input_functor.hpp>
#include <chopper/sketch/sketch_file.hpp>

#include <hibf/sketch/compute_sketches.hpp>

#include "../api/api_test.hpp" // for string_from_file
#include "cli_test.hpp"

TEST_F(cli_test, chopper_layout_from_sketch_file)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const input_filename{tmp_dir.path() / "data.sketches"};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};

    // create sketch_file
    {
        chopper::sketch::sketch_file sout{};

        sout.filenames = {{seq1_filename}, {seq2_filename}, {seq3_filename}, {seq4_filename}};

        sout.chopper_config.k = 19;
        sout.chopper_config.hibf_config.sketch_bits = 12;
        sout.chopper_config.hibf_config.input_fn =
            chopper::input_functor{sout.filenames, false, sout.chopper_config.k, sout.chopper_config.window_size};
        sout.chopper_config.hibf_config.number_of_user_bins = sout.filenames.size();

        seqan::hibf::sketch::compute_sketches(sout.chopper_config.hibf_config, sout.hll_sketches);

        std::ofstream os{input_filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(sout);
    }

    cli_test_result result = execute_app("chopper",
                                         "--threads",
                                         "2",
                                         "--input",
                                         input_filename.c_str(),
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
                                    + input_filename.string()
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

    // ------------------------------------------------------------------------------------------------
    // Tests for options. Needs valid sketch file as input. Hence, we reuse the sketch file from above.
    // ------------------------------------------------------------------------------------------------

    cli_test_result result2 = execute_app("chopper",
                                          "--threads 2",
                                          "--kmer 20",
                                          "--input",
                                          input_filename.c_str(),
                                          "--tmax 64",
                                          "--output",
                                          binning_filename.c_str());

    EXPECT_EQ(result2.exit_code, 0);
    EXPECT_EQ(result2.out, std::string{});
    EXPECT_EQ(
        result2.err,
        std::string{"[WARNING] Given k-mer size (20) differs from k-mer size in the sketch file (20). The results may "
                    "be suboptimal. If this was a conscious decision, you can ignore this warning.\n"});

    cli_test_result result3 = execute_app("chopper",
                                          "--threads 2",
                                          "--input",
                                          input_filename.c_str(),
                                          "--tmax 64",
                                          "--output",
                                          binning_filename.c_str(),
                                          "--sketch-bits 12");
    EXPECT_NE(result3.exit_code, 0);
    EXPECT_EQ(result3.out, std::string{});
    EXPECT_EQ(result3.err, std::string{"[ERROR] You cannot set --sketch-bits when using a sketch file as input.\n"});
}
