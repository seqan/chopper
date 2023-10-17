// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ranges>
#include <string>
#include <vector>

#include <chopper/layout/execute.hpp>

#include "../api_test.hpp"
#include "print_debug_file.hpp"

TEST(execute_test, few_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = (num == 1) ? 880 : 475; // Estimate are 990.71 and 504.88
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = 8;
    config.hibf_config.tmax = 64;
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<std::string> filenames{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"};

    chopper::layout::execute(config, filenames);

    std::string const expected_file{"@CHOPPER_USER_BINS\n"
                                    "@0 seq0\n"
                                    "@1 seq1\n"
                                    "@2 seq2\n"
                                    "@3 seq3\n"
                                    "@4 seq4\n"
                                    "@5 seq5\n"
                                    "@6 seq6\n"
                                    "@7 seq7\n"
                                    "@CHOPPER_USER_BINS_END\n"
                                    "@CHOPPER_CONFIG\n"
                                    "@{\n"
                                    "@    \"chopper_config\": {\n"
                                    "@        \"version\": 2,\n"
                                    "@        \"data_file\": {\n"
                                    "@            \"value0\": \"\"\n"
                                    "@        },\n"
                                    "@        \"debug\": false,\n"
                                    "@        \"sketch_directory\": {\n"
                                    "@            \"value0\": \"\"\n"
                                    "@        },\n"
                                    "@        \"k\": 19,\n"
                                    "@        \"disable_sketch_output\": true,\n"
                                    "@        \"precomputed_files\": false,\n"
                                    "@        \"output_filename\": {\n"
                                    "@            \"value0\": \""
                                    + layout_file.string()
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
                                      "@        \"number_of_user_bins\": 8,\n"
                                      "@        \"number_of_hash_functions\": 2,\n"
                                      "@        \"maximum_false_positive_rate\": 0.05,\n"
                                      "@        \"threads\": 1,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 64,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": true,\n"
                                      "@        \"disable_rearrangement\": true\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:42\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
                                      "7\t0\t6\n"
                                      "6\t6\t6\n"
                                      "5\t12\t6\n"
                                      "4\t18\t6\n"
                                      "3\t24\t6\n"
                                      "2\t30\t6\n"
                                      "0\t36\t6\n"
                                      "1\t42\t22\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file) << actual_file;
}

TEST(execute_test, set_default_tmax)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = (num == 1) ? 1000 : 500;
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{}; // tmax == 0 triggers to set default to the sqrt(#samples)
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = 8;
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<std::string> filenames{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"};

    chopper::layout::execute(config, filenames);

    EXPECT_EQ(config.hibf_config.tmax, 64u);
}

TEST(execute_test, many_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    std::vector<std::string> many_filenames;

    for (size_t i{0}; i < 96u; ++i)
        many_filenames.push_back(seqan3::detail::to_string("seq", i));

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 101 * ((num + 20) / 20);
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.hibf_config.tmax = 64;
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = many_filenames.size();
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    chopper::layout::execute(config, many_filenames);

    std::string const expected_file{"@CHOPPER_USER_BINS\n"
                                    "@0 seq0\n"
                                    "@1 seq1\n"
                                    "@2 seq2\n"
                                    "@3 seq3\n"
                                    "@4 seq4\n"
                                    "@5 seq5\n"
                                    "@6 seq6\n"
                                    "@7 seq7\n"
                                    "@8 seq8\n"
                                    "@9 seq9\n"
                                    "@10 seq10\n"
                                    "@11 seq11\n"
                                    "@12 seq12\n"
                                    "@13 seq13\n"
                                    "@14 seq14\n"
                                    "@15 seq15\n"
                                    "@16 seq16\n"
                                    "@17 seq17\n"
                                    "@18 seq18\n"
                                    "@19 seq19\n"
                                    "@20 seq20\n"
                                    "@21 seq21\n"
                                    "@22 seq22\n"
                                    "@23 seq23\n"
                                    "@24 seq24\n"
                                    "@25 seq25\n"
                                    "@26 seq26\n"
                                    "@27 seq27\n"
                                    "@28 seq28\n"
                                    "@29 seq29\n"
                                    "@30 seq30\n"
                                    "@31 seq31\n"
                                    "@32 seq32\n"
                                    "@33 seq33\n"
                                    "@34 seq34\n"
                                    "@35 seq35\n"
                                    "@36 seq36\n"
                                    "@37 seq37\n"
                                    "@38 seq38\n"
                                    "@39 seq39\n"
                                    "@40 seq40\n"
                                    "@41 seq41\n"
                                    "@42 seq42\n"
                                    "@43 seq43\n"
                                    "@44 seq44\n"
                                    "@45 seq45\n"
                                    "@46 seq46\n"
                                    "@47 seq47\n"
                                    "@48 seq48\n"
                                    "@49 seq49\n"
                                    "@50 seq50\n"
                                    "@51 seq51\n"
                                    "@52 seq52\n"
                                    "@53 seq53\n"
                                    "@54 seq54\n"
                                    "@55 seq55\n"
                                    "@56 seq56\n"
                                    "@57 seq57\n"
                                    "@58 seq58\n"
                                    "@59 seq59\n"
                                    "@60 seq60\n"
                                    "@61 seq61\n"
                                    "@62 seq62\n"
                                    "@63 seq63\n"
                                    "@64 seq64\n"
                                    "@65 seq65\n"
                                    "@66 seq66\n"
                                    "@67 seq67\n"
                                    "@68 seq68\n"
                                    "@69 seq69\n"
                                    "@70 seq70\n"
                                    "@71 seq71\n"
                                    "@72 seq72\n"
                                    "@73 seq73\n"
                                    "@74 seq74\n"
                                    "@75 seq75\n"
                                    "@76 seq76\n"
                                    "@77 seq77\n"
                                    "@78 seq78\n"
                                    "@79 seq79\n"
                                    "@80 seq80\n"
                                    "@81 seq81\n"
                                    "@82 seq82\n"
                                    "@83 seq83\n"
                                    "@84 seq84\n"
                                    "@85 seq85\n"
                                    "@86 seq86\n"
                                    "@87 seq87\n"
                                    "@88 seq88\n"
                                    "@89 seq89\n"
                                    "@90 seq90\n"
                                    "@91 seq91\n"
                                    "@92 seq92\n"
                                    "@93 seq93\n"
                                    "@94 seq94\n"
                                    "@95 seq95\n"
                                    "@CHOPPER_USER_BINS_END\n"
                                    "@CHOPPER_CONFIG\n"
                                    "@{\n"
                                    "@    \"chopper_config\": {\n"
                                    "@        \"version\": 2,\n"
                                    "@        \"data_file\": {\n"
                                    "@            \"value0\": \"\"\n"
                                    "@        },\n"
                                    "@        \"debug\": false,\n"
                                    "@        \"sketch_directory\": {\n"
                                    "@            \"value0\": \"\"\n"
                                    "@        },\n"
                                    "@        \"k\": 19,\n"
                                    "@        \"disable_sketch_output\": true,\n"
                                    "@        \"precomputed_files\": false,\n"
                                    "@        \"output_filename\": {\n"
                                    "@            \"value0\": \""
                                    + layout_file.string()
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
                                      "@        \"number_of_user_bins\": 96,\n"
                                      "@        \"number_of_hash_functions\": 2,\n"
                                      "@        \"maximum_false_positive_rate\": 0.05,\n"
                                      "@        \"threads\": 1,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 64,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": true,\n"
                                      "@        \"disable_rearrangement\": true\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:26\n"
                                      "#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:14\n"
                                      "#LOWER_LEVEL_IBF_1 fullest_technical_bin_idx:14\n"
                                      "#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:14\n"
                                      "#LOWER_LEVEL_IBF_3 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_4 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_5 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_6 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_7 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_8 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_9 fullest_technical_bin_idx:22\n"
                                      "#LOWER_LEVEL_IBF_26 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_27 fullest_technical_bin_idx:0\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
                                      "15\t0;0\t1;14\n"
                                      "16\t0;14\t1;10\n"
                                      "17\t0;24\t1;10\n"
                                      "18\t0;34\t1;10\n"
                                      "19\t0;44\t1;10\n"
                                      "0\t0;54\t1;10\n"
                                      "9\t1;0\t1;14\n"
                                      "10\t1;14\t1;10\n"
                                      "11\t1;24\t1;10\n"
                                      "12\t1;34\t1;10\n"
                                      "13\t1;44\t1;10\n"
                                      "14\t1;54\t1;10\n"
                                      "3\t2;0\t1;14\n"
                                      "4\t2;14\t1;10\n"
                                      "5\t2;24\t1;10\n"
                                      "6\t2;34\t1;10\n"
                                      "7\t2;44\t1;10\n"
                                      "8\t2;54\t1;10\n"
                                      "33\t3;0\t1;26\n"
                                      "32\t3;26\t1;26\n"
                                      "1\t3;52\t1;6\n"
                                      "2\t3;58\t1;6\n"
                                      "30\t4;0\t1;22\n"
                                      "35\t4;22\t1;21\n"
                                      "34\t4;43\t1;21\n"
                                      "38\t5;0\t1;22\n"
                                      "37\t5;22\t1;21\n"
                                      "36\t5;43\t1;21\n"
                                      "20\t6;0\t1;22\n"
                                      "31\t6;22\t1;21\n"
                                      "39\t6;43\t1;21\n"
                                      "23\t7;0\t1;22\n"
                                      "22\t7;22\t1;21\n"
                                      "21\t7;43\t1;21\n"
                                      "27\t8;0\t1;22\n"
                                      "26\t8;22\t1;21\n"
                                      "25\t8;43\t1;21\n"
                                      "29\t9;0\t1;22\n"
                                      "24\t9;22\t1;21\n"
                                      "28\t9;43\t1;21\n"
                                      "40\t10\t1\n"
                                      "41\t11\t1\n"
                                      "42\t12\t1\n"
                                      "43\t13\t1\n"
                                      "44\t14\t1\n"
                                      "45\t15\t1\n"
                                      "46\t16\t1\n"
                                      "47\t17\t1\n"
                                      "48\t18\t1\n"
                                      "59\t19\t1\n"
                                      "58\t20\t1\n"
                                      "57\t21\t1\n"
                                      "56\t22\t1\n"
                                      "55\t23\t1\n"
                                      "54\t24\t1\n"
                                      "53\t25\t1\n"
                                      "51\t26;0\t1;32\n"
                                      "52\t26;32\t1;32\n"
                                      "49\t27;0\t1;32\n"
                                      "50\t27;32\t1;32\n"
                                      "79\t28\t1\n"
                                      "78\t29\t1\n"
                                      "77\t30\t1\n"
                                      "76\t31\t1\n"
                                      "75\t32\t1\n"
                                      "74\t33\t1\n"
                                      "73\t34\t1\n"
                                      "70\t35\t1\n"
                                      "71\t36\t1\n"
                                      "72\t37\t1\n"
                                      "60\t38\t1\n"
                                      "61\t39\t1\n"
                                      "62\t40\t1\n"
                                      "63\t41\t1\n"
                                      "64\t42\t1\n"
                                      "65\t43\t1\n"
                                      "66\t44\t1\n"
                                      "67\t45\t1\n"
                                      "68\t46\t1\n"
                                      "69\t47\t1\n"
                                      "88\t48\t1\n"
                                      "94\t49\t1\n"
                                      "93\t50\t1\n"
                                      "92\t51\t1\n"
                                      "91\t52\t1\n"
                                      "90\t53\t1\n"
                                      "89\t54\t1\n"
                                      "87\t55\t1\n"
                                      "86\t56\t1\n"
                                      "85\t57\t1\n"
                                      "84\t58\t1\n"
                                      "83\t59\t1\n"
                                      "82\t60\t1\n"
                                      "81\t61\t1\n"
                                      "80\t62\t1\n"
                                      "95\t63\t1\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file) << actual_file << std::endl;
}
