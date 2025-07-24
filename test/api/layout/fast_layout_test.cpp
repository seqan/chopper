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

#include <hibf/sketch/compute_sketches.hpp>

#include "../api_test.hpp"
#include "print_debug_file.hpp"

TEST(execute_test, many_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    std::vector<std::vector<std::string>> many_filenames;

    for (size_t i{0}; i < 96u; ++i)
        many_filenames.push_back({seqan3::detail::to_string("seq", i)});

    // Creates sizes of the following series
    // [100,101,...,120,222,223,...,241,343,344,...362,464,465,...,483,585,586,...,600]
    // See also https://godbolt.org/z/9517eaaaG
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 1001 * ((num + 20) / 20) + num;
        for (auto hash : std::views::iota(num * 600, num * 600 + desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.hibf_config.tmax = 64;
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = many_filenames.size();
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

    chopper::layout::execute(config, many_filenames, sketches, minHash_sketches);

    std::string const expected_file{
    "@CHOPPER_USER_BINS\n"
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
    "@        \"window_size\": 19,\n"
    "@        \"disable_sketch_output\": true,\n"
    "@        \"precomputed_files\": false,\n"
    "@        \"output_filename\": {\n"
    "@            \"value0\": \"" +  layout_file.string() + "\"\n"
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
    "@        \"number_of_user_bins\": 96,\n"
    "@        \"number_of_hash_functions\": 2,\n"
    "@        \"maximum_fpr\": 0.05,\n"
    "@        \"relaxed_fpr\": 0.3,\n"
    "@        \"threads\": 1,\n"
    "@        \"sketch_bits\": 12,\n"
    "@        \"tmax\": 64,\n"
    "@        \"empty_bin_fraction\": 0.0,\n"
    "@        \"alpha\": 1.2,\n"
    "@        \"max_rearrangement_ratio\": 0.5,\n"
    "@        \"disable_estimate_union\": true,\n"
    "@        \"disable_rearrangement\": true\n"
    "@    }\n"
    "@}\n"
    "@HIBF_CONFIG_END\n"
    "#TOP_LEVEL_IBF fullest_technical_bin_idx:23\n"
    "#LOWER_LEVEL_IBF_1 fullest_technical_bin_idx:24\n"
    "#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:28\n"
    "#LOWER_LEVEL_IBF_3 fullest_technical_bin_idx:6\n"
    "#LOWER_LEVEL_IBF_4 fullest_technical_bin_idx:15\n"
    "#LOWER_LEVEL_IBF_5 fullest_technical_bin_idx:0\n"
    "#LOWER_LEVEL_IBF_6 fullest_technical_bin_idx:24\n"
    "#LOWER_LEVEL_IBF_8 fullest_technical_bin_idx:24\n"
    "#LOWER_LEVEL_IBF_9 fullest_technical_bin_idx:35\n"
    "#LOWER_LEVEL_IBF_12 fullest_technical_bin_idx:7\n"
    "#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS\n"
    "0	5;0	1;2\n"
    "1	9;0	1;2\n"
    "2	3;0	1;2\n"
    "3	12;0	1;2\n"
    "4	12;2	1;2\n"
    "5	6;0	1;2\n"
    "6	9;2	1;2\n"
    "7	9;4	1;2\n"
    "8	3;2	1;2\n"
    "9	3;4	1;2\n"
    "10	4;0	1;2\n"
    "11	6;2	1;2\n"
    "12	2;0	1;2\n"
    "13	2;2	1;2\n"
    "14	2;4	1;2\n"
    "15	4;2	1;2\n"
    "16	12;4	1;3\n"
    "17	6;4	1;2\n"
    "18	5;2	1;3\n"
    "19	9;6	1;2\n"
    "20	8;0	1;12\n"
    "21	8;12	1;12\n"
    "22	6;6	1;9\n"
    "23	6;15	1;9\n"
    "24	9;8	1;9\n"
    "25	9;17	1;9\n"
    "26	4;4	1;11\n"
    "27	4;15	1;11\n"
    "28	1;0	1;12\n"
    "29	1;12	1;12\n"
    "30	2;6	1;11\n"
    "31	2;17	1;11\n"
    "32	12;7	1;13\n"
    "33	6;24	1;9\n"
    "34	5;5	1;14\n"
    "35	9;26	1;9\n"
    "36	3;6	1;13\n"
    "37	22	1\n"
    "38	21	1\n"
    "39	20	1\n"
    "40	15	1\n"
    "41	18	1\n"
    "42	19	1\n"
    "43	13	1\n"
    "44	17	1\n"
    "45	14	1\n"
    "46	10	1\n"
    "47	16	1\n"
    "48	11	1\n"
    "49	8;24	1;40\n"
    "50	12;20	1;44\n"
    "51	6;33	1;31\n"
    "52	5;19	1;45\n"
    "53	9;35	1;29\n"
    "54	3;19	1;45\n"
    "55	4;26	1;38\n"
    "56	7	1\n"
    "57	1;24	1;40\n"
    "58	0	1\n"
    "59	2;28	1;36\n"
    "60	63	1\n"
    "61	62	1\n"
    "62	60	1\n"
    "63	61	1\n"
    "64	59	1\n"
    "65	58	1\n"
    "66	57	1\n"
    "67	55	1\n"
    "68	56	1\n"
    "69	54	1\n"
    "70	52	1\n"
    "71	53	1\n"
    "72	51	1\n"
    "73	49	1\n"
    "74	48	1\n"
    "75	50	1\n"
    "76	45	1\n"
    "77	47	1\n"
    "78	46	1\n"
    "79	44	1\n"
    "80	39	1\n"
    "81	40	1\n"
    "82	41	1\n"
    "83	42	1\n"
    "84	38	1\n"
    "85	37	1\n"
    "86	43	1\n"
    "87	36	1\n"
    "88	35	1\n"
    "89	34	1\n"
    "90	29	2\n"
    "91	31	2\n"
    "92	27	2\n"
    "93	23	2\n"
    "94	33	1\n"
    "95	25	2\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file) << actual_file << std::endl;
}
