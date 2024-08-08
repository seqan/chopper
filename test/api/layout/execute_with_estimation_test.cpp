// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/sketch/read_hll_files_into.hpp>

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include "../api_test.hpp"

TEST(execute_estimation_test, few_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = (num == 1) ? 1000 : 500;
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.hibf_config.tmax = 64;
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = 8;
    config.determine_best_tmax = true;
    config.disable_sketch_output = true;
    config.output_filename = layout_file;
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<std::vector<std::string>>
        filenames{{"seq0"}, {"seq1"}, {"seq2"}, {"seq3"}, {"seq4"}, {"seq5"}, {"seq6"}, {"seq7"}};

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches);

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    chopper::layout::execute(config, filenames, sketches, union_estimation_timer, rearrangement_timer, dp_algorithm_timer);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 8
## number of hash functions = 2
## maximum false positive rate = 0.05
## relaxed false positive rate = 0.3
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.00	1.00	1.00	15.7KiB
# Best t_max (regarding expected query runtime): 64
)expected_cout");
}

TEST(execute_estimation_test, many_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    std::vector<std::vector<std::string>> many_filenames;

    for (size_t i{0}; i < 96u; ++i)
        many_filenames.push_back({seqan3::detail::to_string("seq", i)});

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 100 * ((num + 20) / 20);
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.determine_best_tmax = true;
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.hibf_config.tmax = 1024;
    config.hibf_config.input_fn = simulated_input;
    config.hibf_config.number_of_user_bins = many_filenames.size();
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches);

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    chopper::layout::execute(config, many_filenames, sketches, union_estimation_timer, rearrangement_timer, dp_algorithm_timer);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 96
## number of hash functions = 2
## maximum false positive rate = 0.05
## relaxed false positive rate = 0.3
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.25	1.00	1.25	75.1KiB
128	1.22	1.24	0.68	0.85	51.4KiB
256	1.33	1.33	0.77	1.02	57.5KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

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
                                    "@        \"window_size\": 19,\n"
                                    "@        \"disable_sketch_output\": true,\n"
                                    "@        \"precomputed_files\": false,\n"
                                    "@        \"output_filename\": {\n"
                                    "@            \"value0\": \""
                                    + layout_file.string()
                                    + "\"\n"
                                      "@        },\n"
                                      "@        \"determine_best_tmax\": true,\n"
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
                                      "@        \"maximum_fpr\": 0.05,\n"
                                      "@        \"relaxed_fpr\": 0.3,\n"
                                      "@        \"threads\": 1,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 128,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": true,\n"
                                      "@        \"disable_rearrangement\": true\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:96\n"
                                      "#LOWER_LEVEL_IBF_14 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_15 fullest_technical_bin_idx:0\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
#ifdef _LIBCPP_VERSION // seqan::hibf::sketch::toolbox::sort_by_cardinalities is not stable
                                      "0\t0\t1\n"
                                      "1\t1\t1\n"
                                      "2\t2\t1\n"
                                      "3\t3\t1\n"
                                      "4\t4\t1\n"
                                      "5\t5\t1\n"
                                      "6\t6\t1\n"
                                      "7\t7\t1\n"
                                      "8\t8\t1\n"
                                      "9\t9\t1\n"
                                      "10\t10\t1\n"
                                      "11\t11\t1\n"
                                      "12\t12\t1\n"
                                      "13\t13\t1\n"
                                      "16\t14;0\t1;22\n"
                                      "15\t14;22\t1;21\n"
                                      "14\t14;43\t1;21\n"
                                      "19\t15;0\t1;22\n"
                                      "18\t15;22\t1;21\n"
                                      "17\t15;43\t1;21\n"
                                      "20\t16\t1\n"
                                      "21\t17\t1\n"
                                      "22\t18\t1\n"
                                      "23\t19\t1\n"
                                      "24\t20\t1\n"
                                      "25\t21\t1\n"
                                      "26\t22\t1\n"
                                      "27\t23\t1\n"
                                      "28\t24\t1\n"
                                      "29\t25\t1\n"
                                      "30\t26\t1\n"
                                      "31\t27\t1\n"
                                      "32\t28\t1\n"
                                      "33\t29\t1\n"
                                      "34\t30\t1\n"
                                      "35\t31\t1\n"
                                      "38\t32\t1\n"
                                      "39\t33\t1\n"
                                      "36\t34\t1\n"
                                      "37\t35\t1\n"
                                      "55\t36\t1\n"
                                      "54\t37\t1\n"
                                      "53\t38\t1\n"
                                      "52\t39\t1\n"
                                      "51\t40\t1\n"
                                      "50\t41\t1\n"
                                      "49\t42\t1\n"
                                      "57\t43\t1\n"
                                      "47\t44\t1\n"
                                      "46\t45\t1\n"
                                      "45\t46\t1\n"
                                      "44\t47\t1\n"
                                      "43\t48\t1\n"
                                      "42\t49\t1\n"
                                      "41\t50\t1\n"
                                      "40\t51\t1\n"
                                      "58\t52\t1\n"
                                      "59\t53\t1\n"
                                      "56\t54\t1\n"
                                      "48\t55\t1\n"
                                      "60\t56\t2\n"
                                      "61\t58\t2\n"
                                      "62\t60\t2\n"
                                      "63\t62\t2\n"
                                      "64\t64\t2\n"
                                      "65\t66\t2\n"
                                      "66\t68\t2\n"
                                      "67\t70\t2\n"
                                      "68\t72\t2\n"
                                      "69\t74\t2\n"
                                      "70\t76\t2\n"
                                      "71\t78\t2\n"
                                      "72\t80\t2\n"
                                      "73\t82\t2\n"
                                      "74\t84\t2\n"
                                      "75\t86\t2\n"
                                      "76\t88\t2\n"
                                      "79\t90\t2\n"
                                      "78\t92\t2\n"
                                      "77\t94\t2\n"
                                      "80\t96\t2\n"
                                      "81\t98\t2\n"
                                      "82\t100\t2\n"
                                      "83\t102\t2\n"
                                      "84\t104\t2\n"
                                      "85\t106\t2\n"
                                      "86\t108\t2\n"
                                      "87\t110\t2\n"
                                      "88\t112\t2\n"
                                      "89\t114\t2\n"
                                      "90\t116\t2\n"
                                      "91\t118\t2\n"
                                      "92\t120\t2\n"
                                      "93\t122\t2\n"
                                      "94\t124\t2\n"
                                      "95\t126\t2\n"};
#else
                                      "0\t0\t1\n"
                                      "19\t1\t1\n"
                                      "18\t2\t1\n"
                                      "17\t3\t1\n"
                                      "16\t4\t1\n"
                                      "15\t5\t1\n"
                                      "14\t6\t1\n"
                                      "13\t7\t1\n"
                                      "12\t8\t1\n"
                                      "11\t9\t1\n"
                                      "10\t10\t1\n"
                                      "9\t11\t1\n"
                                      "8\t12\t1\n"
                                      "7\t13\t1\n"
                                      "4\t14;0\t1;22\n"
                                      "5\t14;22\t1;21\n"
                                      "6\t14;43\t1;21\n"
                                      "1\t15;0\t1;22\n"
                                      "2\t15;22\t1;21\n"
                                      "3\t15;43\t1;21\n"
                                      "32\t16\t1\n"
                                      "33\t17\t1\n"
                                      "34\t18\t1\n"
                                      "35\t19\t1\n"
                                      "30\t20\t1\n"
                                      "36\t21\t1\n"
                                      "37\t22\t1\n"
                                      "38\t23\t1\n"
                                      "39\t24\t1\n"
                                      "31\t25\t1\n"
                                      "20\t26\t1\n"
                                      "21\t27\t1\n"
                                      "22\t28\t1\n"
                                      "23\t29\t1\n"
                                      "25\t30\t1\n"
                                      "26\t31\t1\n"
                                      "27\t32\t1\n"
                                      "28\t33\t1\n"
                                      "24\t34\t1\n"
                                      "29\t35\t1\n"
                                      "40\t36\t1\n"
                                      "41\t37\t1\n"
                                      "42\t38\t1\n"
                                      "43\t39\t1\n"
                                      "44\t40\t1\n"
                                      "45\t41\t1\n"
                                      "46\t42\t1\n"
                                      "47\t43\t1\n"
                                      "48\t44\t1\n"
                                      "59\t45\t1\n"
                                      "58\t46\t1\n"
                                      "57\t47\t1\n"
                                      "56\t48\t1\n"
                                      "55\t49\t1\n"
                                      "54\t50\t1\n"
                                      "53\t51\t1\n"
                                      "52\t52\t1\n"
                                      "51\t53\t1\n"
                                      "50\t54\t1\n"
                                      "49\t55\t1\n"
                                      "79\t56\t2\n"
                                      "78\t58\t2\n"
                                      "77\t60\t2\n"
                                      "76\t62\t2\n"
                                      "75\t64\t2\n"
                                      "74\t66\t2\n"
                                      "73\t68\t2\n"
                                      "70\t70\t2\n"
                                      "71\t72\t2\n"
                                      "72\t74\t2\n"
                                      "60\t76\t2\n"
                                      "61\t78\t2\n"
                                      "62\t80\t2\n"
                                      "63\t82\t2\n"
                                      "64\t84\t2\n"
                                      "65\t86\t2\n"
                                      "66\t88\t2\n"
                                      "67\t90\t2\n"
                                      "68\t92\t2\n"
                                      "69\t94\t2\n"
                                      "88\t96\t2\n"
                                      "94\t98\t2\n"
                                      "93\t100\t2\n"
                                      "92\t102\t2\n"
                                      "91\t104\t2\n"
                                      "90\t106\t2\n"
                                      "89\t108\t2\n"
                                      "87\t110\t2\n"
                                      "86\t112\t2\n"
                                      "85\t114\t2\n"
                                      "84\t116\t2\n"
                                      "83\t118\t2\n"
                                      "82\t120\t2\n"
                                      "81\t122\t2\n"
                                      "80\t124\t2\n"
                                      "95\t126\t2\n"};
#endif
    std::string const actual_file{string_from_file(layout_file)};
    EXPECT_EQ(actual_file, expected_file);
}

TEST(execute_estimation_test, many_ubs_force_all)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    std::vector<std::vector<std::string>> many_filenames;
    std::vector<size_t> many_kmer_counts;

    for (size_t i{0}; i < 96u; ++i)
        many_filenames.push_back({seqan3::detail::to_string("seq", i)});

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 100 * ((num + 20) / 20);
        for (auto hash : std::views::iota(0u, desired_kmer_count))
            it = hash;
    };

    chopper::configuration config{};
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.disable_sketch_output = true;
    config.output_filename = layout_file;
    config.hibf_config.input_fn = simulated_input, config.hibf_config.number_of_user_bins = many_filenames.size(),
    config.hibf_config.tmax = 256;
    config.hibf_config.disable_estimate_union = true; // also disables rearrangement

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches);

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    chopper::layout::execute(config, many_filenames, sketches, union_estimation_timer, rearrangement_timer, dp_algorithm_timer);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 96
## number of hash functions = 2
## maximum false positive rate = 0.05
## relaxed false positive rate = 0.3
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.25	1.00	1.25	75.1KiB
128	1.22	1.24	0.68	0.85	51.4KiB
256	1.33	1.33	0.77	1.02	57.5KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 128,"), std::string::npos);
}

struct dna4_traits3 : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type3 = seqan3::sequence_file_input<dna4_traits3,
                                                        seqan3::fields<seqan3::field::seq>,
                                                        seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

TEST(execute_estimation_test, with_rearrangement)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const sketches_dir{tmp_dir.path() / "test"};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};
    size_t const kmer_size{15};

    std::vector<std::vector<std::string>> filenames{};
    std::vector<std::string> hll_filenames;
    std::vector<size_t> expected_kmer_counts;

    for (size_t i{0}; i < 49u; ++i)
    {
        filenames.push_back({data("seq1.fa").string()});
        filenames.push_back({data("seq2.fa").string()});
        filenames.push_back({data("seq3.fa").string()});
        filenames.push_back({data("small.fa").string()});

        hll_filenames.push_back("seq1.hll");
        hll_filenames.push_back("seq2.hll");
        hll_filenames.push_back("seq3.hll");
        hll_filenames.push_back("small.hll");

        expected_kmer_counts.push_back(387);
        expected_kmer_counts.push_back(470);
        expected_kmer_counts.push_back(465);
        expected_kmer_counts.push_back(578);
    }

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    auto data_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        for (std::string const & filename : filenames[num])
        {
            sequence_file_type3 fin{filename};

            for (auto && [seq] : fin)
            {
                for (auto hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{kmer_size}))
                    it = hash_value;
            }
        }
    };

    chopper::configuration config{};
    config.k = kmer_size;
    config.sketch_directory = sketches_dir;
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.output_filename = layout_file;
    // config.output_verbose_statistics = true;
    config.hibf_config.threads = 1;
    config.hibf_config.tmax = 256;
    config.hibf_config.input_fn = data_input;
    config.hibf_config.number_of_user_bins = filenames.size();

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches);

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    chopper::layout::execute(config, filenames, sketches, union_estimation_timer, rearrangement_timer, dp_algorithm_timer);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    std::string const expected = []()
    {
        std::string result =
            R"expected_cout(## ### Parameters ###
## number of user bins = 196
## number of hash functions = 2
## maximum false positive rate = 0.05
## relaxed false positive rate = 0.3
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
)expected_cout";

#ifdef _LIBCPP_VERSION // seqan::hibf::sketch::toolbox::sort_by_cardinalities is not stable
        result +=
            R"expected_cout(64	1.00	2.08	1.00	2.08	117.1KiB
128	1.22	1.75	1.10	1.93	128.7KiB
256	1.33	1.52	1.18	1.81	138.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout";
#else
        result +=
            R"expected_cout(64	1.00	2.22	1.00	2.22	117.1KiB
128	1.22	1.95	1.15	2.23	134.3KiB
256	1.33	1.52	1.18	1.81	138.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout";
#endif

        return result;
    }();

    EXPECT_EQ(written_file, expected);

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 256,"), std::string::npos);
}
