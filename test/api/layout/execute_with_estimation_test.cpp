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
        size_t const desired_kmer_count = (num == 1) ? 1700 : 1200;
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

    chopper::layout::execute(config, filenames, sketches);

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
64	1.00	1.00	1.00	1.00	38.2KiB
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

    // Creates sizes of the following series
    // [801,802,...,820,922,923,...,941,1043,1044,...,1062,1164,1165,...,1183,1285,1286,...,1300]
    // See also https://godbolt.org/z/9517eaaaG
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 101 * ((num + 20) / 20) + num + 700;
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

    chopper::layout::execute(config, many_filenames, sketches);

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
64	1.00	1.36	1.00	1.36	269.0KiB
128	1.22	1.42	1.00	1.42	269.8KiB
# Best t_max (regarding expected query runtime): 64
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
                                      "@        \"version\": 3,\n"
                                      "@        \"number_of_user_bins\": 96,\n"
                                      "@        \"number_of_hash_functions\": 2,\n"
                                      "@        \"maximum_fpr\": 0.05,\n"
                                      "@        \"relaxed_fpr\": 0.3,\n"
                                      "@        \"threads\": 1,\n"
                                      "@        \"sketch_bits\": 12,\n"
                                      "@        \"tmax\": 64,\n"
                                      "@        \"empty_bin_fraction\": 0.0,\n"
                                      "@        \"track_occupancy\": false,\n"
                                      "@        \"alpha\": 1.2,\n"
                                      "@        \"max_rearrangement_ratio\": 0.5,\n"
                                      "@        \"disable_estimate_union\": true,\n"
                                      "@        \"disable_rearrangement\": true\n"
                                      "@    }\n"
                                      "@}\n"
                                      "@HIBF_CONFIG_END\n"
                                      "#TOP_LEVEL_IBF fullest_technical_bin_idx:63\n"
                                      "#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_1 fullest_technical_bin_idx:52\n"
                                      "#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:52\n"
                                      "#LOWER_LEVEL_IBF_3 fullest_technical_bin_idx:52\n"
                                      "#LOWER_LEVEL_IBF_4 fullest_technical_bin_idx:31\n"
                                      "#LOWER_LEVEL_IBF_5 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_6 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_7 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_8 fullest_technical_bin_idx:0\n"
                                      "#LOWER_LEVEL_IBF_9 fullest_technical_bin_idx:0\n"
                                      "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
                                      "1\t0;0\t1;32\n"
                                      "0\t0;32\t1;32\n"
                                      "6\t1;0\t1;13\n"
                                      "5\t1;13\t1;13\n"
                                      "4\t1;26\t1;13\n"
                                      "3\t1;39\t1;13\n"
                                      "2\t1;52\t1;12\n"
                                      "11\t2;0\t1;13\n"
                                      "10\t2;13\t1;13\n"
                                      "9\t2;26\t1;13\n"
                                      "8\t2;39\t1;13\n"
                                      "7\t2;52\t1;12\n"
                                      "16\t3;0\t1;13\n"
                                      "15\t3;13\t1;13\n"
                                      "14\t3;26\t1;13\n"
                                      "13\t3;39\t1;13\n"
                                      "12\t3;52\t1;12\n"
                                      "21\t4;0\t1;16\n"
                                      "20\t4;16\t1;15\n"
                                      "19\t4;31\t1;11\n"
                                      "18\t4;42\t1;11\n"
                                      "17\t4;53\t1;11\n"
                                      "25\t5;0\t1;16\n"
                                      "24\t5;16\t1;16\n"
                                      "23\t5;32\t1;16\n"
                                      "22\t5;48\t1;16\n"
                                      "29\t6;0\t1;16\n"
                                      "28\t6;16\t1;16\n"
                                      "27\t6;32\t1;16\n"
                                      "26\t6;48\t1;16\n"
                                      "33\t7;0\t1;16\n"
                                      "32\t7;16\t1;16\n"
                                      "31\t7;32\t1;16\n"
                                      "30\t7;48\t1;16\n"
                                      "37\t8;0\t1;16\n"
                                      "36\t8;16\t1;16\n"
                                      "35\t8;32\t1;16\n"
                                      "34\t8;48\t1;16\n"
                                      "41\t9;0\t1;18\n"
                                      "40\t9;18\t1;18\n"
                                      "39\t9;36\t1;14\n"
                                      "38\t9;50\t1;14\n"
                                      "42\t10\t1\n"
                                      "43\t11\t1\n"
                                      "44\t12\t1\n"
                                      "45\t13\t1\n"
                                      "46\t14\t1\n"
                                      "47\t15\t1\n"
                                      "48\t16\t1\n"
                                      "49\t17\t1\n"
                                      "50\t18\t1\n"
                                      "51\t19\t1\n"
                                      "52\t20\t1\n"
                                      "53\t21\t1\n"
                                      "54\t22\t1\n"
                                      "55\t23\t1\n"
                                      "56\t24\t1\n"
                                      "57\t25\t1\n"
                                      "58\t26\t1\n"
                                      "59\t27\t1\n"
                                      "60\t28\t1\n"
                                      "61\t29\t1\n"
                                      "62\t30\t1\n"
                                      "63\t31\t1\n"
                                      "64\t32\t1\n"
                                      "65\t33\t1\n"
                                      "66\t34\t1\n"
                                      "67\t35\t1\n"
                                      "68\t36\t1\n"
                                      "69\t37\t1\n"
                                      "70\t38\t1\n"
                                      "71\t39\t1\n"
                                      "72\t40\t1\n"
                                      "73\t41\t1\n"
                                      "74\t42\t1\n"
                                      "75\t43\t1\n"
                                      "76\t44\t1\n"
                                      "77\t45\t1\n"
                                      "78\t46\t1\n"
                                      "79\t47\t1\n"
                                      "80\t48\t1\n"
                                      "81\t49\t1\n"
                                      "82\t50\t1\n"
                                      "83\t51\t1\n"
                                      "84\t52\t1\n"
                                      "85\t53\t1\n"
                                      "86\t54\t1\n"
                                      "87\t55\t1\n"
                                      "88\t56\t1\n"
                                      "89\t57\t1\n"
                                      "90\t58\t1\n"
                                      "91\t59\t1\n"
                                      "92\t60\t1\n"
                                      "93\t61\t1\n"
                                      "94\t62\t1\n"
                                      "95\t63\t1\n"};
    std::string const actual_file{string_from_file(layout_file)};
    EXPECT_EQ(actual_file, expected_file) << actual_file;
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

    // There are 20 files with counts of {800,900,1000,1100} each. There are 16 files with count 1200.
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        size_t const desired_kmer_count = 100 * ((num + 20) / 20) + 700;
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

    chopper::layout::execute(config, many_filenames, sketches);

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
64	1.00	1.37	1.00	1.37	256.7KiB
128	1.22	1.43	0.99	1.42	255.3KiB
256	1.33	1.33	0.81	1.07	206.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout");

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 256,"), std::string::npos);
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

    for (size_t i{0}; i < 49u; ++i)
    {
        filenames.push_back({data("seq1.fa").string()});
        filenames.push_back({data("seq2.fa").string()});
        filenames.push_back({data("seq3.fa").string()});
        filenames.push_back({data("small.fa").string()});
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

    chopper::layout::execute(config, filenames, sketches);

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
            R"expected_cout(64	1.00	2.26	1.00	2.26	110.5KiB
128	1.22	1.96	1.08	2.11	118.9KiB
256	1.33	1.66	1.18	1.97	130.9KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout";
#else
        result +=
            R"expected_cout(64	1.00	2.34	1.00	2.34	110.1KiB
128	1.22	2.01	1.08	2.17	118.9KiB
256	1.33	1.76	1.39	2.44	152.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout";
#endif
        return result;
    }();

    EXPECT_EQ(written_file, expected);

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 256,"), std::string::npos);
}
