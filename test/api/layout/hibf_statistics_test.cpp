// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <ranges>
#include <string>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/hibf_statistics.hpp>

#include "../api_test.hpp"

TEST(byte_size_to_formatted_str, storage_unit)
{
    // `ULL` means unsigned long long; ensures that shifting works since `55` will usually be a 32bit-int
    // `<< 10` is the same as `* 1024`, `<< 20` the same as `* 1024 * 1024`, and so on
    EXPECT_EQ("8Bytes", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL));
    EXPECT_EQ("8.0KiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 10));
    EXPECT_EQ("8.0MiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 20));
    EXPECT_EQ("8.0GiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 30));
    EXPECT_EQ("8.0TiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 40));
    EXPECT_EQ("8.0PiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 50));
    EXPECT_EQ("8.0EiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(8ULL << 60));
}

TEST(byte_size_to_formatted_str, rounding)
{
    EXPECT_EQ("5.8GiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(6'174'015'488ULL));
    EXPECT_EQ("5.7GiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(6'174'015'487ULL));
    // This is 1 bytes short of 1MiB. Rounding should change the unit from KiB to MiB.
    EXPECT_EQ("1.0MiB", chopper::layout::hibf_statistics::byte_size_to_formatted_str(1'048'575ULL));
}

TEST(byte_size_to_formatted_str, edge_cases)
{
    EXPECT_EQ("0Bytes", chopper::layout::hibf_statistics::byte_size_to_formatted_str(0ULL));
    EXPECT_EQ("16.0EiB",
              chopper::layout::hibf_statistics::byte_size_to_formatted_str(std::numeric_limits<size_t>::max()));
}

TEST(hibf_statistics, only_merged_on_top_level)
{
    // parameters for this test HIBF
    size_t const num_top_level_bins = 4u;

    chopper::configuration config{}; // default config
    config.hibf_config.tmax = 64u;
    config.hibf_config.disable_estimate_union = true; /* also disable rearrangement */

    std::vector<std::string> filenames{"s1", "s2"};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{{}, {}};
    std::vector<size_t> kmer_counts{50, 50};

    chopper::layout::hibf_statistics stats(config, sketches, kmer_counts);

    for (size_t i = 0; i < num_top_level_bins; ++i)
    {
        stats.hibf_layout.max_bins.emplace_back(std::vector<size_t>{0u}, 0u);
        stats.hibf_layout.max_bins.emplace_back(std::vector<size_t>{1u}, 0u);

        stats.hibf_layout.user_bins.emplace_back(std::vector<size_t>{0u}, 0u, 2u, 0u);
        stats.hibf_layout.user_bins.emplace_back(std::vector<size_t>{1u}, 0u, 2u, 1u);
    }

    testing::internal::CaptureStdout();

    stats.print_header_to(std::cout);
    size_t max_64{};
    stats.print_summary_to(max_64, std::cout);
    std::cout.flush();

    std::string summary = testing::internal::GetCapturedStdout();
    std::string expected_cout =
        R"expected_cout(## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
## uncorr_size : The expected size of an tmax-HIBF without FPR correction
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size	uncorr_size	level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
64	1.00	8.00	1.00	8.00	979Bytes	790Bytes	:0:1	:1:2	:395Bytes:584Bytes	:395Bytes:395Bytes	:2:16	:2:8	:0.00:100.00	:-:2	:-:2.00	:-:1.46	:-:1.48
)expected_cout";

    EXPECT_EQ(summary, expected_cout);
}

TEST(execute_test, chopper_layout_statistics)
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

    chopper::configuration config{.data_file = "not needed",
                                  .disable_sketch_output = true,
                                  .output_filename = layout_file.c_str(),
                                  .output_verbose_statistics = true,
                                  .hibf_config = {.input_fn = simulated_input,
                                                  .number_of_user_bins = many_filenames.size(),
                                                  .tmax = 64,
                                                  .disable_estimate_union = true /* also disable rearrangement */}};

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    chopper::layout::execute(config, many_filenames);
    std::string layout_result_stdout = testing::internal::GetCapturedStdout();
    std::string layout_result_stderr = testing::internal::GetCapturedStderr();

    std::string expected_cout =
        R"expected_cout(## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
## uncorr_size : The expected size of an tmax-HIBF without FPR correction
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size	uncorr_size	level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
64	1.00	1.25	1.00	1.25	75.5KiB	46.8KiB	:0:1	:1:12	:38.8KiB:36.7KiB	:38.8KiB:8.0KiB	:64:768	:64:64	:81.25:100.00	:1:32	:1.00:17.45	:1.00:6.20	:1.00:4.84
)expected_cout";

    EXPECT_EQ(layout_result_stdout, expected_cout) << layout_result_stdout;
    EXPECT_EQ(layout_result_stderr, std::string{});
}

TEST(execute_test, chopper_layout_statistics_determine_best_bins)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};
    std::filesystem::path const stats_file{binning_filename.string() + ".stats"};

    std::vector<std::string> filenames{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7", "seq8", "seq9"};

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    auto simulated_input = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        std::vector<size_t> kmer_counts{10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};
        for (auto hash : std::views::iota(0u, kmer_counts[num]))
            it = hash;
    };

    chopper::configuration config{.data_file = "not needed",
                                  .disable_sketch_output = true,
                                  .output_filename = binning_filename.c_str(),
                                  .determine_best_tmax = true,
                                  .force_all_binnings = true,
                                  .output_verbose_statistics = true,
                                  .hibf_config = {.input_fn = simulated_input,
                                                  .number_of_user_bins = filenames.size(),
                                                  .tmax = 128,
                                                  .disable_estimate_union = true /* also disable rearrangement */}};

    chopper::layout::execute(config, filenames);

    std::string expected_cout =
        R"expected_cout(## ### Parameters ###
## number of user bins = 10
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
## uncorr_size : The expected size of an tmax-HIBF without FPR correction
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size	uncorr_size	level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
64	1.00	1.00	1.00	1.00	3.1MiB	2.2MiB	:0	:1	:3.1MiB	:2.2MiB	:64	:64	:100.00	:15	:6.40	:4.20	:3.06
128	1.22	1.22	1.40	1.72	4.3MiB	2.2MiB	:0	:1	:4.3MiB	:2.2MiB	:128	:128	:100.00	:31	:12.80	:6.10	:4.39
# Best t_max (regarding expected query runtime): 64
)expected_cout";

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file, expected_cout) << written_file;
}
