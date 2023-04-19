#include <gtest/gtest.h>

#include <iostream>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/compute_fp_correction.hpp>
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
    size_t const cardinality = 100u;
    size_t const top_level_num_contained_user_bins = 2u;
    size_t const lower_level_split_bin_span = 1u;

    chopper::configuration config{}; // default config
    config.tmax = 64u;
    config.disable_estimate_union = true; /* also disable rearrangement */
    chopper::data_store data{};
    data.fp_correction = chopper::layout::compute_fp_correction(config.false_positive_rate,
                                                                config.num_hash_functions,
                                                                lower_level_split_bin_span);

    std::vector<std::string> filenames{"s1", "s2"};
    std::vector<chopper::sketch::hyperloglog> sketches{{}, {}};
    std::vector<size_t> kmer_counts{50, 50};

    chopper::layout::hibf_statistics stats(config, data.fp_correction, sketches, kmer_counts);

    for (size_t i = 0; i < num_top_level_bins; ++i)
    {
        chopper::layout::hibf_statistics::bin & bin =
            stats.top_level_ibf.bins.emplace_back(chopper::layout::hibf_statistics::bin_kind::merged,
                                                  cardinality,
                                                  top_level_num_contained_user_bins,
                                                  1u, // merged bin always is a single technical bin
                                                  std::vector<size_t>{0, 1}
            );

        for (size_t j = 0; j < top_level_num_contained_user_bins; ++j)
        {
            bin.child_level.bins.emplace_back(chopper::layout::hibf_statistics::bin_kind::split,
                                              cardinality,
                                              1u, // split bin always contains only a single user bin
                                              lower_level_split_bin_span,
                                              std::vector<size_t>{j % 2});
        }
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
64	1.00	16.00	1.00	16.00	1.2KiB	1.2KiB	:0:1	:1:4	:395Bytes:790Bytes	:395Bytes:790Bytes	:4:8	:4:2	:0.00:100.00	:-:1	:-:1.00	:-:1.00	:-:1.00
)expected_cout";

    EXPECT_EQ(summary, expected_cout);
}

TEST(execute_test, chopper_layout_statistics)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    std::vector<std::string> many_filenames;
    std::vector<size_t> many_kmer_counts;

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    for (size_t i{0}; i < 96u; ++i)
    {
        many_filenames.push_back(seqan3::detail::to_string("seq", i));
        many_kmer_counts.push_back(100 * ((i + 20) / 20));
    }

    chopper::configuration config{.data_file = "not needed",
                                  .disable_sketch_output = true,
                                  .output_filename = layout_file.c_str(),
                                  .tmax = 64,
                                  .disable_estimate_union = true /* also disable rearrangement */,
                                  .output_verbose_statistics = true};

    chopper::layout::layout hibf_layout{};

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = many_kmer_counts};

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    chopper::layout::execute(config, many_filenames, data);
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
64	1.00	1.26	1.00	1.26	73.3KiB	44.9KiB	:0:1	:1:12	:37.0KiB:36.2KiB	:37.0KiB:7.8KiB	:64:768	:64:64	:81.25:100.00	:1:32	:1.00:17.45	:1.00:6.20	:1.00:4.87
)expected_cout";

    EXPECT_EQ(layout_result_stdout, expected_cout) << layout_result_stdout;
    EXPECT_EQ(layout_result_stderr, std::string{});
}

TEST(execute_test, chopper_layout_statistics_determine_best_bins)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};
    std::filesystem::path const stats_file{binning_filename.string() + ".stats"};

    chopper::configuration config{.data_file = "not needed",
                                  .disable_sketch_output = true,
                                  .output_filename = binning_filename.c_str(),
                                  .tmax = 128,
                                  .disable_estimate_union = true /* also disable rearrangement */,
                                  .determine_best_tmax = true,
                                  .force_all_binnings = true,
                                  .output_verbose_statistics = true};

    chopper::layout::layout hibf_layout{};
    std::vector<std::string> filenames{"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7", "seq8", "seq9"};

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .hibf_layout = &hibf_layout,
                             .kmer_counts = {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000}};

    chopper::layout::execute(config, filenames, data);

    std::string expected_cout =
        R"expected_cout(## ### Parameters ###
## number of user bins = 10
## number of hash functions = 2
## false positive rate = 0.05
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
64	1.00	1.00	1.00	1.00	1.6MiB	1.2MiB	:0	:1	:1.6MiB	:1.2MiB	:64	:64	:100.00	:16	:6.40	:4.35	:3.06
128	1.22	1.22	1.40	1.71	2.3MiB	1.2MiB	:0	:1	:2.3MiB	:1.2MiB	:128	:128	:100.00	:33	:12.80	:6.29	:4.38
# Best t_max (regarding expected query runtime): 64
)expected_cout";

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file, expected_cout) << written_file;
}
