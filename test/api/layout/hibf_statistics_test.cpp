#include <gtest/gtest.h>

#include <iostream>

#include <chopper/configuration.hpp>
#include <chopper/layout/data_store.hpp>
#include <chopper/layout/hibf_statistics.hpp>

#include "../api_test.hpp"

TEST(hibf_statistics, only_merged_on_top_level)
{
    // parameters for this test HIBF
    size_t const num_top_level_bins = 4u;
    size_t const cardinality = 100u;
    size_t const top_level_num_contained_user_bins = 2u;
    size_t const lower_level_split_bin_span = 1u;

    chopper::configuration config{}; // default config
    chopper::layout::data_store data{};
    data.compute_fp_correction(config.false_positive_rate, config.num_hash_functions, lower_level_split_bin_span);
    std::vector<size_t> kmer_counts{50, 50};

    chopper::layout::hibf_statistics stats(config, data.fp_correction, kmer_counts);

    for (size_t i = 0; i < num_top_level_bins; ++i)
    {
        chopper::layout::hibf_statistics::bin & bin = stats.top_level_ibf.bins.emplace_back(
            chopper::layout::hibf_statistics::bin_kind::merged,
            cardinality,
            top_level_num_contained_user_bins,
            1u // merged bin always is a single technical bin
        );

        for (size_t j = 0; j < top_level_num_contained_user_bins; ++j)
        {
            bin.child_level.bins.emplace_back(
                chopper::layout::hibf_statistics::bin_kind::split,
                cardinality,
                1u, // split bin always contains only a single user bin
                lower_level_split_bin_span
            );
        }
    }

    testing::internal::CaptureStdout();

    stats.print_header();
    size_t max_64{};
    stats.print_summary(max_64);
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
