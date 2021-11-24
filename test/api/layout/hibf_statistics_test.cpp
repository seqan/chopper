#include <gtest/gtest.h>

#include <iostream>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/configuration.hpp>
#include <chopper/layout/data_store.hpp>

#include "../api_test.hpp"

TEST(hibf_statistics, only_merged_on_top_level)
{
    // parameters for this test HIBF
    size_t const num_top_level_bins = 4u;
    size_t const cardinality = 100u;
    size_t const top_level_num_contained_user_bins = 2u;
    size_t const lower_level_split_bin_span = 1u;

    chopper::layout::configuration config; // default config
    chopper::layout::data_store data;
    data.compute_fp_correction(config.fp_rate, config.num_hash_functions, lower_level_split_bin_span);

    chopper::layout::hibf_statistics stats(config, data.fp_correction);

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

    stats.print_summary();
    std::cout.flush();

    std::string summary = testing::internal::GetCapturedStdout();
    std::string expected_cout =
R"expected_cout(level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
0	1	395 Bytes	395 Bytes	4	4	0.00	-	-	-	-
1	4	790 Bytes	790 Bytes	8	2	100.00	1	1.00	1.00	1.00
#Total HIBF size:1 KiB
#Total HIBF size no correction:1 KiB

)expected_cout";

    EXPECT_EQ(summary, expected_cout);
}
