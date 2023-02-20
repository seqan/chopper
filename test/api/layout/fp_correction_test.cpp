#include <gtest/gtest.h>

#include "../api_test.hpp"
#include <chopper/data_store.hpp>

TEST(fp_correction_test, one_bin)
{
    chopper::data_store data;
    data.compute_fp_correction(0.05, 2u, 8u); //fpr=0.05, #hash=2, t_max=8

    std::vector<size_t> const values{9123, 123, 12, 87123, 8123, 4660};

    // Splitting into 1 bin, i.e. not splitting, should not change the bin size.
    for (size_t const value : values)
        EXPECT_EQ(value, value * data.fp_correction[1]);
}

TEST(fp_correction_test, example_split)
{
    chopper::data_store data;
    data.compute_fp_correction(0.01, 5u, 256u); //fpr=0.01, #hash=5, t_max=256

    double const abs_error = 0.00001;
    EXPECT_NEAR(data.fp_correction[1], 1.0, abs_error);
    EXPECT_NEAR(data.fp_correction[2], 1.192316, abs_error);
    EXPECT_NEAR(data.fp_correction[4], 1.412390, abs_error);
    EXPECT_NEAR(data.fp_correction[8], 1.664459, abs_error);
    EXPECT_NEAR(data.fp_correction[16], 1.953384, abs_error);
    EXPECT_NEAR(data.fp_correction[32], 2.284738, abs_error);
    EXPECT_NEAR(data.fp_correction[64], 2.664909, abs_error);
    EXPECT_NEAR(data.fp_correction[128], 3.101225, abs_error);
    EXPECT_NEAR(data.fp_correction[256], 3.602093, abs_error);

    ASSERT_EQ(data.fp_correction.size(), 257u);
    for (size_t i{1u}; i < 256u; ++i)
        ASSERT_LE(data.fp_correction[i], data.fp_correction[i + 1u]);
}
