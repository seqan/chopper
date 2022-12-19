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
