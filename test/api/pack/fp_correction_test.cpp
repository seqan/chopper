#include <gtest/gtest.h>

#include <chopper/pack/pack_data.hpp>

#include "../api_test.hpp"

TEST(fp_correction_test, one_bin)
{
    pack_data data;
    data.compute_fp_correction(0.05, 2u, 8u); //fpr=0.05, #hash=2, t_max=8

    std::vector<size_t> const values{9123, 123, 12, 87123, 8123, 4660};

    for (size_t const value : values)
        EXPECT_EQ(value, value * data.fp_correction[1]);
}
