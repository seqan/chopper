#include <gtest/gtest.h>

#include <chopper/pack/ibf_query_cost.hpp>

#include "../api_test.hpp"

TEST(ibf_query_cost_test, exact)
{
    double value{};
    for (size_t i{64}; i <= ibf_query_cost::maximum_t_max; i *= 2)
    {
        double result = ibf_query_cost::exact(i);
        EXPECT_GT(result, value);
        value = result;
    }

    EXPECT_NO_THROW(ibf_query_cost::exact(ibf_query_cost::maximum_t_max));
    ASSERT_EQ(ibf_query_cost::exact(ibf_query_cost::maximum_t_max), 90.42);
    EXPECT_THROW(ibf_query_cost::exact(ibf_query_cost::maximum_t_max + 1), std::invalid_argument);
}

TEST(ibf_query_cost_test, interpolated)
{
    for (size_t i{64}; i <= ibf_query_cost::maximum_t_max; i *= 2)
        EXPECT_EQ(ibf_query_cost::interpolated(i), ibf_query_cost::exact(i));

    double value{};
    for (size_t i{67}; i < ibf_query_cost::maximum_t_max; i *= 2)
    {
        double result = ibf_query_cost::interpolated(i);
        EXPECT_GT(result, value);
        EXPECT_GT(result, ibf_query_cost::exact(1ULL << (std::bit_width(i) - 1))); // std::bit_floor not in seqan3
        EXPECT_LT(result, ibf_query_cost::exact(std::bit_ceil(i)));
        value = result;
    }

    EXPECT_NO_THROW(ibf_query_cost::interpolated(ibf_query_cost::maximum_t_max));
    EXPECT_THROW(ibf_query_cost::interpolated(ibf_query_cost::maximum_t_max + 1), std::invalid_argument);
}
