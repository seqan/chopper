#include <gtest/gtest.h>

#include <chopper/layout/ibf_query_cost.hpp>

#include "../api_test.hpp"

TEST(ibf_query_cost_test, exact)
{
    double value{};
    for (size_t i{64}; i <= chopper::layout::ibf_query_cost::maximum_t_max; i *= 2)
    {
        double result = chopper::layout::ibf_query_cost::exact(i, 0.0125);
        EXPECT_GT(result, value);
        value = result;
    }

    EXPECT_NO_THROW(chopper::layout::ibf_query_cost::exact(chopper::layout::ibf_query_cost::maximum_t_max, 0.0125));
    ASSERT_EQ(chopper::layout::ibf_query_cost::exact(chopper::layout::ibf_query_cost::maximum_t_max, 0.0125), 90.42);
    EXPECT_THROW(chopper::layout::ibf_query_cost::exact(chopper::layout::ibf_query_cost::maximum_t_max + 1, 0.0125), std::invalid_argument);
}

TEST(ibf_query_cost_test, interpolated)
{
    for (size_t i{64}; i <= chopper::layout::ibf_query_cost::maximum_t_max; i *= 2)
        EXPECT_EQ(chopper::layout::ibf_query_cost::interpolated(i, 0.0125), chopper::layout::ibf_query_cost::exact(i, 0.0125));

    double value{};
    for (size_t i{67}; i < chopper::layout::ibf_query_cost::maximum_t_max; i *= 2)
    {
        double result = chopper::layout::ibf_query_cost::interpolated(i, 0.0125);
        EXPECT_GT(result, value);
        EXPECT_GT(result, chopper::layout::ibf_query_cost::exact(1ULL << (std::bit_width(i) - 1), 0.0125)); // std::bit_floor not in seqan3
        EXPECT_LT(result, chopper::layout::ibf_query_cost::exact(std::bit_ceil(i), 0.0125));
        value = result;
    }

    EXPECT_NO_THROW(chopper::layout::ibf_query_cost::interpolated(chopper::layout::ibf_query_cost::maximum_t_max, 0.0125));
    EXPECT_THROW(chopper::layout::ibf_query_cost::interpolated(chopper::layout::ibf_query_cost::maximum_t_max + 1, 0.0125), std::invalid_argument);
}
