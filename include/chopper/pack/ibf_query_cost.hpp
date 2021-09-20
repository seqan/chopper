#pragma once

#include <array>
#include <seqan3/std/bit>
#include <cassert>
#include <stdexcept>

class ibf_query_cost
{
public:
    ibf_query_cost() = default;
    ibf_query_cost(ibf_query_cost const &) = default;
    ibf_query_cost & operator=(ibf_query_cost const &) = default;
    ibf_query_cost(ibf_query_cost &&) = default;
    ibf_query_cost & operator=(ibf_query_cost &&) = default;
    ~ibf_query_cost() = default;

    constexpr static double get_exact(size_t const t_max)
    {
        if (contains(t_max))
            return cost_factor[position(t_max)];
        else
            throw std::invalid_argument("No exact data available for this t_max.");
    }

    constexpr static double get_interpolated(size_t const t_max)
    {
        if (t_max <= 64u)
        {
            return cost_factor[0];
        }
        else if (t_max > maximum_t_max)
        {
            throw std::invalid_argument("No data availabe for a t_max this large.");
        }
        else if (contains(t_max))
        {
            return cost_factor[position(t_max)];
        }
        else
        {
            size_t const upper_bound{std::bit_ceil(t_max)};
            size_t const lower_bound{upper_bound >> 1};
            double const upper_value{cost_factor[position(upper_bound)]};
            double const lower_value{cost_factor[position(lower_bound)]};

            double const interpolated_value{lower_value + (upper_value - lower_value) * (t_max - lower_bound)
                                                                                      / lower_bound};
            assert(interpolated_value <= upper_value);
            return interpolated_value;
        }
    }

private:
    constexpr static inline const size_t maximum_t_max{65536};

    // (19,19) mean c_{T_max} (see paper)
    constexpr static inline const std::array<double, 11> cost_factor
    {
        //64  128   256   512  1024  2048  4096  8192   16384  32768  65536
        {1.0, 1.1, 1.32, 1.61, 2.69, 4.48, 7.53, 13.65, 23.86, 45.66, 90.42}
    };

    constexpr static bool contains(size_t const value)
    {
        bool const is_power_of_two{std::has_single_bit(value)};
        int const set_bit_index{std::countr_zero(value) + 1};
        return is_power_of_two && set_bit_index >= 6 && set_bit_index <= 16;
    }

    constexpr static size_t position(size_t const value)
    {
        assert(contains(value));
        return std::countr_zero(value) - 6;
    }
};
