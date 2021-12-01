#pragma once

#include <array>
#include <seqan3/std/bit>
#include <cassert>
#include <map>
#include <stdexcept>

namespace chopper::layout
{

class ibf_query_cost
{
public:
    constexpr static inline const size_t maximum_t_max{65536};

    ibf_query_cost() = default;
    ibf_query_cost(ibf_query_cost const &) = default;
    ibf_query_cost & operator=(ibf_query_cost const &) = default;
    ibf_query_cost(ibf_query_cost &&) = default;
    ibf_query_cost & operator=(ibf_query_cost &&) = default;
    ~ibf_query_cost() = default;

    static double exact(size_t const t_max, double const fpr)
    {
        auto it = find_closest_fpr(fpr);

        if (contains(t_max))
            return it->second[position(t_max)];
        else
            throw std::invalid_argument("No exact data available for this t_max.");
    }

    static double interpolated(size_t const t_max, double const fpr)
    {
        auto it = find_closest_fpr(fpr);

        if (t_max <= 64u)
        {
            return it->second[0];
        }
        else if (t_max > maximum_t_max)
        {
            throw std::invalid_argument("No data availabe for a t_max this large.");
        }
        else if (contains(t_max))
        {
            return it->second[position(t_max)];
        }
        else
        {
            size_t const upper_bound{std::bit_ceil(t_max)};
            size_t const lower_bound{upper_bound >> 1};
            double const upper_value{it->second[position(upper_bound)]};
            double const lower_value{it->second[position(lower_bound)]};

            double const interpolated_value{lower_value + (upper_value - lower_value) * (t_max - lower_bound)
                                                                                      / lower_bound};
            assert(interpolated_value <= upper_value);
            return interpolated_value;
        }
    }

private:

    // dummy table to see if refactoring works before entering new values
    static inline const std::map<double, std::array<double, 11>> cost_factors
    {   /* FPR, cost factors relative to a 64 IBF */
        {0.0001, {1.0, 0.87, 1.07, 1.58, 1.86, 2.3, 3.56, 4.78, 6.89, 11.61, 22.7}},
        {0.0005, {1.0, 1.2, 1.04, 1.39, 1.91, 2.87, 3.88, 7.28, 7.99, 14.3, 28.34}},
        {0.0025, {1.0, 0.9, 1.08, 1.29, 1.91, 2.46, 3.84, 5.49, 7.93, 14.17, 28.16}},
        {0.0125, {1.0, 1.2, 1.36, 2.01, 2.64, 3.37, 5.98, 9.92, 16.29, 30.17, 59.08}},
        {0.0625, {1.0, 0.96, 1.2, 1.71, 2.27, 3.54, 6.42, 10.69, 19.31, 37.12, 72.5}},
        {0.3125, {1.0, 1.43, 1.83, 2.92, 4.68, 8.25, 16.58, 27.7, 52.3, 102.65, 202.55}}
    };

    static std::map<double, std::array<double, 11>>::const_iterator find_closest_fpr(double const fpr)
    {
        if (auto it = cost_factors.find(fpr); it != cost_factors.end()) // fpr is found exaclty in map
            return it;

        // otherwise search for the closest one in the map
        auto lower_it = cost_factors.lower_bound(fpr);
        auto upper_it = cost_factors.upper_bound(fpr);

        assert(lower_it != cost_factors.end() || upper_it != cost_factors.end());

        if (lower_it == cost_factors.end())
            return upper_it;

        if (upper_it == cost_factors.end())
            return lower_it;

        if (std::abs(lower_it->first - fpr) < std::abs(upper_it->first - fpr))
            return lower_it;
        else
            return upper_it;
    }

    constexpr static bool contains(size_t const value)
    {
        bool const is_power_of_two{std::has_single_bit(value)};
        int const trailing_zeros{std::countr_zero(value)};
        return is_power_of_two && trailing_zeros >= 6 && trailing_zeros <= 16;
    }

    constexpr static size_t position(size_t const value)
    {
        assert(contains(value));
        return std::countr_zero(value) - 6;
    }
};

} // namespace chopper::layout
