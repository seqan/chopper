// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <array>
#include <bit>
#include <cassert>
#include <map>
#include <stdexcept>

namespace chopper::layout
{

class ibf_query_cost
{
public:
    static inline constexpr const size_t maximum_t_max{65536};

    ibf_query_cost() = default;
    ibf_query_cost(ibf_query_cost const &) = default;
    ibf_query_cost & operator=(ibf_query_cost const &) = default;
    ibf_query_cost(ibf_query_cost &&) = default;
    ibf_query_cost & operator=(ibf_query_cost &&) = default;
    ~ibf_query_cost() = default;

    static double exact(size_t const t_max, double const fpr);

    static double interpolated(size_t const t_max, double const fpr);

private:
    /*!\brief The cost factor to penalize a search in an IBF with more then 64 bins.
     *
     * The table contains experimentally derived cost factors that were collected when measuring the runtime of an
     * (vanilla) IBF with 64/128/... technical bins given a specific FPR. Each measurement was normalized by the
     * runtime of an 64 bin IBF.
     *
     * Each run was conducted 5 times and the mean was taken over the measurements. Low FPR rates were observed to have
     * a rather high variance on the runs.
     *
     * See also `test/benchmark/benchmark_data/query_cost`.
     */
    static inline const std::map<double, std::array<double, 11>> cost_factors{
        /* FPR, cost factors relative to a 64 IBF */
        {0.0001, {1.0000, 1.0602, 1.3492, 1.3524, 1.5645, 1.9595, 3.4143, 5.4849, 6.8115, 10.9489, 19.8932}},
        {0.0005, {1.0000, 1.0534, 1.1068, 1.2821, 1.5151, 1.7112, 3.6442, 4.7700, 6.9978, 12.2086, 22.5374}},
        {0.0025, {1.0000, 1.0015, 1.0031, 1.0876, 1.4027, 1.6920, 3.3014, 4.8019, 7.6273, 13.5664, 24.1108}},
        {0.0125, {1.0000, 1.0071, 1.1713, 1.3430, 1.8335, 2.6955, 5.3925, 8.6168, 15.0510, 28.3340, 54.4134}},
        {0.0500, {1.0000, 1.2241, 1.3336, 1.6827, 2.4608, 3.7554, 7.3573, 12.4689, 23.2699, 45.0874, 86.5339}},
        {0.0625, {1.0000, 1.1011, 1.2670, 1.5964, 2.4030, 3.6996, 7.1772, 12.4852, 23.3882, 44.7427, 87.8259}},
        {0.3125, {1.0000, 1.2818, 1.5493, 2.2546, 3.7804, 6.5428, 12.9410, 24.4539, 47.6262, 93.4733, 185.1019}}};

    static std::map<double, std::array<double, 11>>::const_iterator find_closest_fpr(double const fpr);

    static constexpr bool contains(size_t const value)
    {
        bool const is_power_of_two{std::has_single_bit(value)};
        int const trailing_zeros{std::countr_zero(value)};
        return is_power_of_two && trailing_zeros >= 6 && trailing_zeros <= 16;
    }

    static constexpr size_t position(size_t const value)
    {
        assert(contains(value));
        return std::countr_zero(value) - 6;
    }
};

} // namespace chopper::layout
