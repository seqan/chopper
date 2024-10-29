// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>

#include <chopper/layout/ibf_query_cost.hpp>

namespace chopper::layout
{

double ibf_query_cost::exact(size_t const t_max, double const fpr)
{
    auto it = find_closest_fpr(fpr);

    if (contains(t_max))
        return it->second[position(t_max)];
    else
        throw std::invalid_argument("No exact data available for this t_max.");
}

double ibf_query_cost::interpolated(size_t const t_max, double const fpr)
{
    auto it = find_closest_fpr(fpr);

    if (t_max <= 64u)
    {
        return it->second[0];
    }
    else if (t_max > maximum_t_max)
    {
        throw std::invalid_argument("No data available for a t_max this large.");
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

        double const interpolated_value{lower_value
                                        + (upper_value - lower_value) * (t_max - lower_bound) / lower_bound};
        assert(interpolated_value <= upper_value);
        return interpolated_value;
    }
}

std::map<double, std::array<double, 11>>::const_iterator ibf_query_cost::find_closest_fpr(double const fpr)
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

} // namespace chopper::layout
