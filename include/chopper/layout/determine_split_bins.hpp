// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::determine_split_bins.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <chopper/configuration.hpp>

namespace chopper::layout
{

std::pair<size_t, size_t> determine_split_bins(chopper::configuration const & config,
    std::vector<size_t> const & positions,
    std::vector<size_t> const & cardinalities,
    size_t const num_technical_bins,
    size_t const num_user_bins,
    std::vector<std::vector<size_t>> & partitions);

}