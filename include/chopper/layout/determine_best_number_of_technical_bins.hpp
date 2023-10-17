// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

std::pair<seqan::hibf::layout::layout, std::vector<seqan::hibf::sketch::hyperloglog>>
determine_best_number_of_technical_bins(chopper::configuration & config);

}
