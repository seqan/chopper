// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <utility>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

seqan::hibf::layout::layout
determine_best_number_of_technical_bins(chopper::configuration & config,
                                        std::vector<size_t> const & kmer_counts,
                                        std::vector<seqan::hibf::sketch::hyperloglog> const & sketches);

}
