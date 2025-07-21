// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout
{

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches);

} // namespace chopper::layout
