// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <chopper/configuration.hpp>

namespace chopper::sketch
{

//!\brief Checks the `filenames` for consistent files, either precomputed or sequence files.
void check_filenames(std::vector<std::string> const & filenames, configuration & config);

} // namespace chopper::sketch
