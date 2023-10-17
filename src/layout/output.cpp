// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <iostream>

#include <chopper/layout/output.hpp>
#include <chopper/prefixes.hpp>

namespace chopper::layout
{

void write_user_bins_to(std::vector<std::string> const & filenames, std::ostream & stream)
{
    stream << chopper::prefix::meta_chopper_user_bins_start << '\n';
    size_t counter{};
    for (auto const & filename : filenames)
        stream << seqan::hibf::prefix::meta_header << counter++ << ' ' << filename << '\n';
    stream << chopper::prefix::meta_chopper_user_bins_end << '\n';
}

} // namespace chopper::layout
