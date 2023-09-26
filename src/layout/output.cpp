// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cstddef>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <chopper/layout/output.hpp>
#include <chopper/prefixes.hpp>

#include <hibf/layout/prefixes.hpp>

namespace chopper::layout
{

void write_user_bins_to(std::vector<std::vector<std::string>> const & filenames, std::ostream & stream)
{
    stream << chopper::prefix::meta_chopper_user_bins_start << '\n';
    size_t counter{};
    for (auto const & filenames_of_user_bin : filenames)
    {
        // the below will write lines like this:
        // @0 file1.fa file2.fa
        // @1 fileABC.fa
        stream << seqan::hibf::prefix::meta_header << counter++;
        for (std::string const & filename : filenames_of_user_bin)
            stream << ' ' << filename;
        stream << '\n';
    }
    stream << chopper::prefix::meta_chopper_user_bins_end << '\n';
}

} // namespace chopper::layout
