// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <fstream>   // for basic_istream, ifstream
#include <stdexcept> // for logic_error
#include <string>    // for char_traits, allocator, operator+, getline, string
#include <vector>    // for vector

#include <chopper/configuration.hpp>             // for config
#include <chopper/stats/parse_layout_header.hpp> // for parse_layout_header
#include <chopper/stats/parse_layout_line.hpp>   // for parse_layout_line

#include <hibf/detail/layout/layout.hpp> // for layout

namespace chopper::stats
{

inline hibf::layout::layout read_layout_file(configuration & hibf_config,
                                             std::vector<std::vector<std::string>> & filenames,
                                             std::string const & layout_filename)
{
    hibf::layout::layout hibf_layout{};

    std::ifstream layout_file{layout_filename};

    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + layout_filename + " for reading"}; // GCOVR_EXCL_LINE

    // parse header
    // -------------------------------------------------------------------------
    parse_layout_header(layout_file, hibf_config, hibf_layout);

    std::string current_line;
    while (std::getline(layout_file, current_line))
        hibf_layout.user_bins.emplace_back(parse_layout_line(current_line, filenames));

    return hibf_layout;
}

} // namespace chopper::stats
