// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <charconv>    // for from_chars, from_chars_result
#include <cstddef>     // for size_t
#include <ranges>      // for all_t, operator==, single_view, operator|, _Partial
#include <string>      // for string
#include <string_view> // for string_view
#include <utility>     // for move
#include <vector>      // for vector

#include <hibf/detail/layout/layout.hpp> // for layout

namespace chopper::stats
{

inline seqan::hibf::layout::layout::user_bin
parse_layout_line(std::string const & current_line, std::vector<std::vector<std::string>> & user_bin_filenames)
{
    seqan::hibf::layout::layout::user_bin result{};

    // initialize parsing
    std::string_view const buffer{current_line};
    auto const buffer_end{buffer.end()};
    auto field_end = buffer.begin();
    while (field_end != buffer_end && *field_end != '\t')
        ++field_end;

    // parse filenames
    std::string_view const filenames_str{buffer.begin(), field_end};
    std::vector<std::string> filename_list;
    for (auto const && filename : filenames_str | std::views::split(';'))
    {
        auto const common_view = filename | std::views::common;
        filename_list.emplace_back(common_view.begin(), common_view.end());
    }

    // update input idx and append filesnames to build_data
    result.idx = user_bin_filenames.size();
    user_bin_filenames.push_back(std::move(filename_list));

    size_t tmp{};

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.previous_TB_indices.push_back(tmp);
    }
    while (field_end != buffer_end && *field_end != '\t');

    result.storage_TB_id = result.previous_TB_indices.back();
    result.previous_TB_indices.pop_back();

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.number_of_technical_bins = tmp; // only the last number really counts
    }
    while (field_end != buffer_end && *field_end != '\t');

    return result;
}

} // namespace chopper::stats
