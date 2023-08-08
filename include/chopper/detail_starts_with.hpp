// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <string_view>

namespace chopper::detail
{

inline bool starts_with(std::string const & target, std::string_view const & query)
{
    size_t index{};
    while (index < target.size() && index < query.size() && target[index] == query[index])
        ++index;
    return index == query.size();
}

} // namespace chopper::detail
