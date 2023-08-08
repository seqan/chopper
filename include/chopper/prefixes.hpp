// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <string_view>

namespace chopper::prefix
{

constexpr std::string_view chopper{"chopper"};

constexpr std::string_view header{"#"};

constexpr std::string_view header_config{"#"};

constexpr std::string_view high_level{"HIGH_LEVEL_IBF"};

constexpr std::string_view first_header_line{"#HIGH_LEVEL_IBF"};
static_assert(first_header_line.starts_with(header));
static_assert(first_header_line.ends_with(high_level));

constexpr std::string_view merged_bin{"MERGED_BIN"};

constexpr std::string_view split_bin{"SPLIT_BIN"};

} // namespace chopper::prefix
