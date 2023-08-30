// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <string_view>

#include <hibf/layout/prefixes.hpp>

namespace chopper::prefix
{

constexpr std::string_view meta_chopper_user_bins_start{"@CHOPPER_USER_BINS"};
static_assert(meta_chopper_user_bins_start.starts_with(seqan::hibf::prefix::meta_header));

constexpr std::string_view meta_chopper_user_bins_end{"@CHOPPER_USER_BINS_END"};
static_assert(meta_chopper_user_bins_end.starts_with(seqan::hibf::prefix::meta_header));

constexpr std::string_view meta_chopper_config_start{"@CHOPPER_CONFIG"};
static_assert(meta_chopper_config_start.starts_with(seqan::hibf::prefix::meta_header));

constexpr std::string_view meta_chopper_config_end{"@CHOPPER_CONFIG_END"};
static_assert(meta_chopper_config_end.starts_with(seqan::hibf::prefix::meta_header));

} // namespace chopper::prefix
