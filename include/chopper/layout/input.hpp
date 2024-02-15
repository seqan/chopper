// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <iosfwd>
#include <string>
#include <tuple>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>

namespace chopper::layout
{

std::vector<std::vector<std::string>> read_filenames_from(std::istream & stream);

std::tuple<std::vector<std::vector<std::string>>, configuration, std::vector<seqan::hibf::layout::layout>>
read_layouts_file(std::istream & stream);

} // namespace chopper::layout
