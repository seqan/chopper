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

/*!\brief Reads a layout file: the file names of the user bins, the chopper configuration and the layout.
 * \param[in] stream The content of a layout file.
 * \returns The file names, the configuration and the layout.
 * \details
 * The layout is not validated. Validate it with seqan::hibf::layout::layout::validate and the returned
 * `configuration::hibf_config` before using it.
 */
std::tuple<std::vector<std::vector<std::string>>, configuration, seqan::hibf::layout::layout>
read_layout_file(std::istream & stream);

} // namespace chopper::layout
