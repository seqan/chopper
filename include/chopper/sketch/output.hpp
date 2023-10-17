// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::sketch
{

void write_count_file_line(std::pair<std::string, std::vector<std::string>> const & cluster,
                           uint64_t const weight,
                           std::ofstream & fout);

void write_sketch_file(std::string const & filename,
                       seqan::hibf::sketch::hyperloglog const & sketch,
                       configuration const & config);

} // namespace chopper::sketch
