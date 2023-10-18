// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::sketch
{

void read_hll_files_into(std::filesystem::path const & hll_dir,
                         std::vector<std::string> const & target_filenames,
                         std::vector<seqan::hibf::sketch::hyperloglog> & target);

} // namespace chopper::sketch
