// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <fstream>

#include <chopper/sketch/read_hll_files_into.hpp>

namespace chopper::sketch
{

void read_hll_files_into(std::filesystem::path const & hll_dir,
                         std::vector<std::string> const & target_filenames,
                         std::vector<seqan::hibf::sketch::hyperloglog> & target)
{
    assert(std::filesystem::exists(hll_dir) && !std::filesystem::is_empty(hll_dir)); // checked in chopper_layout

    target.reserve(target_filenames.size());

    try
    {
        for (auto const & filename : target_filenames)
        {
            std::filesystem::path path = hll_dir / std::filesystem::path(filename).stem();
            path += ".hll";
            std::ifstream hll_fin(path, std::ios::binary);

            if (!hll_fin.good())
                throw std::runtime_error{"Could not open file " + path.string()};

            // the sketch bits will be automatically read from the files
            target.emplace_back().load(hll_fin);
        }
    }
    catch (std::runtime_error const & err)
    {
        std::string const chopper_msg{"[CHOPPER LAYOUT ERROR] Something went wrong trying to read the HyperLogLog"
                                      " sketches from files:\n"};
        throw std::runtime_error{chopper_msg + err.what()};
    }
}

} // namespace chopper::sketch
