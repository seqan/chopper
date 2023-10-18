// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/sketch/read_data_file.hpp>

namespace chopper::sketch
{

void read_data_file(configuration const & config, std::vector<std::string> & filenames)
{
    std::ifstream fin{config.data_file.string()};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open data file " + config.data_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        auto tab_pos = line.find('\t');

        if (tab_pos == std::string::npos)
        {
            std::string const filename{line.begin(), line.end()};
            filenames.push_back(filename);
        }
        else
        {
            std::string const filename{line.begin(), line.begin() + tab_pos};
            filenames.push_back(filename);
        }
    }
}

} // namespace chopper::sketch
