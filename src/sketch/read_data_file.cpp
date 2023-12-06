// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <filesystem>
#include <fstream>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/utility/range/to.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/read_data_file.hpp>

namespace chopper::sketch
{

void read_data_file(configuration const & config, std::vector<std::vector<std::string>> & filenames)
{
    std::ifstream fin{config.data_file.string()};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open data file " + config.data_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        std::vector<std::string> names;

        auto const tab_pos = line.find('\t');
        std::string_view const filename_sv{line.begin(),
                                           (tab_pos != std::string::npos) ? line.begin() + tab_pos : line.end()};

        // multiple filenames may be separated by ' '
        for (auto && name : std::views::split(filename_sv, ' '))
        {
            auto common_view = std::views::common(name);
            names.emplace_back(common_view.begin(), common_view.end());
        }

        filenames.push_back(std::move(names));
    }
}

} // namespace chopper::sketch
