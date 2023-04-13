#pragma once

#include <fstream>

#include <robin_hood.h>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>

namespace chopper::sketch
{

inline void read_data_file(configuration const & config, std::vector<std::string> & filenames)
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
