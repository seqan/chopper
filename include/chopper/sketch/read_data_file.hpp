#pragma once

#include <fstream>

#include <robin_hood.h>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>

namespace chopper::sketch
{

inline void read_data_file(std::string const & input_filename, std::vector<std::string> & filenames)
{
    std::ifstream fin{input_filename};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open data file " + input_filename + " for reading."};

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
