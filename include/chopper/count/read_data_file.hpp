#pragma once

#include <fstream>

#include <robin_hood.h>

#include <chopper/configuration.hpp>

namespace chopper::count
{

inline auto read_data_file(configuration const & config)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters; // result

    std::ifstream fin{config.data_file};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open file " + config.data_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        auto column_begin = line.begin();
        auto column_end = line.begin();

        // first column is always the filename
        while (column_end != line.end() && *column_end != '\t') ++column_end; // advance to next tab
        std::string filename{column_begin, column_end};

        // TODO check file extension

        if (config.column_index_to_cluster > 1)
        {
            // advance to column of interest
            for (size_t i = 1; i < config.column_index_to_cluster; ++i)
            {
                while (column_begin != line.end() && *column_begin != '\t') ++column_begin;
                ++column_begin;
            }
            column_end = column_begin;
            while (column_end != line.end() && *column_end != '\t') ++column_end;

            std::string key(column_begin, column_end);

            filename_clusters[key].push_back(filename);
        }
        else
        {
            // use filename itself as key
            filename_clusters[filename].push_back(filename);
        }

    }

    return filename_clusters;
}

} // namespace chopper::count
