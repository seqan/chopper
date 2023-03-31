#pragma once

#include <fstream>

#include <robin_hood.h>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>

namespace chopper::sketch
{

inline auto read_data_file_with_clustering(configuration const & config)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters; // result

    std::ifstream fin{config.data_file.string()};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open data file " + config.data_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        auto column_begin = line.begin();
        auto column_end = line.begin();

        // first column is always the filename
        while (column_end != line.end() && *column_end != '\t')
            ++column_end; // advance to next tab
        std::string filename{column_begin, column_end};

        // TODO check file extension

        if (config.column_index_to_cluster > 1)
        {
            // advance to column of interest
            for (size_t i = 1; i < config.column_index_to_cluster; ++i)
            {
                while (column_begin != line.end() && *column_begin != '\t')
                    ++column_begin;
                ++column_begin;
            }
            column_end = column_begin;
            while (column_end != line.end() && *column_end != '\t')
                ++column_end;

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

inline void read_data_file(configuration const & config, data_store & data)
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
            data.filenames.push_back(filename);
            data.extra_information_strings.push_back(std::string{});
        }
        else
        {
            std::string const filename{line.begin(), line.begin() + tab_pos};
            std::string const extra_info{line.begin() + tab_pos + 1, line.end()};
            data.filenames.push_back(filename);
            data.extra_information_strings.push_back(extra_info);
        }
    }
}

} // namespace chopper::sketch
