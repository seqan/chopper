#pragma once

#include <seqan3/std/charconv>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <chopper/configuration.hpp>
#include <chopper/layout/data_store.hpp>

namespace chopper::layout
{

inline auto read_filename_data_file(data_store & data, configuration const & config)
{
    std::ifstream file_in{config.count_filename};

    if (!file_in.good())
        throw std::runtime_error{"[CHOPPER LAYOUT ERROR] Could not open file " + config.count_filename.string()};

    std::string line;
    while (std::getline(file_in, line) && line[0] == '#'); // skip comments

    do
    {
        // read filename
        char const * buffer = line.c_str();
        auto ptr = &buffer[0];
        auto const buffer_end = ptr + line.size();

        if (line.empty())
            continue;

        while (ptr != buffer_end && *ptr != '\t') ++ptr;
        data.filenames.push_back(std::string(&buffer[0], ptr));

        if (ptr == buffer_end) // only file info, no kmer info
            throw std::runtime_error{"[CHOPPER LAYOUT ERROR] Your file only contains sequence names but no kmer counts."
                                     "Offending line: '" + line + "'."};

        // read kmer_count
        ++ptr; // skip tab
        size_t tmp;
        auto res = std::from_chars(ptr, buffer_end, tmp);
        data.kmer_counts.push_back(tmp);
        ptr = res.ptr;
        data.extra_information.push_back(std::vector<std::string>{});

        while (ptr != buffer_end) // even more information to come
        {
            ++ptr; // skip tab
            auto start = ptr;
            while (ptr != buffer_end && *ptr != '\t') ++ptr;
            data.extra_information.back().push_back(std::string(start, ptr));
        }
    }
    while (std::getline(file_in, line));
}

} // namespace chopper::layout
