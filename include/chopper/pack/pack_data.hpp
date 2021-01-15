#pragma once

#include<string>
#include<vector>

#include <chopper/pack/previous_level.hpp>

struct pack_data
{
    std::vector<std::string> filenames;
    std::vector<size_t> kmer_counts;
    std::vector<std::vector<std::string>> extra_information;

    //!\brief A reference to the output stream to cache the results to.
    std::stringstream * output_buffer{nullptr};
    //!\brief A reference to the stream to cache the header to.
    std::stringstream * header_buffer{nullptr};
    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level previous{};
};
