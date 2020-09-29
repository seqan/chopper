#pragma once

#include<string>
#include<vector>

struct pack_data
{
    std::vector<std::string> filenames;
    std::vector<size_t> kmer_counts;
    std::vector<std::vector<std::string>> extra_information;
};
