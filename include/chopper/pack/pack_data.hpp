#pragma once

#include <seqan3/std/filesystem>

struct pack_data
{
    std::vector<std::filesystem::path> filenames;
    std::vector<size_t> kmer_counts;
    std::vector<std::vector<std::string>> extra_information;
};
