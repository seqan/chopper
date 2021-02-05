#pragma once

#include <chopper/search/technical_binning_directory.hpp>
#include <chopper/detail_hibf_user_bins.hpp>

struct search_data
{
    std::vector<seqan3::technical_binning_directory<>> hibf;
    std::vector<std::vector<int64_t>> hibf_bin_levels;
    hibf_user_bins user_bins;
};
