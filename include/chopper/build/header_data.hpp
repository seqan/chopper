#pragma once

#include <chopper/detail_parse_binning_line.hpp>

struct header_data
{
    size_t hibf_num_technical_bins{};
    std::string hibf_max_bin_id{};
    data_file_record * hibf_max_record{nullptr};
};
