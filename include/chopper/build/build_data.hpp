#pragma once

#include <chopper/build/batch.hpp>
#include <chopper/build/region.hpp>
#include <chopper/detail_parse_binning_line.hpp>

struct build_data
{
    size_t hibf_num_technical_bins{};
    size_t num_libfs{};
    std::string hibf_max_bin_id{};
    size_t hibf_max_bin{};
    data_file_record * hibf_max_record{nullptr};
    batch * hibf_max_batch_record{nullptr};

    std::unordered_map<size_t, size_t> merged_max_bin_map{};
    std::unordered_map<std::string, size_t> merged_bin_map{};
    std::unordered_map<std::string, std::vector<region>> region_map{};
};
