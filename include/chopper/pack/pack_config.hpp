#pragma once

#include <seqan3/std/filesystem>

#include <chopper/pack/previous_level.hpp>

struct pack_config
{
    std::filesystem::path data_file;
    std::filesystem::path output_filename{"binning.out"};
    uint16_t bins{64};
    std::filesystem::path hll_dir{};
    size_t num_hash_functions{2};
    double fp_rate{0.05};
    double alpha{10};
    double max_ratio{0.5};
    size_t num_threads{1u};
    int aggregate_by_column{-1};
    bool union_estimate{false};
    bool rearrange_bins{false};
};
