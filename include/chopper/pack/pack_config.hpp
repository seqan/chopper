#pragma once

#include <seqan3/std/filesystem>

struct pack_config
{
    std::filesystem::path data_file;
    std::filesystem::path output_filename{"binning.out"};
    uint16_t bins{64};
    double alpha{10};
    int aggregate_by_column{-1};
};
