#pragma once

#include <seqan3/std/filesystem>

struct pack_config
{
    std::filesystem::path data_file;
    uint16_t bins{64};
    int aggregate_by_column{-1};
};
