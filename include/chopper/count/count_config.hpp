#pragma once

#include <seqan3/std/filesystem>
#include <thread>

struct count_config
{
    std::filesystem::path data_file{};
    size_t column_index_to_cluster{1};
    size_t num_threads{std::thread::hardware_concurrency()};
    uint8_t k{25};
    unsigned int w{500};
    bool disable_minimizers{false};
};
