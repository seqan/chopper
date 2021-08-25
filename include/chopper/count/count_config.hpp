#pragma once

#include <seqan3/std/filesystem>
#include <thread>

struct count_config
{
    std::filesystem::path data_file{};
    std::filesystem::path output_filename{};
    std::filesystem::path hll_dir{};
    size_t column_index_to_cluster{1u};
    size_t num_threads{1u};
    uint8_t k{25};
    unsigned int w{500};
    uint8_t sketch_bits{12};
    bool disable_minimizers{false};
    bool exclusively_hlls{false};
};
