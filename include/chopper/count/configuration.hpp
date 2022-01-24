#pragma once

#include <filesystem>
#include <thread>

namespace chopper::count
{

struct configuration
{
    std::filesystem::path data_file;
    std::string output_prefix{}; // given by the user
    std::filesystem::path count_filename{};   // set internally
    std::filesystem::path sketch_directory{}; // set internally
    size_t column_index_to_cluster{1u};
    size_t threads{1u};
    uint8_t k{19};
    uint8_t sketch_bits{12};
    bool disable_sketch_output{false};
    //!\brief Whether the input files are precomputed files instead of sequence files.
    bool precomputed_files{false};
};

} // namespace chopper::count
