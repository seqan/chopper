#pragma once

#include <seqan3/std/filesystem>

struct split_config
{
    std::string data_filename;
    std::vector<std::string> seqfiles;
    std::string output_graph_file{"graph.dot"}; // default

    uint8_t kmer_size{25};
    uint16_t window_size{100};

    // traverse config
    std::filesystem::path out_path{"/tmp/traverse_graph.out"};
    int16_t bins{64};
    bool write_out_graph{false};
    bool write_out_weights{false};
    // For a merged low level IBF several splittings will come into the same file.
    // So the traversal output needs to be appended and the bin_index adjusted.
    size_t bin_index_offset{0};
};
