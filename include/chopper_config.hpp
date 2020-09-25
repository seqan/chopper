#pragma once

struct chopper_config
{
    std::vector<std::string> seqfiles;
    std::string output_graph_file{"graph.dot"}; // default

    uint8_t kmer_size{25};
    uint16_t window_size{100};
};
