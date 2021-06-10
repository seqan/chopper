#pragma once

struct search_config
{
    std::string chopper_index_filename{};
    std::string chopper_index_map_filename{};
    std::string query_filename{};
    std::string output_filename{"./search.out"};
    uint8_t k{25};
    uint8_t errors{0};
    uint8_t threads{1};
    bool verbose{false};
    bool write_time{false};
};
