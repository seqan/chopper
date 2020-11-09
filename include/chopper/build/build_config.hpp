#pragma once

struct build_config
{
    std::string traversal_path_prefix{};
    std::string binning_filename{};
    std::string output_prefix{"./"};

    uint8_t k{25};
    uint8_t overlap{200}; // overlap when inserting sequence regions into the IBF
};
