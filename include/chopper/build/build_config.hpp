#pragma once

struct build_config
{
    std::string traversal_path_prefix{};
    std::string binning_filename{};
    std::string output_prefix{"./"};

    size_t hash_funs{2};
    double FPR{0.0001};

    uint8_t k{25};
    uint16_t overlap{200}; // overlap when inserting sequence regions into the IBF
};
