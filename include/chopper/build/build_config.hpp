#pragma once

struct build_config
{
    std::string traversal_path_prefix{};
    std::string binning_filename{};
    std::string output_prefix{"./"};

    size_t high_level_ibf_num_technical_bins{};

    uint8_t k{25};
};
