#pragma once

struct build_config
{
    std::string chopper_split_filename{};
    std::string chopper_pack_filename{};
    std::string binning_filename{};
    std::string output_filename{"out.chopper"};
    uint8_t threads{1};
    bool verbose{false};

    size_t hash_funs{2};
    double FPR{0.0001};

    uint8_t k{25};
    uint16_t overlap{200}; // overlap when inserting sequence regions into the IBF
};
