#pragma once

#include <chopper/build/build_config.hpp>

size_t compute_bin_size(build_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return (size_t)std::ceil((-(double)(config.hash_funs * number_of_kmers_to_be_stored)) /
                   (std::log(1 - std::pow(10.0, std::log10(config.FPR) / config.hash_funs))));
}
