#pragma once

#include <chopper/build/build_config.hpp>

inline size_t compute_bin_size(build_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return std::ceil( - static_cast<double>(number_of_kmers_to_be_stored * config.hash_funs) /
                     std::log(1 - std::exp(std::log(config.FPR) / config.hash_funs)));
}


// -NUM_ELEM*HASHES
// ----------------------  = SIZE
// LN(1-FPR^(1/HASHES))

// -NUM_ELEMS*HASHES
// -----------------------
// LN(1 - e^(LN(FPR) / HASHES) )
