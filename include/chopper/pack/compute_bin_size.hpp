#pragma once

#include <chopper/pack/pack_config.hpp>

inline size_t compute_bin_size(pack_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return std::ceil( - static_cast<double>(number_of_kmers_to_be_stored * config.num_hash_functions) /
                     std::log(1 - std::exp(std::log(config.fp_rate) / config.num_hash_functions)));
}


// -NUM_ELEM*HASHES
// ----------------------  = SIZE
// LN(1-FPR^(1/HASHES))

// -NUM_ELEMS*HASHES
// -----------------------
// LN(1 - e^(LN(FPR) / HASHES) )