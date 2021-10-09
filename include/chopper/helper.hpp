#pragma once

#include <cstddef>
#include <string>

#include <chopper/pack/pack_config.hpp>

/*!\brief Returns the smallest natural number that is greater or equal to `value` and a multiplicative of 64.
* \param[in] value The Input value that is smaller or equal to the return value.
*/
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value)
{
    return ((value + 63) >> 6) << 6;
}

//!\brief Round bytes to the appropriate unit and convert to string with unit
[[nodiscard]] std::string byte_size_to_formatted_str(size_t bytes)
{
    size_t iterations{};
    while (bytes >> 10 && iterations < 3)
    {
        bytes >>= 10;
        ++iterations;
    }

    std::string result{std::to_string(bytes)};
    switch (iterations)
    {
        case 0:
            result += " Bytes";
            break;
        case 1:
            result += " KiB";
            break;
        case 2:
            result += " MiB";
            break;
        case 3:
            result += " GiB";
            break;
    }

    return result;
}

// -NUM_ELEM*HASHES
// ----------------------  = SIZE
// LN(1-FPR^(1/HASHES))

// -NUM_ELEMS*HASHES
// -----------------------
// LN(1 - e^(LN(FPR) / HASHES) )
inline size_t compute_bin_size(pack_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return std::ceil( - static_cast<double>(number_of_kmers_to_be_stored * config.num_hash_functions) /
                     std::log(1 - std::exp(std::log(config.fp_rate) / config.num_hash_functions)));
}
