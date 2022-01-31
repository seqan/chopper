#pragma once

#include <cmath>
#include <cstddef>
#include <string>

namespace chopper
{

/*!\brief Returns the smallest natural number that is greater or equal to `value` and a multiplicative of 64.
* \param[in] value The Input value that is smaller or equal to the return value.
*/
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value)
{
    return ((value + 63) >> 6) << 6;
}

//!\brief Round bytes to the appropriate unit and convert to string with unit.
[[nodiscard]] inline std::string byte_size_to_formatted_str(size_t bytes)
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
            result += "Bytes";
            break;
        case 1:
            result += "KiB";
            break;
        case 2:
            result += "MiB";
            break;
        case 3:
            result += "GiB";
            break;
    }

    return result;
}

} // namespace chopper
