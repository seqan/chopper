#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>

namespace chopper
{

/*!\brief Returns the smallest natural number that is greater or equal to `value` and a multiplicative of 64.
* \param[in] value The Input value that is smaller or equal to the return value.
*/
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value) noexcept
{
    return ((value + 63) >> 6) << 6;
}

//!\brief Round bytes to the appropriate unit and convert to string with unit.
[[nodiscard]] inline std::string byte_size_to_formatted_str(size_t const bytes)
{
    size_t iterations{};
    size_t integer{bytes};

    while (integer >> 10u && iterations < 6u)
    {
        integer >>= 10u;
        ++iterations;
    }

    // While this is a bit more involved, we can avoid using floating point numbers.
    auto first_decimal_position = [&]()
    {
        assert(iterations > 0u);
        size_t decimal{bytes};
        decimal -= integer << (iterations * 10u);             // Substract bytes represented by integer, e.g. -5GiB
        decimal >>= (iterations - 1u) * 10u;                  // Shift to next smallest unit, e.g. 800MiB
        decimal = decimal * 1000u / 1024u;                    // Account for using decimal system, i.e. 800MiB != 0.8GiB
        size_t const diff{decimal - (decimal / 100u) * 100u}; // We want to round up to 1 decimal position
        uint32_t const round_up{diff >= 50u};
        decimal += round_up * 100u - diff;
        decimal /= 100u;
        return decimal;
    };

    auto formatted_string = [&]()
    {
        static constexpr int8_t int_to_char_offset{'0'}; // int 0 as char: char{0 + 48} = '0'
        size_t const decimal = iterations ? first_decimal_position() : 0u;
        assert(decimal <= 10u);

        if (!iterations) // No decimals for Bytes
            return std::to_string(integer);
        else if (decimal < 10u) // No need to round integer part
            return std::to_string(integer) + '.' + static_cast<char>(decimal + int_to_char_offset);
        else // Round integer part, e.g., 5.99 MiB should report 6.0 MiB
        {
            ++integer;
            // Check whether rounding results in a change of unit, e.g. 1023.99MiB to 1.0GiB
            if (integer >> 10u)
            {
                ++iterations;
                integer >>= 10u;
            }
            return std::to_string(integer) + ".0";
        }
    };

    std::string result{formatted_string()};
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
    case 4:
        result += "TiB";
        break;
    case 5:
        result += "PiB";
        break;
    default:
        result += "EiB";
        break;
    }

    return result;
}

} // namespace chopper
