#pragma once

#include <string>

//!\brief Round bytes to the appropriate unit and convert to string with unit
std::string byte_size_to_formatted_str(size_t bytes)
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