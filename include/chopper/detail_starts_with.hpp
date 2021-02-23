#pragma once

#include <string>
#include <string_view>

inline bool starts_with(std::string const & target, std::string_view const & query)
{
    size_t index{};
    while (index < target.size() && index < query.size() && target[index] == query[index])
        ++index;
    return index == query.size();
}
