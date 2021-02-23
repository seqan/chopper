#pragma once

#include <utility>

struct pair_hash
{
    std::size_t operator () (std::pair<int32_t, uint32_t> const & pair) const
    {
        return (static_cast<size_t>(pair.first) << 32) | static_cast<size_t>(pair.second);
    }
};
