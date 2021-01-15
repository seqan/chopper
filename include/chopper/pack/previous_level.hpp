#pragma once

#include <string>

//!\brief Information about the previous IBF level to be passed down to ensure correct output.
struct previous_level
{
    std::string bin_indices;
    std::string num_of_bins;
    std::string estimated_sizes;

    bool empty() const
    {
        assert(bin_indices.empty() == num_of_bins.empty());
        assert(bin_indices.empty() == estimated_sizes.empty());
        assert(estimated_sizes.empty() == num_of_bins.empty());
        return bin_indices.empty();
    }
};
