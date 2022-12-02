#pragma once

#include <cassert>
#include <string>

namespace chopper::layout
{

//!\brief Information about the previous IBF level to be passed down to ensure correct output.
struct previous_level
{
    std::string bin_indices;
    std::string num_of_bins;
    std::string estimated_sizes;
    std::string optimal_score;
    std::string correction;
    std::string tmax;

    bool empty() const
    {
        assert(((((bin_indices.empty() == num_of_bins.empty()) == estimated_sizes.empty()) == optimal_score.empty())
                == correction.empty())
               == tmax.empty());
        return bin_indices.empty();
    }
};

} // namespace chopper::layout
