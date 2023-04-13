#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <chopper/next_multiple_of_64.hpp>

namespace chopper::layout
{

/*!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
    * \sa https://godbolt.org/z/zTj1v9W94
    */
inline std::vector<double>
compute_fp_correction(double const desired_fpr, size_t const num_hash_functions, size_t const requested_max_tb)
{
    std::vector<double> fp_correction{};

    size_t const max_tb = next_multiple_of_64(requested_max_tb);

    fp_correction.resize(max_tb + 1, 0.0);
    fp_correction[1] = 1.0;

    // std::log1p(arg) = std::log(1 + arg). More precise than std::log(1 + arg) if arg is close to zero.
    double const numerator = std::log1p(-std::exp(std::log(desired_fpr) / num_hash_functions));

    for (size_t split = 2u; split <= max_tb; ++split)
    {
        double const log_target_fpr = std::log1p(-std::exp(std::log1p(-desired_fpr) / split));
        fp_correction[split] = numerator / std::log1p(-std::exp(log_target_fpr / num_hash_functions));
        assert(fp_correction[split] >= 1.0);
    }

    return fp_correction;
}

} // namespace chopper::layout
