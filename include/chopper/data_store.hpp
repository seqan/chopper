#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/helper.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/layout.hpp>
#include <chopper/layout/previous_level.hpp>
#include <chopper/sketch/toolbox.hpp>

namespace chopper
{

struct data_store
{
    /*!\name References to global instances of the HIBF.
     * \{
     */
    //!\brief The desired maximum false positive rate of the resulting index.
    double const false_positive_rate{};

    //!\brief An object starting at the top level IBF to collecting statistics about the HIBF on the way.
    layout::hibf_statistics::level * stats{nullptr};

    //!\brief The layout that is build by layout::hierarchical_binning.
    layout::layout * hibf_layout;
    //!\}

    /*!\name Local Storage one IBF in the HIBF.
     *
     * These member variables change on each IBF of the HIBF s.t. the current IBF can be constructed from
     * the current subset of data. The same data is also used for the top level IBF that holds all the data.
     * \{
     */
    //!\brief The file names of the user input. Since the input might be sorted, we need to keep track of the names.
    std::vector<std::string> filenames{};

    //!\brief The kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> kmer_counts{};

    //!\brief The hyperloglog sketches of all input files to estimate their size and similarities.
    std::vector<sketch::hyperloglog> sketches{};

    //!\brief Extra information in given in the input file
    std::vector<std::string> extra_information_strings{};

    //!\brief Extra information in given in the input file
    std::vector<std::vector<std::string>> extra_information{};

    //!\brief The false positive correction based on fp_rate, num_hash_functions and requested_max_tb.
    std::vector<double> fp_correction{};

    //!\brief Stores sketches if needed and provides utility functions for user bin rearrangement or union estimation.
    sketch::toolbox sketch_toolbox{};

    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    layout::previous_level previous{};

    //!\brief Matrix of estimates of merged bin cardinalites
    std::vector<uint64_t> union_estimates{};

    bool user_bins_arranged{false};
    //!\}

    /*!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
     * \sa https://godbolt.org/z/zTj1v9W94
     */
    void compute_fp_correction(double const desired_fpr, size_t const num_hash_functions, size_t const requested_max_tb)
    {
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
    }
};

} // namespace chopper
