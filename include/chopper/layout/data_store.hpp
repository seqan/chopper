#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <chopper/helper.hpp>
#include <chopper/configuration.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/previous_level.hpp>
#include <chopper/sketch/user_bin_sequence.hpp>

namespace chopper::layout
{

struct data_store
{
    //!\brief The file names of the user input. Since the input might be sorted, we need to keep track of the names.
    std::vector<std::string> filenames{};
    //!\brief The kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> kmer_counts{};
    std::vector<std::vector<std::string>> extra_information{};
    std::vector<double> fp_correction{};

    //!\brief Stores sketches if needed and provides utility functions for user bin rearrangement or union estimation.
    sketch::user_bin_sequence sketch_toolbox;

    //!\brief The desired maximum false positive rate of the resulting index.
    double false_positive_rate{};

    //!\brief Matrix of estimates of merged bin cardinalites
    std::vector<uint64_t> union_estimates{};
    bool user_bins_arranged{false};

    //!\brief A reference to the output stream to cache the results to.
    std::stringstream * output_buffer{nullptr};
    //!\brief A reference to the stream to cache the header to.
    std::stringstream * header_buffer{nullptr};
    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level previous{};

    //!\brief An object starting at the top level IBF to collecting statistics about the HIBF on the way.
    hibf_statistics::level * stats{nullptr};

    //!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
    void compute_fp_correction(double const fp_rate, size_t const num_hash_functions, size_t const requested_max_tb)
    {
        size_t const max_tb = next_multiple_of_64(requested_max_tb);

        fp_correction.resize(max_tb + 1, 0.0);
        fp_correction[1] = 1.0;

        double const denominator = std::log(1 - std::exp(std::log(fp_rate) / num_hash_functions));

        for (size_t i = 2; i <= max_tb; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fp_rate, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::exp(std::log(tmp) / num_hash_functions)) / denominator;
            assert(fp_correction[i] >= 1.0);
        }
    }
};

} // namespace chopper::layout
