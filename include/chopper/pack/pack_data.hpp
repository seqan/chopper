#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <chopper/helper.hpp>
#include <chopper/pack/previous_level.hpp>

struct pack_data
{
    //!\brief The file names of the user input. Since the input might be sorted, we need to keep track of the names.
    std::vector<std::string> filenames;
    //!\brief The kmer counts associated with the above files used to pack user bin into technical bins.
    std::vector<size_t> kmer_counts;
    std::vector<std::vector<std::string>> extra_information;
    std::vector<double> fp_correction{};

    //!\brief A reference to the output stream to cache the results to.
    std::stringstream * output_buffer{nullptr};
    //!\brief A reference to the stream to cache the header to.
    std::stringstream * header_buffer{nullptr};
    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level previous{};

    //!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
    void compute_fp_correction(double const fp_rate, size_t const num_hash_functions)
    {
        size_t const max_tb = next_multiple_of_64(kmer_counts.size());
        fp_correction.resize(max_tb + 1, 0.0);

        double const denominator = std::log(1 - std::exp(std::log(fp_rate) / num_hash_functions));

        for (size_t i = 1; i <= max_tb; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fp_rate, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::exp(std::log(tmp) / num_hash_functions)) / denominator;
            assert(fp_correction[i] >= 1.0);
        }
    }
};
