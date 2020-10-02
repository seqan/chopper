#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <chopper/pack/print_matrix.hpp>

struct simple_binning
{
    std::vector<size_t> const & user_bin_kmer_counts;
    std::vector<std::string> const & user_bin_names;
    std::string const & ibf_name;
    std::ostream & output_file;

    size_t const num_user_bins;
    size_t const num_technical_bins;
    size_t const kmer_count_sum;
    size_t const kmer_count_average_per_bin;

    // ctor
    simple_binning(std::vector<size_t> const & input,
                   std::vector<std::string> const & names,
                   std::string const & ibf_name_,
                   std::ostream & file_handle,
                   size_t num_bins = 0) :
        user_bin_kmer_counts{input},
        user_bin_names{names},
        ibf_name{ibf_name_},
        output_file{file_handle},
        num_user_bins{input.size()},
        num_technical_bins{(num_bins == 0) ? ((user_bin_kmer_counts.size() + 63) / 64 * 64) : num_bins},
        kmer_count_sum{std::accumulate(user_bin_kmer_counts.begin(), user_bin_kmer_counts.end(), 0u)},
        kmer_count_average_per_bin{std::max<size_t>(1u, kmer_count_sum / num_technical_bins)}
    {
        std::cout << "#Techincal bins: " << num_technical_bins << std::endl;
        std::cout << "#User bins: " << input.size() << std::endl;
    }

    void dp_algorithm()
    {
        // write result to output_file

        if (num_technical_bins < num_user_bins)
            throw std::logic_error{"this algorithm only works if num_technical_bins >= num_user_bins."};

        std::vector<std::vector<size_t>> matrix(num_technical_bins); // rows
        for (auto & v : matrix)
            v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

        std::vector<std::vector<size_t>> trace(num_technical_bins); // rows
        for (auto & v : trace)
            v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

        size_t const extra_bins = num_technical_bins - num_user_bins + 1;

        // initialize first column (first row is initialized with inf)
        for (size_t i = 0; i < extra_bins; ++i)
        {
            matrix[i][0] = user_bin_kmer_counts[0] / (i + 1);
        }

        // we must iterate column wise
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            int current_weight = user_bin_kmer_counts[j];

            for (size_t i = j; i < j + extra_bins; ++i)
            {
                size_t minimum{std::numeric_limits<size_t>::max()};

                for (size_t i_prime = j - 1; i_prime < i; ++i_prime)
                {
                    size_t score = std::max<int>(current_weight / (i - i_prime), matrix[i_prime][j-1]);

                    // std::cout << "j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                    minimum = (score < minimum) ? (trace[i][j] = i_prime, score) : minimum;
                }

                matrix[i][j] = minimum;
            }
        }

        // print_matrix(matrix, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());
        // print_matrix(trace, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());

        // backtracking
        size_t trace_i = num_technical_bins - 1;
        size_t trace_j = num_user_bins - 1;

        size_t bin_id{};
        while (trace_j > 0)
        {
            size_t next_i = trace[trace_i][trace_j];
            size_t const kmer_count = user_bin_kmer_counts[trace_j];
            size_t const number_of_bins = (trace_i - next_i);
            size_t const kmer_count_per_bin = (kmer_count + number_of_bins - 1) / number_of_bins; // round up

            // ouput_file << IBF_ID,NAME,NUM_TECHNICAL_BINS,ESTIMATED_TB_SIZE
            output_file << ibf_name << '_' << bin_id << '\t'
                        << user_bin_names[trace_j] << '\t'
                        << number_of_bins << '\t'
                        << kmer_count_per_bin << '\n';
            ++bin_id;

            trace_i = trace[trace_i][trace_j];
            --trace_j;
        }
        ++trace_i; // because we want the length not the index. Now trace_i == number_of_bins
        size_t const kmer_count = user_bin_kmer_counts[0];
        size_t const kmer_count_per_bin =  (kmer_count + trace_i - 1) / trace_i;

        // ouput_file << IBF_ID,NAME,NUM_TECHNICAL_BINS,ESTIMATED_TB_SIZE
        output_file << ibf_name << '_' << bin_id << '\t'
                    << user_bin_names[0] << '\t'
                    << trace_i << '\t'
                    << kmer_count_per_bin << '\n';
    }
};
