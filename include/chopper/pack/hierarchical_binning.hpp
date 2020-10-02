#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <seqan3/range/views/to.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/pack/pack_config.hpp>
#include <chopper/pack/print_matrix.hpp>
#include <chopper/pack/simple_binning.hpp>

struct hierarchical_binning
{
    double const alpha{10}; // scale low-level IBF impact

    std::vector<std::string> & names;
    std::vector<size_t> & user_bin_kmer_counts;

    size_t const num_user_bins;
    size_t const num_technical_bins;
    size_t const kmer_count_sum;
    size_t const kmer_count_average_per_bin;

    // ctor
    hierarchical_binning(std::vector<std::string> & names_, std::vector<size_t> & input, pack_config const & config) :
        names{names_},
        user_bin_kmer_counts{input},
        num_user_bins{input.size()},
        num_technical_bins{(config.bins == 0) ? ((user_bin_kmer_counts.size() + 63) / 64 * 64) : config.bins},
        kmer_count_sum{std::accumulate(user_bin_kmer_counts.begin(), user_bin_kmer_counts.end(), 0u)},
        kmer_count_average_per_bin{std::max<size_t>(1u, kmer_count_sum / num_technical_bins)}
    {
        std::cout << "#Techincal bins: " << num_technical_bins << std::endl;
        std::cout << "#User bins: " << input.size() << std::endl;
    }

    template <typename names_type, typename distribution_type>
    void sort_by_distribution(names_type & names, distribution_type & distribution)
    {
        // generate permutation of indices sorted in descinding order by the sequence lengths
        auto permutation = std::views::iota(0u, distribution.size()) | seqan3::views::to<std::vector>;
        assert(permutation.size() == distribution.size());
        auto distribution_compare = [&distribution] (auto const i1, auto const i2)
                                        { return distribution[i2] < distribution[i1]; };
        std::sort(permutation.begin(), permutation.end(), distribution_compare);

        // apply permutation on names and distribution
        for (size_t i = 0; i < permutation.size(); i++)
        {
            auto current = i;
            while (i != permutation[current])
            {
                auto next = permutation[current];
                std::swap(names[current], names[next]);
                std::swap(distribution[current], distribution[next]);
                permutation[current] = current;
                current = next;
            }
            permutation[current] = current;
        }
    }

    void dp_algorithm()
    {
        sort_by_distribution(names, user_bin_kmer_counts);
        // seqan3::debug_stream << std::endl << "Sorted list: " << user_bin_kmer_counts << std::endl << std::endl;

        std::vector<std::vector<size_t>> matrix(num_technical_bins); // rows
        for (auto & v : matrix)
            v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

        std::vector<std::vector<size_t>> ll_matrix(num_technical_bins); // rows
        for (auto & v : ll_matrix)
            v.resize(num_user_bins, 0u); // columns

        std::vector<std::vector<std::pair<size_t, size_t>>> trace(num_technical_bins); // rows
        for (auto & v : trace)
            v.resize(num_user_bins, {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()}); // columns

        // initialize first column (first row is initialized with inf)
        for (size_t i = 0; i < num_technical_bins; ++i)
        {
            matrix[i][0] = user_bin_kmer_counts[0] / (i + 1);
            trace[i][0] = {0u, 0u}; // unnecessary?
        }

        // initialize first row
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            matrix[0][j] = user_bin_kmer_counts[j] + matrix[0][j - 1];
            ll_matrix[0][j] = matrix[0][j];
            trace[0][j] = {0u, j - 1}; // unnecessary?
        }

        // we must iterate column wise
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            size_t const current_weight = user_bin_kmer_counts[j];

            for (size_t i = 1; i < num_technical_bins; ++i)
            {
                size_t minimum{std::numeric_limits<size_t>::max()};
                size_t full_minimum{std::numeric_limits<size_t>::max()};

                for (size_t i_prime = 0; i_prime < i; ++i_prime)
                {
                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                    size_t score = std::max<size_t>(current_weight / (i - i_prime), matrix[i_prime][j-1]);
                    size_t full_score = score * (i + 1) /*#TBs*/ + alpha * ll_matrix[i_prime][j-1];

                    // std::cout << "j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                    if (full_score < full_minimum)
                    {
                        minimum = score;
                        full_minimum = full_score;
                        trace[i][j] = {i_prime, j - 1};
                        ll_matrix[i][j] = ll_matrix[i_prime][j - 1];
                    }
                }

                size_t j_prime{j - 1};
                size_t weight{current_weight};

                // if the user bin j-1 was not split into multiple technical bins!
                // I may merge the current user bin j into the former
                while (j_prime != 0 && ((i - trace[i][j_prime].first) < 2) && weight < minimum)
                {
                    weight += user_bin_kmer_counts[j_prime];
                    --j_prime;

                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                    size_t score = std::max<size_t>(weight, matrix[i - 1][j_prime]);
                    size_t full_score = score * (i + 1) /*#TBs*/ + alpha * (ll_matrix[i - 1][j_prime] + weight);

                    if (full_score < full_minimum)
                    {
                        minimum = score;
                        full_minimum = full_score;
                        trace[i][j] = {i - 1, j_prime};
                        ll_matrix[i][j] = ll_matrix[i - 1][j_prime] + weight;
                    }
                }

                matrix[i][j] = minimum;
            }
        }

        // print_matrix(matrix, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());
        // print_matrix(ll_matrix, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());
        // print_matrix(trace, num_technical_bins, num_user_bins, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));

        // backtracking
        size_t trace_i = num_technical_bins - 1;
        int trace_j = num_user_bins - 1;
        std::cout << "optimum: " << matrix[trace_i][trace_j] << std::endl;
        std::cout << std::endl;

        std::ofstream output_file{"output.binning"};

        output_file << "BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE" << std::endl;

        size_t bin_id{};

        while (trace_j >= 0)
        {
            // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
            size_t next_i = trace[trace_i][trace_j].first;
            size_t next_j = trace[trace_i][trace_j].second;

            size_t kmer_count = user_bin_kmer_counts[trace_j];
            size_t number_of_bins = (trace_i - next_i);

            if (trace_j == 0)
            {
                // we only arrive here if the bin wasn't merged with some before so it is safe to assume
                // that the bin was split (even if only into 1 bin).

                ++trace_i; // because we want the length not the index
                int const kmer_count = user_bin_kmer_counts[0];
                int const average_bin_size = kmer_count / trace_i;

                output_file << "SPLIT_BIN_" << bin_id << '\t'
                            << names[0] << '\t'
                            << trace_i << '\t'
                            << average_bin_size << '\n';

                --trace_j;
                // std::cout << "split " << trace_j << " into " << trace_i << ": " << kmer_count / trace_i << std::endl;
            }
            else if (number_of_bins == 0) // start of merged bin
            {
                std::vector<size_t> merged_bins{kmer_count};
                std::vector<std::string> merged_bin_names{names[trace_j]};
                // std::cout << "merged [" << trace_j;
                while (trace_j > 0 && next_i == trace_i)
                {
                    trace_i = next_i; // unnecessary?
                    --trace_j;
                    kmer_count += user_bin_kmer_counts[trace_j];
                    merged_bins.push_back(user_bin_kmer_counts[trace_j]);
                    merged_bin_names.push_back(names[trace_j]);
                    next_i = trace[trace_i][trace_j].first;
                    // std::cout << "," << trace_j;
                }
                assert(trace_j == 0 || trace_i - next_i == 1);
                assert(kmer_count == std::accumulate(merged_bins.begin(), merged_bins.end(), 0u));

                ++number_of_bins;
                trace_i = next_i;
                --trace_j;

                // now do the binning for the low-level IBF:
                std::string const merged_ibf_name{"MERGED_BIN_" + std::to_string(bin_id)};
                simple_binning algo{merged_bins, merged_bin_names, merged_ibf_name, output_file};
                algo.dp_algorithm();

                // std::cout << "]: " << kmer_count << std::endl;
                // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
            }
            else if (number_of_bins == 1 && next_j != static_cast<size_t>(trace_j) - 1) // merged bin
            {
                std::vector<size_t> merged_bins{kmer_count};
                std::vector<std::string> merged_bin_names{names[trace_j]};
                // std::cout << "merged [" << trace_j;
                --trace_j;
                while (static_cast<size_t>(trace_j) != next_j)
                {
                    kmer_count += user_bin_kmer_counts[trace_j];
                    merged_bins.push_back(user_bin_kmer_counts[trace_j]);
                    merged_bin_names.push_back(names[trace_j]);
                    // std::cout << "," << trace_j;
                    --trace_j;
                }
                trace_i = next_i;
                trace_j = next_j; // unneccessary?

                // now do the binning for the low-level IBF:
                std::string const merged_ibf_name{"COLORFUL_MERGED_BIN_" + std::to_string(bin_id)};
                simple_binning algo{merged_bins, merged_bin_names, merged_ibf_name, output_file};
                algo.dp_algorithm();
                // std::cout << "]: " << kmer_count << std::endl;
            }
            else
            {
                size_t const kmer_count_per_bin = kmer_count / number_of_bins; // round down

                output_file << "SPLIT_BIN_" << bin_id << '\t'
                            << names[trace_j] << '\t'
                            << number_of_bins << '\t'
                            << kmer_count_per_bin << '\n';

                // std::cout << "split " << trace_j << " into " << number_of_bins << ": " << kmer_count_per_bin << std::endl;

                trace_i = trace[trace_i][trace_j].first;
                --trace_j;
            }
            ++bin_id;
        }
    }
};
