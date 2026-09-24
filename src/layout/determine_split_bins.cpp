// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <numeric> // for allocator, string
#include <vector>  // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_split_bins.hpp>

#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/print_matrix.hpp> // for data_store
#include <hibf/misc/divide_and_ceil.hpp>

namespace chopper::layout
{

std::pair<size_t, size_t> determine_split_bins(chopper::configuration const & config,
                                               std::vector<size_t> const & positions,
                                               std::vector<size_t> const & cardinalities,
                                               size_t const num_technical_bins,
                                               size_t const num_user_bins,
                                               std::vector<std::vector<size_t>> & partitions)
{
    if (num_user_bins == 0)
        return {0, 0};

    assert(num_technical_bins > 0u);
    assert(num_user_bins > 0u);
    assert(num_user_bins <= num_technical_bins);

    auto const fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = config.hibf_config.maximum_fpr, //
                                                     .hash_count = config.hibf_config.number_of_hash_functions,
                                                     .t_max = num_technical_bins});

    std::vector<std::vector<size_t>> matrix(num_technical_bins); // rows
    for (auto & v : matrix)
        v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

    std::vector<std::vector<size_t>> trace(num_technical_bins); // rows
    for (auto & v : trace)
        v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

    size_t const extra_bins = num_technical_bins - num_user_bins + 1;

    // initialize first column (first row is initialized with inf)
    double const ub_cardinality = static_cast<double>(cardinalities[positions[0]]);
    for (size_t i = 0; i < extra_bins; ++i)
    {
        size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * fpr_correction[i + 1]);
        matrix[i][0] = seqan::hibf::divide_and_ceil(corrected_ub_cardinality, i + 1u);
    }

    // we must iterate column wise
    for (size_t j = 1; j < num_user_bins; ++j)
    {
        double const ub_cardinality = static_cast<double>(cardinalities[positions[j]]);

        for (size_t i = j; i < j + extra_bins; ++i)
        {
            size_t minimum{std::numeric_limits<size_t>::max()};

            for (size_t i_prime = j - 1; i_prime < i; ++i_prime)
            {
                size_t const corrected_ub_cardinality =
                    static_cast<size_t>(ub_cardinality * fpr_correction[(i - i_prime)]);
                size_t score = std::max<size_t>(seqan::hibf::divide_and_ceil(corrected_ub_cardinality, i - i_prime),
                                                matrix[i_prime][j - 1]);

                // std::cout << "j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                minimum = (score < minimum) ? (trace[i][j] = i_prime, score) : minimum;
            }

            matrix[i][j] = minimum;
        }
    }

    // seqan::hibf::layout::print_matrix(matrix, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());
    //seqan::hibf::layout::print_matrix(trace, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());

    // backtracking
    // first, in the last column, find the row with minimum score (it can happen that the last rows are equally good)
    // the less rows, the better
    size_t trace_i{num_technical_bins - 1};
    size_t best_score{std::numeric_limits<size_t>::max()};
    for (size_t best_i{0}; best_i < num_technical_bins; ++best_i)
    {
        if (matrix[best_i][num_user_bins - 1] < best_score)
        {
            best_score = matrix[best_i][num_user_bins - 1];
            trace_i = best_i;
        }
    }
    size_t const number_of_split_tbs{trace_i + 1}; // trace_i is the position. So +1 for the number

    // now that we found the best trace_i start usual backtracking
    size_t trace_j = num_user_bins - 1;

    // size_t max_id{};
    size_t max_size{};

    size_t bin_id{};

    while (trace_j > 0)
    {
        size_t next_i = trace[trace_i][trace_j];
        size_t const number_of_bins = (trace_i - next_i);
        size_t const cardinality = cardinalities[positions[trace_j]];
        size_t const corrected_cardinality = static_cast<size_t>(cardinality * fpr_correction[number_of_bins]);
        size_t const cardinality_per_bin = seqan::hibf::divide_and_ceil(corrected_cardinality, number_of_bins);

        if (cardinality_per_bin > max_size)
        {
            // max_id = bin_id;
            max_size = cardinality_per_bin;
        }

        for (size_t splits{0}; splits < number_of_bins; ++splits)
        {
            partitions[partitions.size() - 1 - bin_id].push_back(positions[trace_j]);
            ++bin_id;
        }

        trace_i = trace[trace_i][trace_j];
        --trace_j;
    }
    ++trace_i; // because we want the length not the index. Now trace_i == number_of_bins
    size_t const cardinality = cardinalities[positions[0]];
    size_t const corrected_cardinality = static_cast<size_t>(cardinality * fpr_correction[trace_i]);
    // NOLINTNEXTLINE(clang-analyzer-core.DivideZero)
    size_t const cardinality_per_bin = seqan::hibf::divide_and_ceil(corrected_cardinality, trace_i);

    if (cardinality_per_bin > max_size)
    {
        // max_id = bin_id;
        max_size = cardinality_per_bin;
    }

    for (size_t splits{0}; splits < trace_i; ++splits)
    {
        partitions[partitions.size() - 1 - bin_id].push_back(positions[0]);
        ++bin_id;
    }

    return {number_of_split_tbs, max_size};
}

} // namespace chopper::layout
