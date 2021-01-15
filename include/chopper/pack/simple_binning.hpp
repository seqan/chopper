#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <chopper/pack/pack_data.hpp>
#include <chopper/pack/print_matrix.hpp>

/*!\brief Distributes x Technical Bins across y User Bins while minimizing the maximal Technical Bin size
 *
 * # Terminology
 *
 * ## Technical Bin
 * \copydetails simple_binning::num_technical_bins
 *
 * ## User Bin
 * \copydetails simple_binning::num_user_bins
 *
 * ## Notation
 *  | Name    | Description                                                                 |
 *  |---------|---------------------------------------------------------------------------- |
 *  | **x**   | Number of Technical Bins (TB)                                               |
 *  | **y**   | Number of User Bins (UB)                                                    |
 *  | **b_i** | The bin size (kmer content) of Technical Bin $i$                            |
 *  | **c_j** | The kmer content of User Bin $j$                                            |
 *  | **M**   | A DP matrix that tracks the maximum technical bin size \f$\max_{i} (b_i)\f$.|
 *
 *  \attention The number of technical bins **x** must be greater that the number of user bins **y** for this algorithm.
     *            If you want to use less technical bins than user bins, see the hierarchical_binning algorithm.
 *
 * # Algorithm
 *
 * Since the size of the IBF depends on the maximal Technical Bin size, we want to minimize \f$ max_{i} (b_i)\f$.
 *
 * Let \f$r = x - y\f$ be the surplus of TBs
 *
 * ## Initialization
 * <img src="execute_init.png" align="left"  width="10%"/>
 *
 * <br><br>
 * \f$\qquad \forall_{i \in [0,r]} \quad M_{i,0} = \frac{c_0}{i + 1}\f$
 * <br><br><br><br><br><br><br><br><br><br>
 *
 * ## Recursion
 * <img src="execute_recursion.png" align="left"  width="10%"/>
 *
 * <br><br>
 * \f$\forall_{i,j} \quad M_{i,j} = \min_{i' \in [i - r - 1, i - 1]} \max(M_{i',j-1}, \frac{c_j}{i - i'})\f$
 * <br><br><br><br><br><br><br><br><br><br>
 *
 * ## Backtracking
 *
 * Assume we filled a trace matrix **T** during the computation of **M**.
 *
 * We now want to recover the number of bins **n_j** for each User Bin **j**.
 *
 * <img src="execute_backtracking.png" align="left"  width="10%"/>
 *
 * Backtracking pseudo code:
 * ```
 * // Start at the bottom-right cell.
 * j = y - 1;
 * i = x - 1;
 * n = array(y); // array of length y
 *
 * while (j > 0)
 * {
 *     next_i = T[i][j];
 *     n[j] = i - next_i;
 *     i = next_i;
 *     j = j - 1;
 * }
 * ```
 */
struct simple_binning
{
private:
    //!\brief The filenames of the respective user bins.
    std::vector<std::string> const & user_bin_names;
    //!\brief The kmer counts of the respective user bins.
    std::vector<size_t> const & user_bin_kmer_counts;
    //!\brief The identifier of the (colorful) bin that is written into the file.
    std::string const & bin_name;
    //!\brief The output stream to write to.
    std::ostream & output_file;

    /*!\brief The number of User bins.
     *
     * The user may impose a structure on his sequence data in the form of *logical groups* (e.g. species).
     * When querying the IBF, the user is interested in an answer that differentiates between these groups.
     */
    size_t const num_user_bins;
    /*!\brief The number of Technical bins.
     *
     * A *Technical Bin* represents an actual bin in the binning directory.
     * In the IBF, it stores its kmers in a **single Bloom Filter** (which is interleaved with all the other BFs).
     */
    size_t const num_technical_bins;
    //!\brief The total sum of all values in user_bin_kmer_counts.
    size_t const kmer_count_sum;
    //!\brief The average count calculated from kmer_count_sum / num_technical_bins.
    size_t const kmer_count_average_per_bin;

public:
    /*!\brief The constructor from user bin names, their kmer counts and a configuration.
     * \param[in] data The filenames and kmer counts associated with the user bin, as well as the ostream buffer.
     * \param[in] bin_name_ The bin identifier to write into the output file.
     * \param[in] num_bins (optional) The number of technical bins.
     *
     * If the `num_bins` parameter is omitted or set to 0, then number of technical bins used in this algorithm
     * is automatically set to the next multiple of 64 given the number of user bins (e.g. #UB = 88 -> #TB = 124).
     *
     * \attention The number of technical bins must be greater or equal to the number of user bins!
     *            If you want to use less technical bins than user bins, see the hierarchical_binning algorithm.
     */
    simple_binning(pack_data & data, std::string const & bin_name_, size_t num_bins = 0) :
        user_bin_names{data.filenames},
        user_bin_kmer_counts{data.kmer_counts},
        bin_name{bin_name_},
        output_file{*data.output_buffer},
        num_user_bins{data.kmer_counts.size()},
        num_technical_bins{(num_bins == 0) ? ((user_bin_kmer_counts.size() + 63) / 64 * 64) : num_bins},
        kmer_count_sum{std::accumulate(user_bin_kmer_counts.begin(), user_bin_kmer_counts.end(), 0u)},
        kmer_count_average_per_bin{std::max<size_t>(1u, kmer_count_sum / num_technical_bins)}
    {
        std::cout << "#Techincal bins: " << num_technical_bins << std::endl;
        std::cout << "#User bins: " << data.kmer_counts.size() << std::endl;

        if (num_user_bins > num_technical_bins)
        {
            throw std::runtime_error{"You cannot have less technical bins than user bins for this simple binning "
                                     "algorithm. Please see the hierarchical_binning algorithm or increase the number "
                                     "of technical bins."};
        }
    }

    //!\brief Executes the simple binning algorithm and packs user bins into technical bins.
    size_t execute()
    {
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

        size_t max_id{};
        size_t max_size{};

        size_t bin_id{};
        while (trace_j > 0)
        {
            size_t next_i = trace[trace_i][trace_j];
            size_t const kmer_count = user_bin_kmer_counts[trace_j];
            size_t const number_of_bins = (trace_i - next_i);
            size_t const kmer_count_per_bin = (kmer_count + number_of_bins - 1) / number_of_bins; // round up

            // columns: IBF_ID,NAME,NUM_TECHNICAL_BINS,ESTIMATED_TB_SIZE
            output_file << bin_name << '_' << bin_id << '\t'
                        << user_bin_names[trace_j] << '\t'
                        << number_of_bins << '\t'
                        << kmer_count_per_bin << '\n';

            if (kmer_count_per_bin > max_size)
            {
                max_id = bin_id;
                max_size = kmer_count_per_bin;
            }

            bin_id += number_of_bins;

            trace_i = trace[trace_i][trace_j];
            --trace_j;
        }
        ++trace_i; // because we want the length not the index. Now trace_i == number_of_bins
        size_t const kmer_count = user_bin_kmer_counts[0];
        size_t const kmer_count_per_bin =  (kmer_count + trace_i - 1) / trace_i;

        if (kmer_count_per_bin > max_size)
        {
            max_id = bin_id;
            max_size = kmer_count_per_bin;
        }

        // columns: IBF_ID,NAME,NUM_TECHNICAL_BINS,ESTIMATED_TB_SIZE
        output_file << bin_name << '_' << bin_id << '\t'
                    << user_bin_names[0] << '\t'
                    << trace_i << '\t'
                    << kmer_count_per_bin << '\n';

        return max_id;
    }
};
