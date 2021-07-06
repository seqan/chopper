#pragma once

#include <cassert>
#include <cmath>

#include <seqan3/utility/views/to.hpp>

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/pack/pack_config.hpp>
#include <chopper/pack/simple_binning.hpp>
#include <chopper/union/user_bin_sequence.hpp>

class hierarchical_binning
{
private:
    //!\brief The file names of the user input. Since the input might be sorted, we need to keep track of the names.
    std::vector<std::string> & names;
    //!\brief The kmer counts associated with the above files used to pack user bin into technical bins.
    std::vector<size_t> & user_bin_kmer_counts;

    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level const previous;

    //!\brief The user configuration passed down from the command line.
    pack_config const config;

    //!\brief The number of hash functions for the IBFs.
    size_t const num_hash_functions;
    //!\brief The desired false positive rate of the IBFs.
    double const fp_rate;

    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double const alpha;
    //!\brief The number of user bins, initialised with the length of user_bin_kmer_counts.
    size_t const num_user_bins;
    //!\brief The number of technical bins requested by the user.
    size_t const num_technical_bins;
    //!\brief The total sum of all values in user_bin_kmer_counts.
    size_t const kmer_count_sum;
    //!\brief The average count calculated from kmer_count_sum / num_technical_bins.
    size_t const kmer_count_average_per_bin;

    //!\brief If given, the hll sketches are dumped to this directory and restored when they already exist.
    std::filesystem::path const & hll_dir;

    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool const estimate_union;
    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool const rearrange_bins;
    //!\brief The maximal cardinality ratio in the clustering intervals.
    double const max_ratio;
    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t const num_threads;

    //!\brief A reference to the output stream to cache the results to.
    std::stringstream & output_buff;
    //!\brief A reference to the stream to cache the header to.
    std::stringstream & header_buff;

public:
    hierarchical_binning() = delete; //!< Deleted.
    hierarchical_binning(hierarchical_binning const &) = delete; //!< Deleted.
    hierarchical_binning & operator=(hierarchical_binning const &) = delete; //!< Deleted.
    hierarchical_binning(hierarchical_binning &&) = default; //!< Defaulted.
    hierarchical_binning & operator=(hierarchical_binning &&) = default; //!< Defaulted.
    ~hierarchical_binning() = default; //!< Defaulted.

    /*!\brief The constructor from user bin names, their kmer counts and a configuration.
     * \param[in, out] data The data input: filenames associated with the user bin and a kmer count per user bin.
     * \param[in] config_ A configuration object that holds information from the user that influence the computation.
     *
     *
     * Each entry in the names_ and input vector respectively is considered a user bin (both vectors must have the
     * same length).
     */
    hierarchical_binning(pack_data & data, pack_config const & config_) :
        names{data.filenames},
        user_bin_kmer_counts{data.kmer_counts},
        previous{data.previous},
        config{config_},
        num_hash_functions{config.num_hash_functions},
        fp_rate{config.fp_rate},
        alpha{config.alpha},
        num_user_bins{data.kmer_counts.size()},
        num_technical_bins{(config.bins == 0) ? ((user_bin_kmer_counts.size() + 63) / 64 * 64) : config.bins},
        kmer_count_sum{std::accumulate(user_bin_kmer_counts.begin(), user_bin_kmer_counts.end(), 0u)},
        kmer_count_average_per_bin{std::max<size_t>(1u, kmer_count_sum / num_technical_bins)},
        hll_dir{config.hll_dir},
        estimate_union{config.union_estimate},
        rearrange_bins{config.rearrange_bins},
        max_ratio{config.max_ratio},
        num_threads{config.num_threads},
        output_buff{*data.output_buffer},
        header_buff{*data.header_buffer}
    {
        // std::cout << "#Techincal bins: " << num_technical_bins << std::endl;
        // std::cout << "#User bins: " << data.kmer_counts.size() << std::endl;
        // std::cout << "Output file: " << config.output_filename.string() << std::endl;

        if (names.size() != user_bin_kmer_counts.size())
            throw std::runtime_error{"The filenames and kmer counts do not have the same length."};
    }

    //!\brief Executes the hierarchical binning algorithm and packs user bins into technical bins.
    size_t execute()
    {
        static constexpr size_t max_size_t{std::numeric_limits<size_t>::max()};

        sort_by_distribution(names, user_bin_kmer_counts);
        // seqan3::debug_stream << std::endl << "Sorted list: " << user_bin_kmer_counts << std::endl << std::endl;

        // precompute f_h factors that adjust the split bin size 
        // to prevent FPR inflation due to multiple testing 
        std::vector<double> fp_correction(num_technical_bins + 1, 0.0);

        double const h_inverse = 1.0 / static_cast<double>(num_hash_functions);
        double const denominator = std::log(1 - std::pow(fp_rate, h_inverse));

        for (size_t i = 1; i <= num_technical_bins; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fp_rate, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::pow(tmp, h_inverse)) / denominator;
        }

        std::vector<std::vector<uint64_t>> union_estimates;
        // Depending on cli flags given, use HyperLogLog estimates and/or rearrangement algorithms
        if (estimate_union)
        {
            user_bin_sequence bin_sequence(names, user_bin_kmer_counts, hll_dir);

            if (rearrange_bins)
                bin_sequence.rearrange_bins(max_ratio, num_threads);

            bin_sequence.estimate_interval_unions(union_estimates, num_threads);
        }
    	
        // technical bins (outer) = rows; user bins (inner) = columns
        std::vector<std::vector<size_t>> matrix(num_technical_bins,
                                                std::vector<size_t>(num_user_bins, max_size_t));

        // technical bins (outer) = rows; user bins (inner) = columns
        std::vector<std::vector<size_t>> ll_matrix(num_technical_bins,
                                                   std::vector<size_t>(num_user_bins, 0u));

        // technical bins (outer) = rows; user bins (inner) = columns
        std::vector<std::vector<std::pair<size_t, size_t>>> trace(num_technical_bins,
                                                                  std::vector<std::pair<size_t, size_t>>(
                                                                      num_user_bins, {max_size_t, max_size_t}));

        initialization(matrix, ll_matrix, trace, fp_correction, union_estimates);

        recursion(matrix, ll_matrix, trace, fp_correction, union_estimates);

        // print_matrix(matrix, num_technical_bins, num_user_bins, max_size_t);
        // print_matrix(ll_matrix, num_technical_bins, num_user_bins, max_size_t);
        // print_matrix(trace, num_technical_bins, num_user_bins, std::make_pair(max_size_t, max_size_t));

        return backtracking(matrix, ll_matrix, trace);
    }

private:
    /*!\brief Sorts both input vectors (names and distribution) only by looking at the values in `distribution`.
     * \param[in, out] names The names to be sorted in parallel to the `distribution` vector.
     * \param[in, out] distribution The vector to be used to sort both input vectors by.
     */
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

    /*!\brief Initialize the matrices M (hig_level_ibf), L (low_level_ibfs) and T (trace)
     *
     * \image html hierarchical_dp_init.png
     */
    void initialization(std::vector<std::vector<size_t>> & matrix,
                        std::vector<std::vector<size_t>> & ll_matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> & trace,
                        std::vector<double> const & fp_correction,
                        std::vector<std::vector<uint64_t>> const & union_estimates)
    {
        // initialize first column
        double const ub_cardinality = static_cast<double>(user_bin_kmer_counts[0]);
        for (size_t i = 0; i < num_technical_bins; ++i)
        {
            size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * fp_correction[i + 1]);
            matrix[i][0] = corrected_ub_cardinality / (i + 1);
            trace[i][0] = {0u, 0u}; // unnecessary?
        }

        // initialize first row
        if (estimate_union)
        {
            for (size_t j = 1; j < num_user_bins; ++j)
            {
                ll_matrix[0][j] = user_bin_kmer_counts[j] + matrix[0][j - 1];
                matrix[0][j] = union_estimates[j][0];
                trace[0][j] = {0u, j - 1}; // unnecessary?
            }
        }
        else
        {
            for (size_t j = 1; j < num_user_bins; ++j)
            {
                size_t const sum = user_bin_kmer_counts[j] + matrix[0][j - 1];
                matrix[0][j] = sum;
                ll_matrix[0][j] = sum;
                trace[0][j] = {0u, j - 1}; // unnecessary?
            }
        }
    }

    /*!\brief Performs the recursion.
     *
     * \image html hierarchical_dp_recursion.png
     *
     * Explanations to the formula:
     *
     * Remember that M (matrix) stores the maximum technical bin size of the high level IBF (HIBF) and
     * L (ll_matrix) stores the sum of all estimate low level ibfs (LIBF) memory footprints
     * (assuming perfect splitting of kmer_content which can obviously not be achieved but `alpha` may be adjusted
     * in order to counteract this underestimation).
     *
     * Now in order to minimize the memory footprint we...
     *
     * 1. ... (\f$v_ij\f$) firstly check if we should split the bin. <br>
     * Therefore we compute for every possible \f$ i' \f$ above \f$ i \f$, the technical bin size if we split
     * \f$ c_j \f$ into \f$ i - i' \f$ bins (\f$ \frac{c_j}{i - i'} \f$). We only take the maximum of the current
     * maximum where I come from (\f$ M_{i',j-1} \f$) and the new technical bin size computed just now.
     * This maximum is the new current *maximum technical bin size* that needs to be multiplied by the current number
     * technical bins \f$ (i + 1) \f$ in order to estimate the memory footprint of the HIBF.
     * Now that we have the memory footprint of the HIBF we also consider the LIBFs memory footprint but nothing
     * changed here, since we split not merge, so we just take \f$ L_{i',j-1} \f$ scaled by alpha.
     *
     * 2. ... (\f$h_ij\f$) secondly check if we should merge the bin with the ones before.<br>
     * Therefore we start by compute for every possible \f$ j' \f$ to the left of \f$ j \f$, the merged bin weight
     * (\f$ \sum_{g = j'}^{j} c_g \f$) of merging together user bins \f$ [j', ...,  j] \f$. This is only possible
     * iff every bin \f$ [j', ...,  j] \f$ was **not splitted**. If we merge those bins starting with user bin
     * \f$ j' \f$ we start the trace at \f$ M_{i-1,j'-1} \f$. Therefore we need to compute the new maximal technical
     * bin size of the HIBF, which is the maximum of the merged bin weight and where we would come from
     * (\f$ \max(M_{i-1,j'-1}, \sum_{g = j'}^{j} c_g) \f$) and multiple it by the number of technical bins so far
     * \f$ (i + 1) \f$ to get the HIBF memory footprint. The LIBFs memory footprint also changes since we introduce a
     * new merged bin. Namely, we add the weight of the new merged bin (\f$ \sum_{g = j'}^{j} c_g \f$) again to the
     * LIBFs memory footprint from where we would come from \f$ L_{i-1,j'-1} \f$ and scale this by alpha. Just adding
     * the combined merged bin weight neglects the fact, that the merged bin weight has to be distributed within the
     * new low level IBF, potentially causing the effective text ratio and thereby the IBF memory footprint to increase.
     * This we cannot know before hand how the data are, we need to accept this as a flaw in the "optimumal result" of
     * this algorithm. It would be too computational intensive to compute the splitting for every possibility.
     *
     */
    void recursion(std::vector<std::vector<size_t>> & matrix,
                   std::vector<std::vector<size_t>> & ll_matrix,
                   std::vector<std::vector<std::pair<size_t, size_t>>> & trace,
                   std::vector<double> const & fp_correction,
                   std::vector<std::vector<uint64_t>> const & union_estimates)
    {
        // we must iterate column wise
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            size_t const current_weight = user_bin_kmer_counts[j];
            double const ub_cardinality = static_cast<double>(current_weight);

            for (size_t i = 1; i < num_technical_bins; ++i)
            {
                size_t minimum{std::numeric_limits<size_t>::max()};
                size_t full_minimum{std::numeric_limits<size_t>::max()};

                // check vertical cells
                for (size_t i_prime = 0; i_prime < i; ++i_prime)
                {
                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                    size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * fp_correction[(i - i_prime)]);
                    size_t score = std::max<size_t>(corrected_ub_cardinality / (i - i_prime), matrix[i_prime][j-1]);
                    size_t full_score = score * (i + 1) /*#TBs*/ + alpha * ll_matrix[i_prime][j-1];

                    // std::cout << " ++ j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                    if (full_score < full_minimum)
                    {
                        minimum = score;
                        full_minimum = full_score;
                        trace[i][j] = {i_prime, j - 1};
                        ll_matrix[i][j] = ll_matrix[i_prime][j - 1];
                    }
                }

                // seqan3::debug_stream << "current vertical minimum of " << "j:" << j << " i:" << i
                //                      << " -> score:" << full_minimum << " (M_ij=" << minimum << ")"
                //                      << " trace:" << trace[i][j]
                //                      << std::endl;

                // check horizontal cells
                size_t j_prime{j - 1};
                size_t weight{current_weight};

                auto get_weight = [&] ()
                {
                    // if we use the union estimate we plug in that value instead of the sum (weight)
                    // union_estimate[i][j] is the union of {j, ..., i}
                    // the + 1 is necessary because j_prime is decremented directly after weight is updated
                    return estimate_union ? union_estimates[j][j_prime + 1] : weight;
                };

                // if the user bin j-1 was not split into multiple technical bins!
                // I may merge the current user bin j into the former
                while (j_prime != 0 && ((i - trace[i][j_prime].first) < 2) && get_weight() < minimum)
                {
                    weight += user_bin_kmer_counts[j_prime];
                    --j_prime;

                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                    size_t const score = std::max<size_t>(get_weight(), matrix[i - 1][j_prime]);
                    size_t const full_score = score * (i + 1) /*#TBs*/ + alpha * (ll_matrix[i - 1][j_prime] + weight);

                    // seqan3::debug_stream << " -- " << "j_prime:" << j_prime
                    //                      << " -> full_score:" << full_score << " (M_{i-1,j'}=" << score << ")"
                    //                      << std::endl;

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
    }

    //!\brief Backtracks the trace matrix and writes the resulting binning into the output file.
    size_t backtracking(std::vector<std::vector<size_t>> const & matrix,
                        std::vector<std::vector<size_t>> const & ll_matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> const & trace)
    {
        bool const high = previous.empty();

        // backtracking
        size_t trace_i = num_technical_bins - 1;
        int trace_j = num_user_bins - 1;
        // std::cout << "optimum: " << matrix[trace_i][trace_j] << std::endl;
        // std::cout << std::endl;

        if (output_buff.tellp() == 0) // beginning of the file
            output_buff << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES" << std::endl;

        size_t high_level_max_id{};
        size_t high_level_max_size{};

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

                output_buff << names[0] << '\t'
                            << previous.bin_indices  << (high ? "" : ";") << bin_id << '\t'
                            << previous.num_of_bins  << (high ? "" : ";") << trace_i << '\t'
                            << previous.estimated_sizes << (high ? "" : ";") << average_bin_size << '\n';

                if (average_bin_size > high_level_max_size)
                {
                    high_level_max_id = bin_id;
                    high_level_max_size = average_bin_size;
                }

                --trace_j;
                // std::cout << "split " << trace_j << " into " << trace_i << ": " << kmer_count / trace_i << std::endl;
            }
            else if (number_of_bins == 0) // start of merged bin
            {
                pack_data libf_data{};
                libf_data.output_buffer = &output_buff;
                libf_data.header_buffer = &header_buff;

                libf_data.kmer_counts = {kmer_count};
                libf_data.filenames = {names[trace_j]};
                // std::cout << "merged [" << trace_j;
                while (trace_j > 0 && next_i == trace_i)
                {
                    trace_i = next_i; // unnecessary?
                    --trace_j;
                    kmer_count += user_bin_kmer_counts[trace_j];
                    libf_data.kmer_counts.push_back(user_bin_kmer_counts[trace_j]);
                    libf_data.filenames.push_back(names[trace_j]);
                    next_i = trace[trace_i][trace_j].first;
                    // std::cout << "," << trace_j;
                }
                assert(trace_j == 0 || trace_i - next_i == 1);
                assert(kmer_count == std::accumulate(libf_data.kmer_counts.begin(), libf_data.kmer_counts.end(), 0u));

                ++number_of_bins;
                trace_i = next_i;
                --trace_j;

                libf_data.previous = previous;
                libf_data.previous.bin_indices += (high ? "" : ";") + std::to_string(bin_id);
                libf_data.previous.num_of_bins  += (high ? "" : ";") + std::string{"1"};
                libf_data.previous.estimated_sizes += (high ? "" : ";") + std::to_string(kmer_count);

                std::string const merged_ibf_name{std::string{merged_bin_prefix} + "_" + libf_data.previous.bin_indices};

                // now do the binning for the low-level IBF:
                size_t merged_max_bin_id;
                if (libf_data.kmer_counts.size() > num_technical_bins)
                {
                    hierarchical_binning algo{libf_data, config};
                    merged_max_bin_id = algo.execute();
                }
                else
                {
                    simple_binning algo{libf_data};
                    merged_max_bin_id = algo.execute();
                }
                header_buff << "#" << merged_ibf_name << " max_bin_id:" << merged_max_bin_id << '\n';

                if (kmer_count > high_level_max_size)
                {
                    high_level_max_id = bin_id;
                    high_level_max_size = kmer_count;
                }

                // std::cout << "]: " << kmer_count << std::endl;
                // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
            }
            else if (number_of_bins == 1 && next_j != static_cast<size_t>(trace_j) - 1) // merged bin
            {
                pack_data libf_data{};
                libf_data.output_buffer = &output_buff;
                libf_data.header_buffer = &header_buff;

                libf_data.kmer_counts = {kmer_count};
                libf_data.filenames = {names[trace_j]};
                // std::cout << "merged [" << trace_j;
                --trace_j;
                while (static_cast<size_t>(trace_j) != next_j)
                {
                    kmer_count += user_bin_kmer_counts[trace_j];
                    libf_data.kmer_counts.push_back(user_bin_kmer_counts[trace_j]);
                    libf_data.filenames.push_back(names[trace_j]);
                    // std::cout << "," << trace_j;
                    --trace_j;
                }
                trace_i = next_i;
                trace_j = next_j; // unneccessary?

                libf_data.previous = previous;
                libf_data.previous.bin_indices += (high ? "" : ";") + std::to_string(bin_id);
                libf_data.previous.num_of_bins  += (high ? "" : ";") + std::string{"1"};
                libf_data.previous.estimated_sizes += (high ? "" : ";") + std::to_string(kmer_count);

                std::string const merged_ibf_name{std::string{merged_bin_prefix} + "_" + libf_data.previous.bin_indices};

                size_t merged_max_bin_id;
                // now do the binning for the low-level IBF:
                if (libf_data.kmer_counts.size() > num_technical_bins)
                {
                    hierarchical_binning algo{libf_data, config};
                    merged_max_bin_id = algo.execute();
                }
                else
                {
                    simple_binning algo{libf_data};
                    merged_max_bin_id = algo.execute();
                }
                header_buff << "#" << merged_ibf_name << " max_bin_id:" << merged_max_bin_id << '\n';

                if (kmer_count > high_level_max_size)
                {
                    high_level_max_id = bin_id;
                    high_level_max_size = kmer_count;
                }
                // std::cout << "]: " << kmer_count << std::endl;
            }
            else
            {
                size_t const kmer_count_per_bin = kmer_count / number_of_bins; // round down

                output_buff << names[trace_j] << '\t'
                            << previous.bin_indices  << (high ? "" : ";") << bin_id << '\t'
                            << previous.num_of_bins  << (high ? "" : ";") << number_of_bins << '\t'
                            << previous.estimated_sizes << (high ? "" : ";") << kmer_count_per_bin << '\n';

                // std::cout << "split " << trace_j << " into " << number_of_bins << ": " << kmer_count_per_bin << std::endl;

                if (kmer_count_per_bin > high_level_max_size)
                {
                    high_level_max_id = bin_id;
                    high_level_max_size = kmer_count_per_bin;
                }

                trace_i = trace[trace_i][trace_j].first;
                --trace_j;
            }

            bin_id += number_of_bins;
        }

        return high_level_max_id;
    }
};
