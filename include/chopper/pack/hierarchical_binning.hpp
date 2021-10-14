#pragma once

#include <cassert>
#include <cmath>

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/helper.hpp>
#include <chopper/pack/ibf_query_cost.hpp>
#include <chopper/pack/pack_config.hpp>
#include <chopper/pack/simple_binning.hpp>

class hierarchical_binning
{
private:
    //!\brief The user configuration passed down from the command line.
    pack_config const config{};
    //!\brief The data input: filenames associated with the user bin and a kmer count per user bin.
    pack_data * const data{nullptr};

    //!\brief The number of user bins, initialised with the length of user_bin_kmer_counts.
    size_t const num_user_bins{};
    //!\brief The number of technical bins requested by the user.
    size_t const num_technical_bins{};
    //!\brief The cost for querying `num_technical_bins` bins.
    double const interpolated_cost{};

    //!\brief The total query cost of all k-mers.
    double total_query_cost{};

public:
    hierarchical_binning() = default; //!< Defaulted.
    hierarchical_binning(hierarchical_binning const &) = delete; //!< Deleted. Would modify same data.
    hierarchical_binning & operator=(hierarchical_binning const &) = delete; //!< Deleted. Would modify same data.
    hierarchical_binning(hierarchical_binning &&) = default; //!< Defaulted.
    hierarchical_binning & operator=(hierarchical_binning &&) = default; //!< Defaulted.
    ~hierarchical_binning() = default; //!< Defaulted.

    /*!\brief The constructor from user bin names, their kmer counts and a configuration.
     * \param[in, out] data_ The data input: filenames associated with the user bin and a kmer count per user bin.
     * \param[in] config_ A configuration object that holds information from the user that influence the computation.
     *
     * Each entry in the names_ and input vector respectively is considered a user bin (both vectors must have the
     * same length).
     */
    hierarchical_binning(pack_data & data_, pack_config const & config_) :
        config{config_},
        data{std::addressof(data_)},
        num_user_bins{data->kmer_counts.size()},
        num_technical_bins{data->previous.empty() ? config.t_max : needed_technical_bins(num_user_bins)},
        interpolated_cost{ibf_query_cost::interpolated(num_technical_bins)}
    {
        assert(data != nullptr);
        assert(data->output_buffer != nullptr);
        assert(data->header_buffer != nullptr);

        if (config.debug)
        {
            *data->header_buffer << std::fixed << std::setprecision(2);
            *data->output_buffer << std::fixed << std::setprecision(2);
        }

        if (data->filenames.size() != data->kmer_counts.size())
            throw std::runtime_error{"The filenames and kmer counts do not have the same length."};
    }

    //!\brief Executes the hierarchical binning algorithm and packs user bins into technical bins.
    std::tuple<size_t, double> execute()
    {
        assert(data != nullptr);
        assert(data->output_buffer != nullptr);
        assert(data->header_buffer != nullptr);

        static constexpr size_t max_size_t{std::numeric_limits<size_t>::max()};

        data->arrange_user_bins(config);

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

        initialization(matrix, ll_matrix, trace);

        recursion(matrix, ll_matrix, trace);

        // print_matrix(matrix, num_technical_bins, num_user_bins, max_size_t);
        // print_matrix(ll_matrix, num_technical_bins, num_user_bins, max_size_t);
        // print_matrix(trace, num_technical_bins, num_user_bins, std::make_pair(max_size_t, max_size_t));

        return std::make_tuple(backtracking(matrix, trace), total_query_cost);
    }

private:
    /*!\brief Returns the number of technical bins given a number of user bins.
     * \param[in] requested_num_ub The number of user bins.
     */
    [[nodiscard]] size_t needed_technical_bins(size_t const requested_num_ub) const
    {
        return std::min<size_t>(next_multiple_of_64(requested_num_ub), config.t_max);
    }

    /*!\brief Returns the maximum number of needed levels when merging `num_ubs_in_merge` many user bins.
     * \param[in] num_ubs_in_merge The number of user bins in the merge.
     */
    [[nodiscard]] size_t max_merge_levels(size_t const num_ubs_in_merge) const
    {
        size_t const lower_lvl_tbs = needed_technical_bins(num_ubs_in_merge);
        double const levels = std::log(num_ubs_in_merge) / std::log(lower_lvl_tbs);
        return static_cast<size_t>(std::ceil(levels));
    }

    /*!\brief Initialize the matrices M (hig_level_ibf), L (low_level_ibfs) and T (trace)
     *
     * \image html hierarchical_dp_init.png
     */
    void initialization(std::vector<std::vector<size_t>> & matrix,
                        std::vector<std::vector<size_t>> & ll_matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
    {
        assert(data != nullptr);

        // initialize first column
        for (size_t i = 0; i < num_technical_bins; ++i)
        {
            size_t cardinality = data->kmer_counts[0];
            // if this is the top recursive level, we must apply the FPR correction
            if (data->previous.empty())
                cardinality = std::ceil(static_cast<double>(cardinality) * data->fp_correction[i + 1]);

            matrix[i][0] = cardinality / (i + 1);
            trace[i][0] = {0u, 0u}; // unnecessary?
        }

        // initialize first row
        size_t sum = data->kmer_counts[0];
        if (config.estimate_union)
        {
            for (size_t j = 1; j < num_user_bins; ++j)
            {
                sum += data->kmer_counts[j];
                matrix[0][j] = data->union_estimates[j][0];
                ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
                trace[0][j] = {0u, j - 1}; // unnecessary?
            }
        }
        else
        {
            for (size_t j = 1; j < num_user_bins; ++j)
            {
                sum += data->kmer_counts[j];
                matrix[0][j] = sum;
                ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
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
                   std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
    {
        assert(data != nullptr);

        // we must iterate column wise
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            size_t const current_weight = data->kmer_counts[j];

            for (size_t i = 1; i < num_technical_bins; ++i)
            {
                size_t minimum{std::numeric_limits<size_t>::max()};
                size_t full_minimum{std::numeric_limits<size_t>::max()};

                // check vertical cells
                for (size_t i_prime = 0; i_prime < i; ++i_prime)
                {
                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint

                    size_t cardinality = current_weight;
                    // if this is the top recursive level, we must apply the FPR correction
                    if (data->previous.empty())
                        cardinality = std::ceil(static_cast<double>(cardinality) * data->fp_correction[i - i_prime]);
                    
                    size_t score = std::max<size_t>(cardinality / (i - i_prime), matrix[i_prime][j-1]);
                    size_t full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_matrix[i_prime][j-1];

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
                    // union_estimates[i][j] is the union of {j, ..., i}
                    // the + 1 is necessary because j_prime is decremented directly after weight is updated
                    return config.estimate_union ? data->union_estimates[j][j_prime + 1] : weight;
                };

                // if the user bin j-1 was not split into multiple technical bins!
                // I may merge the current user bin j into the former
                while (j_prime != 0 && ((i - trace[i][j_prime].first) < 2) && get_weight() < minimum)
                {
                    weight += data->kmer_counts[j_prime];
                    --j_prime;

                    // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                    // ll_kmers: estimate for the number of k-mers that have to be resolved on lower levels
                    // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                    size_t const score = std::max<size_t>(matrix[i - 1][j_prime], get_weight());
                    size_t const ll_kmers = ll_matrix[i - 1][j_prime] + max_merge_levels(j - j_prime) * weight;
                    size_t const full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_kmers;

                    // seqan3::debug_stream << " -- " << "j_prime:" << j_prime
                    //                      << " -> full_score:" << full_score << " (M_{i-1,j'}=" << score << ")"
                    //                      << std::endl;

                    if (full_score < full_minimum)
                    {
                        minimum = score;
                        full_minimum = full_score;
                        trace[i][j] = {i - 1, j_prime};
                        ll_matrix[i][j] = ll_kmers;
                    }
                }

                matrix[i][j] = minimum;
            }
        }
    }

    //!\brief Backtracks the trace matrix and writes the resulting binning into the output file.
    size_t backtracking(std::vector<std::vector<size_t>> const & matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> const & trace)
    {
        assert(data != nullptr);
        assert(data->output_buffer != nullptr);
        assert(data->header_buffer != nullptr);

        bool const high = data->previous.empty();

        // backtracking
        size_t trace_i = num_technical_bins - 1;
        int trace_j = num_user_bins - 1;
        // std::cout << "optimum: " << matrix[trace_i][trace_j] << std::endl;
        // std::cout << std::endl;

        if (data->output_buffer->tellp() == 0) // beginning of the file
        {
            if (config.debug)
                *data->output_buffer << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX" << std::endl;
            else
                *data->output_buffer << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS" << std::endl;
        }

        size_t high_level_max_id{};
        size_t high_level_max_size{};

        size_t bin_id{};
        size_t const optimal_score{matrix[trace_i][trace_j]};
        double correction{};

        // TODO std::to_chars after https://github.com/seqan/product_backlog/issues/396
        auto to_string_with_precision = [](double const value)
        {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(2) << value;
            return stream.str();
        };

        while (trace_j >= 0)
        {
            // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
            size_t next_i = trace[trace_i][trace_j].first;
            size_t next_j = trace[trace_i][trace_j].second;

            size_t kmer_count = data->kmer_counts[trace_j];
            size_t number_of_bins = (trace_i - next_i);

            correction = data->previous.empty() ? data->fp_correction[std::max<size_t>(1u, number_of_bins)] : 1u;

            if (trace_j == 0)
            {
                // we only arrive here if the bin wasn't merged with some before so it is safe to assume
                // that the bin was split (even if only into 1 bin).

                ++trace_i; // because we want the length not the index
                size_t const kmer_count = data->kmer_counts[0];
                size_t const average_bin_size = kmer_count / trace_i;

                // add query cost for determination of best t_max
                total_query_cost += (data->previous.cost + interpolated_cost) * kmer_count;

                // add split bin to ibf statistics
                data->stats->emplace_back(hibf_statistics::bin_kind::split, average_bin_size, 1ul, trace_i);

                *data->output_buffer << data->filenames[0] << '\t'
                                     << data->previous.bin_indices  << (high ? "" : ";") << bin_id << '\t'
                                     << data->previous.num_of_bins  << (high ? "" : ";") << trace_i;

                if (config.debug)
                {
                    *data->output_buffer << '\t'
                                         << data->previous.estimated_sizes << (high ? "" : ";") << average_bin_size << '\t'
                                         << data->previous.optimal_score << (high ? "" : ";") << optimal_score << '\t'
                                         << data->previous.correction << (high ? "" : ";") << correction << '\t'
                                         << data->previous.tmax << (high ? "" : ";") << num_technical_bins;
                }

                *data->output_buffer << '\n';

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
                libf_data.output_buffer = data->output_buffer;
                libf_data.header_buffer = data->header_buffer;

                libf_data.kmer_counts = {kmer_count};
                libf_data.filenames = {data->filenames[trace_j]};
                size_t num_contained_ubs = 1;
                size_t const j = trace_j;

                // std::cout << "merged [" << trace_j;
                while (trace_j > 0 && next_i == trace_i)
                {
                    trace_i = next_i; // unnecessary?
                    --trace_j;
                    kmer_count += data->kmer_counts[trace_j];
                    libf_data.kmer_counts.push_back(data->kmer_counts[trace_j]);
                    libf_data.filenames.push_back(data->filenames[trace_j]);
                    ++num_contained_ubs;
                    next_i = trace[trace_i][trace_j].first;
                    // std::cout << "," << trace_j;
                }
                assert(trace_j == 0 || trace_i - next_i == 1);
                assert(kmer_count == std::accumulate(libf_data.kmer_counts.begin(), libf_data.kmer_counts.end(), 0u));

                ++number_of_bins;
                trace_i = next_i;
                --trace_j;

                libf_data.previous = data->previous;
                libf_data.previous.bin_indices += (high ? "" : ";") + std::to_string(bin_id);
                libf_data.previous.num_of_bins  += (high ? "" : ";") + std::string{"1"};
                libf_data.previous.cost += interpolated_cost;
                if (config.debug)
                {
                    libf_data.previous.estimated_sizes += (high ? "" : ";") + std::to_string(kmer_count);
                    libf_data.previous.optimal_score += (high ? "" : ";") + std::to_string(optimal_score);
                    libf_data.previous.correction += (high ? "" : ";") + to_string_with_precision(correction);
                    libf_data.previous.tmax += (high ? "" : ";") + std::to_string(num_technical_bins);
                }

                std::string const merged_ibf_name{std::string{merged_bin_prefix} + "_" + libf_data.previous.bin_indices};

                // add merged bin to ibf statistics
                uint64_t const cardinality = config.estimate_union ? data->union_estimates[j][trace_j + 1] : kmer_count;
                hibf_statistics::bin & bin_stats = data->stats->emplace_back(hibf_statistics::bin_kind::merged,
                                                   cardinality, num_contained_ubs, 1ul);

                // now do the binning for the low-level IBF:
                size_t merged_max_bin_id;
                if (libf_data.kmer_counts.size() > config.t_max)
                {
                    // recursively call hierarchical binning if there are still too many UBs
                    libf_data.stats = &bin_stats.child_level;
                    auto const && [bin_id, cost] = hierarchical_binning{libf_data, config}.execute();
                    merged_max_bin_id = bin_id;
                    total_query_cost += cost;
                }
                else
                {
                    // use simple binning to distribute remaining UBs
                    libf_data.stats = &bin_stats.child_level;
                    simple_binning algo{libf_data, 0, config.debug};
                    merged_max_bin_id = algo.execute();
                    total_query_cost += (data->previous.cost + interpolated_cost
                                        + ibf_query_cost::interpolated(algo.get_num_technical_bins()))
                                        * kmer_count;
                }
                *data->header_buffer << "#" << merged_ibf_name << " max_bin_id:" << merged_max_bin_id << '\n';

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
                libf_data.output_buffer = data->output_buffer;
                libf_data.header_buffer = data->header_buffer;

                libf_data.kmer_counts = {kmer_count};
                libf_data.filenames = {data->filenames[trace_j]};
                size_t num_contained_ubs = 1;
                size_t const j = trace_j;

                // std::cout << "merged [" << trace_j;
                --trace_j;
                while (static_cast<size_t>(trace_j) != next_j)
                {
                    kmer_count += data->kmer_counts[trace_j];
                    libf_data.kmer_counts.push_back(data->kmer_counts[trace_j]);
                    libf_data.filenames.push_back(data->filenames[trace_j]);
                    ++num_contained_ubs;
                    // std::cout << "," << trace_j;
                    --trace_j;
                }
                trace_i = next_i;
                trace_j = next_j; // unneccessary?

                libf_data.previous = data->previous;
                libf_data.previous.bin_indices += (high ? "" : ";") + std::to_string(bin_id);
                libf_data.previous.num_of_bins  += (high ? "" : ";") + std::string{"1"};
                libf_data.previous.cost += interpolated_cost;
                if (config.debug)
                {
                    libf_data.previous.estimated_sizes += (high ? "" : ";") + std::to_string(kmer_count);
                    libf_data.previous.optimal_score += (high ? "" : ";") + std::to_string(optimal_score);
                    libf_data.previous.correction += (high ? "" : ";") + to_string_with_precision(correction);
                    libf_data.previous.tmax += (high ? "" : ";") + std::to_string(num_technical_bins);
                }
                std::string const merged_ibf_name{std::string{merged_bin_prefix} + "_" + libf_data.previous.bin_indices};

                // add merged bin to ibf statistics
                uint64_t const cardinality = config.estimate_union ? data->union_estimates[j][trace_j + 1] : kmer_count;
                hibf_statistics::bin & bin_stats = data->stats->emplace_back(hibf_statistics::bin_kind::merged,
                                                cardinality, num_contained_ubs, 1ul);

                size_t merged_max_bin_id;
                // now do the binning for the low-level IBF:
                if (libf_data.kmer_counts.size() > config.t_max)
                {
                    libf_data.stats = &bin_stats.child_level;
                    auto const && [bin_id, cost] = hierarchical_binning{libf_data, config}.execute();
                    merged_max_bin_id = bin_id;
                    total_query_cost += cost;
                }
                else
                {
                    libf_data.stats = &bin_stats.child_level;
                    simple_binning algo{libf_data, 0, config.debug};
                    merged_max_bin_id = algo.execute();
                    total_query_cost += (data->previous.cost + interpolated_cost
                                        + ibf_query_cost::interpolated(algo.get_num_technical_bins()))
                                        * kmer_count;
                }
                *data->header_buffer << "#" << merged_ibf_name << " max_bin_id:" << merged_max_bin_id << '\n';

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

                *data->output_buffer << data->filenames[trace_j] << '\t'
                                     << data->previous.bin_indices  << (high ? "" : ";") << bin_id << '\t'
                                     << data->previous.num_of_bins  << (high ? "" : ";") << number_of_bins;

                total_query_cost += (data->previous.cost + interpolated_cost) * kmer_count;

                // add split bin to ibf statistics
                data->stats->emplace_back(hibf_statistics::bin_kind::split, kmer_count_per_bin, 1ul, number_of_bins);

                if (config.debug)
                {
                    *data->output_buffer << '\t'
                                         << data->previous.estimated_sizes << (high ? "" : ";") << kmer_count_per_bin << '\t'
                                         << data->previous.optimal_score << (high ? "" : ";") << optimal_score << '\t'
                                         << data->previous.correction << (high ? "" : ";") << correction << '\t'
                                         << data->previous.tmax << (high ? "" : ";") << num_technical_bins;
                }

                *data->output_buffer << '\n';

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
