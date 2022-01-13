#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include <chopper/helper.hpp>
#include <chopper/layout/configuration.hpp>
#include <chopper/layout/ibf_query_cost.hpp>

namespace chopper::layout
{

class hibf_statistics
{
public:
    hibf_statistics() = default; //!< Defaulted.
    hibf_statistics(hibf_statistics const & b) = default; //!< Defaulted.
    hibf_statistics & operator=(hibf_statistics const &) = default; //!< Defaulted.
    hibf_statistics(hibf_statistics && b) = default; //!< Defaulted.
    hibf_statistics & operator=(hibf_statistics &&) = default; //!< Defaulted.
    ~hibf_statistics() = default; //!< Defaulted.

    /*!\brief Construct an empty HIBF with an empty top level IBF
     * \param[in] config_ User configuration for the HIBF.
     * \param[in] fp_correction_ The false positive correction factors to use for the statistics.
     * \param[in] kmer_counts The original user bin weights (kmer counts).
     */
    hibf_statistics(configuration const & config_,
                    std::vector<double> const & fp_correction_,
                    std::vector<size_t> const & kmer_counts) :
        config{config_},
        fp_correction{&fp_correction_},
        total_kmer_count{std::accumulate(kmer_counts.begin(), kmer_counts.end(), size_t{})}
    {}

    struct bin; // forward declaration

    //!\brief A representation of an IBF level that gathers information about bins in an IBF.
    struct level
    {
        //!\brief The bins of the current IBF level. May be split or merged bins.
        std::vector<bin> bins;

        //!\brief The query cost to arrive at this IBF (updated before backtracking respective DP).
        double current_query_cost{0.0};
    };

    //!\brief The kind of bin that is stored.
    enum class bin_kind
    {
        split, //!< A single user bin, split into 1 or more bins (even though 1 is not technically split).
        merged //!< Multiple user bins are merged into a single technical bin.
    };

    //!\brief Represents a (set) of user bins (see ibf_statistics::bin_kind).
    class bin
    {
    public:
        bin_kind const kind; //!< Either a split or merged bin.
        size_t const cardinality; //!< The size/weight of the bin (either a kmer count or hll sketch estimation).
        size_t const num_contained_ubs; //!< [MERGED] How many UBs are merged within this TB.
        size_t const num_spanning_tbs; //!< [SPLIT] How many TBs are used for this sindle UB.

        level child_level; //!< [MERGED] The lower level ibf statistics.

        bin() = default; //!< Defaulted.
        bin(bin const & b) = default; //!< Defaulted.
        bin & operator=(bin const &) = default; //!< Defaulted.
        bin(bin && b) = default; //!< Defaulted.
        bin & operator=(bin &&) = default; //!< Defaulted.
        ~bin() = default; //!< Defaulted.

        bin(bin_kind const kind_,
            size_t const card,
            size_t const contained_ubs,
            size_t const spanning_tbs) :
            kind{kind_},
            cardinality{card},
            num_contained_ubs{contained_ubs},
            num_spanning_tbs{spanning_tbs}
        {
            assert((kind == bin_kind::split  && num_contained_ubs == 1u) ||
                   (kind == bin_kind::merged && num_spanning_tbs  == 1u));
        }
    };

    //!\brief Gather all statistics to have all members ready.
    void finalize()
    {
        gather_statistics(top_level_ibf, 0);

        expected_HIBF_query_cost = total_query_cost / total_kmer_count;
    }

    //!\brief Prints a tab-separated summary of the statistics of this HIBF to the command line.
    void print_summary(size_t & t_max_64_memory)
    {
        if (summaries.empty())
            finalize();

        if (t_max_64_memory == 0)
            t_max_64_memory = total_hibf_size_in_byte();

        double const relative_memory_size = total_hibf_size_in_byte() /
                                            static_cast<double>(t_max_64_memory);
        double const query_time_memory_usage_prod = expected_HIBF_query_cost * relative_memory_size;

        std::cout << std::fixed << std::setprecision(2);

        std::cout << "#T_Max:" << config.tmax << '\n'
                  << "#C_{T_Max}:" << chopper::layout::ibf_query_cost::interpolated(config.tmax, config.false_positive_rate) << '\n'
                  << "#relative expected HIBF query time cost (l):" << expected_HIBF_query_cost << '\n' /*relative to a 64 bin IBF*/
                  << "#relative HIBF memory usage (m):" << relative_memory_size << '\n' /*relative to the 64 T_Max HIBF*/
                  << "#l*m:" << query_time_memory_usage_prod << '\n';

        // print column names
        std::cout << "level\tnum_ibfs\tlevel_size\tlevel_size_no_corr\ttotal_num_tbs"
                     "\tavg_num_tbs\tsplit_tb_percentage\tmax_split_tb\tavg_split_tb\tmax_factor\tavg_factor\n";

        size_t total_size{};
        size_t total_size_no_corr{};

        // go through each level and collect and output the statistics
        for (auto const & [level, s] : summaries)
        {
            size_t const level_size = std::reduce(s.ibf_mem_size.begin(), s.ibf_mem_size.end());
            size_t const level_size_no_corr = std::reduce(s.ibf_mem_size_no_corr.begin(), s.ibf_mem_size_no_corr.end());

            total_size += level_size;
            total_size_no_corr += level_size_no_corr;

            size_t const total_num_tbs = std::reduce(s.num_tbs.begin(), s.num_tbs.end());

            size_t const total_num_split_tbs = std::reduce(s.num_split_tbs.begin(), s.num_split_tbs.end());
            double const split_tb_percentage = 100.0 * static_cast<double>(total_num_split_tbs) / total_num_tbs;

            size_t const max_split_bin_span = *std::max_element(s.max_split_tb_span.begin(), s.max_split_tb_span.end());

            std::cout << level << '\t'
                      << s.num_ibfs << '\t'
                      << to_formatted_BF_size(level_size) << '\t'
                      << to_formatted_BF_size(level_size_no_corr) << '\t'
                      << total_num_tbs << '\t'
                      << total_num_tbs / s.num_ibfs << '\t'
                      << split_tb_percentage << '\t';

            // if there are no split bins on this level, the following statistics don't make sense
            if (max_split_bin_span != 0)
            {
                size_t const total_num_split_ubs = std::reduce(s.num_split_ubs.begin(), s.num_split_ubs.end());
                double const avg_split_bin = static_cast<double>(total_num_split_tbs)
                                           / static_cast<double>(total_num_split_ubs);
                size_t const total_split_tb_kmers = std::reduce(s.split_tb_kmers.begin(), s.split_tb_kmers.end());
                double const avg_factor = static_cast<double>(std::reduce(s.split_tb_corr_kmers.begin(),
                                                                          s.split_tb_corr_kmers.end()))
                                        / static_cast<double>(total_split_tb_kmers);

                std::cout << max_split_bin_span << '\t'
                          << avg_split_bin << '\t'
                          << (*fp_correction)[max_split_bin_span] << '\t'
                          << avg_factor << '\n';
            }
            else
            {
                std::cout << "-\t-\t-\t-\n";
            }
        }

        std::cout << "#Total HIBF size:" << to_formatted_BF_size(total_size) << '\n'
                  << "#Total HIBF size no correction:" << to_formatted_BF_size(total_size_no_corr) << "\n\n";
    }

    //!\brief Return the total corrected size of the HIBF in bytes
    size_t total_hibf_size_in_byte()
    {
        if (summaries.empty())
            gather_statistics(top_level_ibf, 0);

        size_t total_size{};

        // go through each level and collect the memory sizes
        for (auto const & [level, summary] : summaries)
        {
            (void) level;

            total_size += std::reduce(summary.ibf_mem_size.begin(), summary.ibf_mem_size.end());
        }

        return compute_bin_size(total_size) / 8;
    }

    //!\brief The top level IBF of this HIBF, often starting point for recursions.
    level top_level_ibf;

    //!\brief The estimated query cost of every single kmer in this HIBF.
    double total_query_cost{0.0};

    //!\brief The estimated query cost relative to the total k-mer count in the data set.
    double expected_HIBF_query_cost{0.0};
private:
    //!\brief Copy of the user configuration for this HIBF.
    configuration const config{};

    //!\brief The false positive correction factors to use for the statistics.
    std::vector<double> const * const fp_correction{nullptr};

    //!\brief The original kmer count of all user bins.
    size_t const total_kmer_count{};

    //!\brief Statistics for all IBFs on a certain level of the HIBF.
    struct level_summary
    {
        size_t num_ibfs{};

        std::vector<size_t> num_tbs{};
        std::vector<size_t> num_ubs{};

        std::vector<size_t> num_split_tbs{};
        std::vector<size_t> num_merged_tbs{};

        std::vector<size_t> num_split_ubs{};
        std::vector<size_t> num_merged_ubs{};

        std::vector<size_t> max_split_tb_span{};
        std::vector<size_t> split_tb_corr_kmers{};
        std::vector<size_t> split_tb_kmers{};

        std::vector<size_t> max_ubs_in_merged{};

        std::vector<size_t> ibf_mem_size{};
        std::vector<size_t> ibf_mem_size_no_corr{};
    };

    //!\brief The gathered summary of statistics for each level of this HIBF.
    std::map<size_t, level_summary> summaries;

    /*!\brief Computes the bin size in bits.
    *
    * -NUM_ELEM*HASHES
    * ----------------------  = SIZE
    * LN(1-FPR^(1/HASHES))
    *
    * -NUM_ELEMS*HASHES
    * -----------------------
    * LN(1 - e^(LN(FPR) / HASHES) )
    */
    size_t compute_bin_size(size_t const number_of_kmers_to_be_stored) const
    {
        return std::ceil( - static_cast<double>(number_of_kmers_to_be_stored * config.num_hash_functions) /
               std::log(1 - std::exp(std::log(config.false_positive_rate) / config.num_hash_functions)));
    }

    /*!\brief Compute the Bloom Filter size from `number_of_kmers_to_be_stored` and
     *        return it as a formatted string with the appropriate unit.
     * \param[in] number_of_kmers_to_be_stored
     */
    std::string to_formatted_BF_size(size_t const number_of_kmers_to_be_stored) const
    {
        size_t const size_in_bytes = compute_bin_size(number_of_kmers_to_be_stored) / 8;
        return byte_size_to_formatted_str(size_in_bytes);
    }

    /*!\brief Recursively gather all the statistics from the bins.
     * \param[in] curr_level The current IBF from which the statistics will be extracted.
     * \param[in] level_summary_index The index of `curr_level` in `summeries`.
     */
    void gather_statistics(level const & curr_level, size_t const level_summary_index)
    {
        level_summary & summary = summaries[level_summary_index];
        summary.num_ibfs += 1;

        size_t max_cardinality{}, max_cardinality_no_corr{}, num_tbs{}, num_ubs{}, num_split_tbs{},
               num_merged_tbs{}, num_split_ubs{}, num_merged_ubs{}, max_split_tb_span{},
               split_tb_kmers{}, max_ubs_in_merged{}, split_tb_corr_kmers{};

        for (bin const & current_bin : curr_level.bins)
        {
            size_t const cardinality_per_split_bin = (current_bin.cardinality + current_bin.num_spanning_tbs - 1) /
                                                     current_bin.num_spanning_tbs; // round up
            size_t const corrected_cardinality = std::ceil(cardinality_per_split_bin *
                                                           (*fp_correction)[current_bin.num_spanning_tbs]);
            max_cardinality = std::max(max_cardinality, corrected_cardinality);
            max_cardinality_no_corr = std::max(max_cardinality_no_corr, cardinality_per_split_bin);

            num_tbs += current_bin.num_spanning_tbs;
            num_ubs += current_bin.num_contained_ubs;

            if (current_bin.kind == bin_kind::split)
            {
                num_split_tbs += current_bin.num_spanning_tbs;
                num_split_ubs += 1;
                split_tb_corr_kmers += corrected_cardinality * current_bin.num_spanning_tbs;
                split_tb_kmers += cardinality_per_split_bin * current_bin.num_spanning_tbs;
                max_split_tb_span = std::max(max_split_tb_span, current_bin.num_spanning_tbs);
                total_query_cost += curr_level.current_query_cost * current_bin.cardinality;
            }
            else
            {
                num_merged_tbs += 1;
                num_merged_ubs += current_bin.num_contained_ubs;
                max_ubs_in_merged = std::max(max_ubs_in_merged, current_bin.num_contained_ubs);

                gather_statistics(current_bin.child_level, level_summary_index + 1);
            }
        }

        summary.num_tbs.push_back(num_tbs);
        summary.num_ubs.push_back(num_ubs);

        summary.num_split_tbs.push_back(num_split_tbs);
        summary.num_merged_tbs.push_back(num_merged_tbs);

        summary.num_split_ubs.push_back(num_split_ubs);
        summary.num_merged_ubs.push_back(num_merged_ubs);

        summary.max_split_tb_span.push_back(max_split_tb_span);
        summary.split_tb_corr_kmers.push_back(split_tb_corr_kmers);
        summary.split_tb_kmers.push_back(split_tb_kmers);

        summary.max_ubs_in_merged.push_back(max_ubs_in_merged);

        summary.ibf_mem_size.push_back(max_cardinality * num_tbs);
        summary.ibf_mem_size_no_corr.push_back(max_cardinality_no_corr * num_tbs);
    }
};

} // namespace chopper::layout
