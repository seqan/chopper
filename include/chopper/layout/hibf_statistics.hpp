#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include <chopper/helper.hpp>
#include <chopper/configuration.hpp>
#include <chopper/layout/ibf_query_cost.hpp>
#include <chopper/sketch/hyperloglog.hpp>
#include <chopper/sketch/user_bin_sequence.hpp>

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

        //!\brief A pointer to the filenames of the user input sequences.
        std::vector<std::string> filenames;
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
        compute_total_query_cost(top_level_ibf);

        gather_statistics(top_level_ibf, 0);

        expected_HIBF_query_cost = total_query_cost / total_kmer_count;
    }

    //!\brief Prints a column names of the summary to the command line.
    static void print_header(bool const verbose = true)
    {
        // print column names explanation in header
        std::cout << "## ### Notation ###\n"
                  << "## X-IBF = An IBF with X number of bins.\n"
                  << "## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.\n";

        std::cout << "## ### Column Description ###\n"
                     "## tmax : The maximum number of technical bin on each level\n"
                     "## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF\n"
                     "## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF\n"
                     "## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF\n"
                     "## (l*m)_tmax : Computed by l_tmax * m_tmax\n"
                     "## size : The expected total size of an tmax-HIBF\n"
                  << ((verbose) ? "## uncorr_size : The expected size of an tmax-HIBF without FPR correction\n" : "");

        // print column names
        std::cout << "# tmax"     << '\t'
                  << "c_tmax"     << '\t'
                  << "l_tmax"     << '\t'
                  << "m_tmax"     << '\t'
                  << "(l*m)_tmax" << '\t'
                  << "size";

        if (verbose) // uncorrected size and add level statistics
        {
            std::cout << '\t'
                      << "uncorr_size"         << '\t'
                      << "level"               << '\t'
                      << "num_ibfs"            << '\t'
                      << "level_size"          << '\t'
                      << "level_size_no_corr"  << '\t'
                      << "total_num_tbs"       << '\t'
                      << "avg_num_tbs"         << '\t'
                      << "split_tb_percentage" << '\t'
                      << "max_split_tb"        << '\t'
                      << "avg_split_tb"        << '\t'
                      << "max_factor"          << '\t'
                      << "avg_factor";
        }

        std::cout << '\n';
    }

    //!\brief Prints a tab-separated summary of the statistics of this HIBF to the command line.
    void print_summary(size_t & t_max_64_memory, bool const verbose = true)
    {
        if (summaries.empty())
            finalize();

        if (t_max_64_memory == 0)
            t_max_64_memory = total_hibf_size_in_byte();

        double const relative_memory_size = total_hibf_size_in_byte() /
                                            static_cast<double>(t_max_64_memory);
        double const query_time_memory_usage_prod = expected_HIBF_query_cost * relative_memory_size;

        std::cout << std::fixed << std::setprecision(2);

        std::string level_str, num_ibfs_str, level_size_str, level_size_no_corr_str, total_num_tbs_str, avg_num_tbs_str,
                    split_tb_percentage_str, max_split_tb_str, avg_split_tb_str, max_factor_str, avg_factor_str;

        size_t total_size{};
        size_t total_size_no_corr{};

        // go through each level and collect and output the statistics
        auto to_string_with_precision = [] (auto num)
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << num;
            return ss.str();
        };

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

            level_str               += ":" + to_string_with_precision(level);
            num_ibfs_str            += ":" + to_string_with_precision(s.num_ibfs);
            level_size_str          += ":" + to_formatted_BF_size(level_size);
            level_size_no_corr_str  += ":" + to_formatted_BF_size(level_size_no_corr);
            total_num_tbs_str       += ":" + to_string_with_precision(total_num_tbs);
            avg_num_tbs_str         += ":" + to_string_with_precision(total_num_tbs / s.num_ibfs);
            split_tb_percentage_str += ":" + to_string_with_precision(split_tb_percentage);

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

                max_split_tb_str += ":" + to_string_with_precision(max_split_bin_span);
                avg_split_tb_str += ":" + to_string_with_precision(avg_split_bin);
                max_factor_str   += ":" + to_string_with_precision((*fp_correction)[max_split_bin_span]);
                avg_factor_str   += ":" + to_string_with_precision(avg_factor);
            }
            else
            {
                max_split_tb_str += ":-";
                avg_split_tb_str += ":-";
                max_factor_str   += ":-";
                avg_factor_str   += ":-";
            }
        }

        std::cout << std::fixed << std::setprecision(2);

        std::cout /*        tmax */ << config.tmax << '\t'
                  /*      c_tmax */ << chopper::layout::ibf_query_cost::interpolated(config.tmax, config.false_positive_rate) << '\t'
                  /*      l_tmax */ << expected_HIBF_query_cost << '\t' /*relative to a 64 bin IBF*/
                  /*      m_tmax */ << relative_memory_size << '\t' /*relative to the 64 T_Max HIBF*/
                  /*   (l*m)tmax */ << query_time_memory_usage_prod << '\t'
                  /*  corr. size */ << to_formatted_BF_size(total_size) << ((verbose) ? '\t' : '\n');

        if (verbose)
        {
            // uncorrected FPR
            std::cout /*uncorr. size */ << to_formatted_BF_size(total_size_no_corr) << '\t';

            // per level statistics:
            std::cout /* level               */ << level_str << '\t'
                      /* num_ibfs            */ << num_ibfs_str << '\t'
                      /* level_size          */ << level_size_str << '\t'
                      /* level_size_no_corr  */ << level_size_no_corr_str << '\t'
                      /* total_num_tbs       */ << total_num_tbs_str << '\t'
                      /* avg_num_tbs         */ << avg_num_tbs_str << '\t'
                      /* split_tb_percentage */ << split_tb_percentage_str << '\t'
                      /* max_split_tb        */ << max_split_tb_str << '\t'
                      /* avg_split_tb        */ << avg_split_tb_str << '\t'
                      /* max_factor          */ << max_factor_str << '\t'
                      /* avg_factor          */ << avg_factor_str << '\n';
        }
    }

    //!\brief Return the total corrected size of the HIBF in bytes
    size_t total_hibf_size_in_byte()
    {
        if (summaries.empty())
            finalize();

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

    //!\brief Computes the estimated query cost
    void compute_total_query_cost(level & curr_level)
    {
        // Compute number of technical bins in current level (<= tmax)
        size_t number_of_tbs{0};
        size_t level_kmer_count{0};
        size_t index{0};
        std::vector<size_t> merged_bin_indices{};
        std::vector<sketch::hyperloglog> merged_bin_sketches{};

        for (bin const & current_bin : curr_level.bins)
        {
            if (current_bin.kind == bin_kind::merged)
            {
                ++number_of_tbs;
                merged_bin_indices.push_back(index);

                if (config.estimate_union)
                {
                    // compute merged_bin_sketch
                    assert(!current_bin.child_level.filenames.empty());
                    std::vector<sketch::hyperloglog> sketches;
                    sketch::user_bin_sequence::read_hll_files_into(config.sketch_directory,
                                                                   current_bin.child_level.filenames,
                                                                   sketches);
                    sketch::hyperloglog hll = sketches[0];
                    for (size_t i = 1; i < sketches.size(); ++i)
                        hll.merge(sketches[i]);

                    merged_bin_sketches.push_back(std::move(hll));
                }
            }
            else if (current_bin.kind == bin_kind::split) // bin_kind::split
            {
                number_of_tbs += current_bin.num_spanning_tbs;
                level_kmer_count += current_bin.cardinality;
            }
            ++index;
        }
        assert(number_of_tbs <= config.tmax);

        // Add cost of querying the current IBF
        // (how costly is querying number_of_tbs (e.g. 128 tbs) compared to 64 tbs given the current FPR)
        curr_level.current_query_cost += ibf_query_cost::interpolated(number_of_tbs, config.false_positive_rate);

        // Add costs of querying the HIBF for each kmer in this level.
        total_query_cost += curr_level.current_query_cost * level_kmer_count;

        // update query cost of all merged bins
        for (size_t i = 0; i < merged_bin_indices.size(); ++i)
        {
            auto & current_bin = curr_level.bins[merged_bin_indices[i]];

            // Pass on cost of querying the current level
            current_bin.child_level.current_query_cost = curr_level.current_query_cost;

            // If merged bins share kmers, we need to penalize this
            // because querying a kmer will result in multi level look-ups.
            if (config.estimate_union)
            {
                double const current_estimate = merged_bin_sketches[i].estimate();

                for (size_t j = i + 1; j < merged_bin_indices.size(); ++j)
                {
                    sketch::hyperloglog tmp = merged_bin_sketches[i]; // copy needed, s.t. current is not modified
                    double union_estimate = tmp.merge_and_estimate_SIMD(merged_bin_sketches[j]);
                    // Jaccard distance estimate
                    double distance = 2.0 - (current_estimate + merged_bin_sketches[j].estimate()) / union_estimate;
                    // Since the sizes are estimates, the distance might be slighlty above 1.0 or below 0.0
                    // but we need to avoid nagetive numbers
                    distance = std::min(std::max(distance, 0.0), 1.0);

                    current_bin.child_level.current_query_cost += (1.0 - distance);
                }
            }
        }

        // call function recursively for each merged bin
        for (size_t i : merged_bin_indices)
            compute_total_query_cost(curr_level.bins[i].child_level);
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
