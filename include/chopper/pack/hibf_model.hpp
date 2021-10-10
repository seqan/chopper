#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include <chopper/helper.hpp>
#include <chopper/pack/pack_config.hpp>

#include <robin_hood.h>

class hibf_model
{
public:
    struct bin;
    using ibf = std::vector<bin>;

    enum class bin_kind
    {
        split, merged
    };

    /*!\brief Represents either a single user bin split into multiple technical bins
     * or multiple user bins merged into a single technical bin.
     */
    struct bin
    {
        bin_kind const kind;
        size_t const cardinality;
        size_t const num_contained_ubs;
        size_t const num_spanning_tbs;
        ibf child_ibf;

        bin() = delete; //!< Deleted. Enforce user supplied values for member initialization.
        bin(bin const & b) = default; //!< Defaulted.
        bin & operator=(bin const &) = default; //!< Defaulted.
        bin(bin && b) = default; //!< Defaulted.
        bin & operator=(bin &&) = default; //!< Defaulted.
        ~bin() = default; //!< Defaulted.

        bin(bin_kind const kind_,
            size_t const cardinality_,
            size_t const num_contained_ubs_,
            size_t const num_spanning_tbs_) : 
            kind{kind_},
            cardinality{cardinality_},
            num_contained_ubs{num_contained_ubs_},
            num_spanning_tbs{num_spanning_tbs_}
        { 
            assert((kind == bin_kind::split  && num_contained_ubs == 1ul) || 
                   (kind == bin_kind::merged && num_spanning_tbs  == 1ul));
        }
    };

    //!\brief Statistics for all IBFs on a certain level of the HIBF model.
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

private:
    //!\brief User configuration for this HIBF model.
    pack_config const & config;
    //!\brief The false positive correction factors to use for the statistics.
    std::vector<double> const & fp_correction;

    //!\brief The top level IBF of this HIBF model, often starting point for recursions.
    ibf top_level_ibf;
    //!\brief The gathered summary of statistics for each level of this HIBF.
    robin_hood::unordered_map<size_t, level_summary> summary;

public:
    hibf_model() = delete; //!< Deleted. Would be malformed.
    hibf_model(hibf_model const & b) = delete; //!< Defaulted.
    hibf_model & operator=(hibf_model const &) = delete; //!< Defaulted.
    hibf_model(hibf_model && b) = default; //!< Defaulted.
    hibf_model & operator=(hibf_model &&) = default; //!< Defaulted.
    ~hibf_model() = default; //!< Defaulted.

    /*!\brief Construct an empty HIBF model with an empty top level IBF
     * \param[in] config_ User configuration for the HIBF model.
     * \param[in] fp_correction_ The false positive correction factors to use for the statistics.
     */
    hibf_model(pack_config const & config_, std::vector<double> const & fp_correction_) :
        config{config_},
        fp_correction{fp_correction_}
    {}

    //!\brief Returns a reference to the top level IBF of this HIBF model.
    [[nodiscard]] ibf & get_top_level_ibf()
    {
        return top_level_ibf;
    }

    //!\brief Prints a tab-indented summary of the statistics of this HIBF model to the command line
    void print_summary()
    {
        if (summary.empty())
        {
            gather_statistics(top_level_ibf, 0);
        }

        std::cout << "\n\tStatistics summary:\n\tlevel\tnum_ibfs\tlevel_size\tlevel_size_no_corr\ttotal_num_tbs"
                     "\tavg_num_tbs\tsplit_tb_percentage\tmax_split_tb\tavg_split_tb\tmax_factor\tavg_factor\n";

        size_t total_size{};
        size_t total_size_no_corr{};

        // go through each level and collect and output the statistics 
        for (auto const & [level, s] : summary)
        {
            size_t const level_size = std::reduce(s.ibf_mem_size.begin(), s.ibf_mem_size.end());
            size_t const level_size_no_corr = std::reduce(s.ibf_mem_size_no_corr.begin(), s.ibf_mem_size_no_corr.end());
            
            total_size += level_size;
            total_size_no_corr += level_size_no_corr;

            size_t const total_num_tbs = std::reduce(s.num_tbs.begin(), s.num_tbs.end());

            size_t const total_num_split_tbs = std::reduce(s.num_split_tbs.begin(), s.num_split_tbs.end());
            double const split_tb_percentage = 100 * (double)total_num_split_tbs / (double)total_num_tbs;
            
            size_t const max_split_bin_span = *std::max_element(s.max_split_tb_span.begin(), s.max_split_tb_span.end());

            std::cout << '\t' << level << '\t'
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
                          << fp_correction[max_split_bin_span] << '\t'
                          << avg_factor << '\n';
            }
            else 
            {
                std::cout << "-\t-\t-\t-\n";
            }
        }

        std::cout << "\tTotal HIBF size: " << to_formatted_BF_size(total_size)
                  << "\n\tTotal HIBF size no correction: " << to_formatted_BF_size(total_size_no_corr) << "\n\n";
    }

private:
    /*!\brief Compute the Bloom Filter size from `number_of_kmers_to_be_stored` and 
     * return it as a formatted string with the appropriate unit.
     * \param[in] number_of_kmers_to_be_stored
     */
    std::string to_formatted_BF_size(size_t const number_of_kmers_to_be_stored) const
    {
        size_t const size = compute_bin_size(config, number_of_kmers_to_be_stored) / 8;
        return byte_size_to_formatted_str(size);
    }

    /*!\brief Recursively gather all the statistics from the bins in the model.
     * \param[in] curr_ibf The current IBF from which the statistics will be extracted.
     * \param[in] level The level of `curr_ibf` in the HIBF.
     */
    void gather_statistics(ibf const & curr_ibf, size_t const level)
    {
        level_summary & s = summary[level];
        s.num_ibfs += 1;
        
        size_t max_cardinality{}, max_cardinality_no_corr{}, num_tbs{}, num_ubs{}, num_split_tbs{}, 
               num_merged_tbs{}, num_split_ubs{}, num_merged_ubs{}, max_split_tb_span{},
               split_tb_kmers{}, max_ubs_in_merged{}, split_tb_corr_kmers{};

        for (bin const & b : curr_ibf)
        {
            size_t const corrected_cardinality = std::ceil(b.cardinality * fp_correction[b.num_spanning_tbs]);
            max_cardinality = std::max(max_cardinality, corrected_cardinality);
            max_cardinality_no_corr = std::max(max_cardinality_no_corr, b.cardinality);

            num_tbs += b.num_spanning_tbs;
            num_ubs += b.num_contained_ubs;

            if (b.kind == bin_kind::split)
            {
                num_split_tbs += b.num_spanning_tbs;
                num_split_ubs += 1;
                split_tb_corr_kmers += corrected_cardinality * b.num_spanning_tbs;
                split_tb_kmers += b.cardinality * b.num_spanning_tbs;
                max_split_tb_span = std::max(max_split_tb_span, b.num_spanning_tbs);
            }
            else 
            {
                num_merged_tbs += 1;
                num_merged_ubs += b.num_contained_ubs;
                max_ubs_in_merged = std::max(max_ubs_in_merged, b.num_contained_ubs);

                gather_statistics(b.child_ibf, level + 1);
            }
        }

        s.num_tbs.push_back(num_tbs);
        s.num_ubs.push_back(num_ubs);

        s.num_split_tbs.push_back(num_split_tbs);
        s.num_merged_tbs.push_back(num_merged_tbs);

        s.num_split_ubs.push_back(num_split_ubs);
        s.num_merged_ubs.push_back(num_merged_ubs);

        s.max_split_tb_span.push_back(max_split_tb_span);
        s.split_tb_corr_kmers.push_back(split_tb_corr_kmers);
        s.split_tb_kmers.push_back(split_tb_kmers);

        s.max_ubs_in_merged.push_back(max_ubs_in_merged);

        s.ibf_mem_size.push_back(max_cardinality * num_tbs);
        s.ibf_mem_size_no_corr.push_back(max_cardinality_no_corr * num_tbs);
    }
};
