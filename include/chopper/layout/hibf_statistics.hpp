// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <numeric>
#include <string>
#include <typeindex>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace std
{

template <>
struct hash<std::vector<size_t>>
{
    size_t operator()(std::vector<size_t> const & vec) const
    {
        return std::accumulate(vec.begin(),
                               vec.end(),
                               vec.size(),
                               [](std::size_t hash, size_t value)
                               {
                                   return hash ^= value + 0x9e3779b9ULL + (hash << 6) + (hash >> 2);
                               });
    }
};

} // namespace std

namespace chopper::layout
{

class hibf_statistics
{
public:
    hibf_statistics() = delete;                                    //!< Deleted. Holds reference members.
    hibf_statistics(hibf_statistics const & b) = delete;           //!< Deleted. Holds const member.
    hibf_statistics & operator=(hibf_statistics const &) = delete; //!< Deleted. Holds const member.
    hibf_statistics(hibf_statistics && b) = delete;                //!< Deleted. Holds const member.
    hibf_statistics & operator=(hibf_statistics &&) = delete;      //!< Deleted. Holds const member.
    ~hibf_statistics() = default;                                  //!< Defaulted.

    /*!\brief Construct an empty HIBF with an empty top level IBF
     * \param[in] config_ User configuration for the HIBF.
     * \param[in] sketches_ The sketches of the input.
     * \param[in] kmer_counts The original user bin weights (kmer counts).
     */
    hibf_statistics(configuration const & config_,
                    std::vector<seqan::hibf::sketch::hyperloglog> const & sketches_,
                    std::vector<size_t> const & kmer_counts);

    //!\brief Represents a (set) of user bins (see ibf_statistics::bin_kind).
    class bin;

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

    //!\brief Gather all statistics to have all members ready.
    void finalize();

    //!\brief Prints a column names of the summary to the command line.
    static void print_header_to(std::ostream & stream, bool const verbose = true);

    //!\brief Prints a tab-separated summary of the statistics of this HIBF to the command line.
    void print_summary_to(size_t & t_max_64_memory, std::ostream & stream, bool const verbose = true);

    //!\brief Return the total corrected size of the HIBF in bytes
    size_t total_hibf_size_in_byte();

    //!\brief Round bytes to the appropriate unit and convert to string with unit.
    [[nodiscard]] static std::string byte_size_to_formatted_str(size_t const bytes);

    //!\brief The top level IBF of this HIBF, often starting point for recursions.
    level top_level_ibf;

    //!\brief The estimated query cost of every single kmer in this HIBF.
    double total_query_cost{0.0};

    //!\brief The estimated query cost relative to the total k-mer count in the data set.
    double expected_HIBF_query_cost{0.0};

    //!\brief A reference to the input counts.
    seqan::hibf::layout::layout hibf_layout;

private:
    //!\brief Copy of the user configuration for this HIBF.
    configuration const config{};

    //!\brief The false positive correction factors to use for the statistics.
    std::vector<double> const fp_correction{};

    //!\brief A reference to the input sketches.
    std::vector<seqan::hibf::sketch::hyperloglog> const & sketches;

    //!\brief A reference to the input counts.
    std::vector<size_t> const & counts;

    //!\brief The original kmer count of all user bins.
    size_t const total_kmer_count{};

    //!\brief Statistics for all IBFs on a certain level of the HIBF.
    struct level_summary;

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
    size_t compute_bin_size(size_t const number_of_kmers_to_be_stored) const;

    /*!\brief Compute the Bloom Filter size from `number_of_kmers_to_be_stored` and
     *        return it as a formatted string with the appropriate unit.
     * \param[in] number_of_kmers_to_be_stored
     */
    std::string to_formatted_BF_size(size_t const number_of_kmers_to_be_stored) const;

    void collect_bins();

    void compute_cardinalities(level & curr_level);

    //!\brief Computes the estimated query cost
    void compute_total_query_cost(level & curr_level);

    /*!\brief Recursively gather all the statistics from the bins.
     * \param[in] curr_level The current IBF from which the statistics will be extracted.
     * \param[in] level_summary_index The index of `curr_level` in `summeries`.
     */
    void gather_statistics(level const & curr_level, size_t const level_summary_index);
};

class hibf_statistics::bin
{
public:
    bin_kind kind;            //!< Either a split or merged bin.
    size_t cardinality;       //!< The size/weight of the bin (either a kmer count or hll sketch estimation).
    size_t num_contained_ubs; //!< [MERGED] How many UBs are merged within this TB.
    size_t num_spanning_tbs;  //!< [SPLIT] How many TBs are used for this sindle UB.
    std::vector<size_t> user_bin_indices; //!< The user bin indices of this bin.
    size_t tb_index;                      // The (first) technical bin idx this bin is stored in.
    level child_level;                    //!< [MERGED] The lower level ibf statistics.
    size_t child_level_idx;               //!< [MERGED] The lower level ibf statistics.

    bin() = default;                        //!< Defaulted.
    bin(bin const & b) = default;           //!< Defaulted.
    bin & operator=(bin const &) = default; //!< Defaulted.
    bin(bin && b) = default;                //!< Defaulted.
    bin & operator=(bin &&) = default;      //!< Defaulted.
    ~bin() = default;                       //!< Defaulted.

    bin(bin_kind const kind_, size_t const spanning_tbs, std::vector<size_t> const & user_bin_indices_) :
        kind{kind_},
        num_contained_ubs{user_bin_indices_.size()},
        num_spanning_tbs{spanning_tbs},
        user_bin_indices{user_bin_indices_}
    {
        assert((kind == bin_kind::split && num_contained_ubs == 1u)
               || (kind == bin_kind::merged && num_spanning_tbs == 1u));
    }
};

struct hibf_statistics::level_summary
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

} // namespace chopper::layout
