// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

// clang-format off
#include <chopper/workarounds.hpp>
// clang-format on

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wrestrict"
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#include <iostream>
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/ibf_query_cost.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

hibf_statistics::hibf_statistics(configuration const & config_,
                                 std::vector<seqan::hibf::sketch::hyperloglog> const & sketches_,
                                 std::vector<size_t> const & kmer_counts) :
    config{config_},
    fp_correction{
        seqan::hibf::layout::compute_fpr_correction({.fpr = config_.hibf_config.maximum_fpr,
                                                     .hash_count = config_.hibf_config.number_of_hash_functions,
                                                     .t_max = config_.hibf_config.tmax})},
    merged_fpr_correction_factor{seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config_.hibf_config.maximum_fpr,
         .relaxed_fpr = config_.hibf_config.relaxed_fpr,
         .hash_count = config_.hibf_config.number_of_hash_functions})},
    sketches{sketches_},
    counts{kmer_counts},
    total_kmer_count{std::accumulate(kmer_counts.begin(), kmer_counts.end(), size_t{})}
{}

void hibf_statistics::finalize()
{
    collect_bins();

    compute_cardinalities(top_level_ibf);

    compute_total_query_cost(top_level_ibf);

    gather_statistics(top_level_ibf, 0);

    expected_HIBF_query_cost = total_query_cost / total_kmer_count;
}

//!\brief Prints a column names of the summary to the command line.
void hibf_statistics::print_header_to(std::ostream & stream, bool const verbose)
{
    // print column names explanation in header
    stream << "## ### Notation ###\n"
           << "## X-IBF = An IBF with X number of bins.\n"
           << "## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.\n";

    stream << "## ### Column Description ###\n"
              "## tmax : The maximum number of technical bin on each level\n"
              "## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF\n"
              "## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF\n"
              "## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF\n"
              "## (l*m)_tmax : Computed by l_tmax * m_tmax\n"
              "## size : The expected total size of an tmax-HIBF\n"
           << ((verbose) ? "## uncorr_size : The expected size of an tmax-HIBF without FPR correction\n" : "");

    // print column names
    stream << "# tmax" << '\t' << "c_tmax" << '\t' << "l_tmax" << '\t' << "m_tmax" << '\t' << "(l*m)_tmax" << '\t'
           << "size";

    if (verbose) // uncorrected size and add level statistics
    {
        stream << '\t' << "uncorr_size" << '\t' << "level" << '\t' << "num_ibfs" << '\t' << "level_size" << '\t'
               << "level_size_no_corr" << '\t' << "total_num_tbs" << '\t' << "avg_num_tbs" << '\t'
               << "split_tb_percentage" << '\t' << "max_split_tb" << '\t' << "avg_split_tb" << '\t' << "max_factor"
               << '\t' << "avg_factor";
    }

    stream << '\n';
}

void hibf_statistics::print_summary_to(size_t & t_max_64_memory, std::ostream & stream, bool const verbose)
{
    if (summaries.empty())
        finalize();

    if (t_max_64_memory == 0)
        t_max_64_memory = total_hibf_size_in_byte();

    double const relative_memory_size = total_hibf_size_in_byte() / static_cast<double>(t_max_64_memory);
    double const query_time_memory_usage_prod = expected_HIBF_query_cost * relative_memory_size;

    stream << std::fixed << std::setprecision(2);

    std::string level_str, num_ibfs_str, level_size_str, level_size_no_corr_str, total_num_tbs_str, avg_num_tbs_str,
        split_tb_percentage_str, max_split_tb_str, avg_split_tb_str, max_factor_str, avg_factor_str;

    size_t total_size{};
    size_t total_size_no_corr{};

    // go through each level and collect and output the statistics
    auto to_string_with_precision = [](auto num)
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

#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wrestrict"
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY

        level_str += ":" + to_string_with_precision(level);
        num_ibfs_str += ":" + to_string_with_precision(s.num_ibfs);
        level_size_str += ":" + to_formatted_BF_size(level_size);
        level_size_no_corr_str += ":" + to_formatted_BF_size(level_size_no_corr);
        total_num_tbs_str += ":" + to_string_with_precision(total_num_tbs);
        avg_num_tbs_str += ":" + to_string_with_precision(total_num_tbs / s.num_ibfs);
        split_tb_percentage_str += ":" + to_string_with_precision(split_tb_percentage);

        // if there are no split bins on this level, the following statistics don't make sense
        if (max_split_bin_span != 0)
        {
            size_t const total_num_split_ubs = std::reduce(s.num_split_ubs.begin(), s.num_split_ubs.end());
            double const avg_split_bin =
                static_cast<double>(total_num_split_tbs) / static_cast<double>(total_num_split_ubs);
            size_t const total_split_tb_kmers = std::reduce(s.split_tb_kmers.begin(), s.split_tb_kmers.end());
            double const avg_factor =
                static_cast<double>(std::reduce(s.split_tb_corr_kmers.begin(), s.split_tb_corr_kmers.end()))
                / static_cast<double>(total_split_tb_kmers);

            max_split_tb_str += ":" + to_string_with_precision(max_split_bin_span);
            avg_split_tb_str += ":" + to_string_with_precision(avg_split_bin);
            max_factor_str += ":" + to_string_with_precision((fp_correction)[max_split_bin_span]);
            avg_factor_str += ":" + to_string_with_precision(avg_factor);
        }
        else
        {
            max_split_tb_str += ":-";
            avg_split_tb_str += ":-";
            max_factor_str += ":-";
            avg_factor_str += ":-";
        }
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMCPY
    }

    stream << std::fixed << std::setprecision(2);

    stream /*        tmax */ << config.hibf_config.tmax
                             << '\t'
                             /*      c_tmax */
                             << chopper::layout::ibf_query_cost::interpolated(config.hibf_config.tmax,
                                                                              config.hibf_config.maximum_fpr)
                             << '\t'
                             /*      l_tmax */
                             << expected_HIBF_query_cost
                             << '\t' /*relative to a 64 bin IBF*/
                                     /*      m_tmax */
                             << relative_memory_size
                             << '\t' /*relative to the 64 T_Max HIBF*/
                                     /*   (l*m)tmax */
                             << query_time_memory_usage_prod
                             << '\t'
                             /*  corr. size */
                             << to_formatted_BF_size(total_size) << ((verbose) ? '\t' : '\n');

    if (verbose)
    {
        // uncorrected FPR
        stream /*uncorr. size */ << to_formatted_BF_size(total_size_no_corr) << '\t';

        // per level statistics:
        stream /* level               */ << level_str
                                         << '\t'
                                         /* num_ibfs            */
                                         << num_ibfs_str
                                         << '\t'
                                         /* level_size          */
                                         << level_size_str
                                         << '\t'
                                         /* level_size_no_corr  */
                                         << level_size_no_corr_str
                                         << '\t'
                                         /* total_num_tbs       */
                                         << total_num_tbs_str
                                         << '\t'
                                         /* avg_num_tbs         */
                                         << avg_num_tbs_str
                                         << '\t'
                                         /* split_tb_percentage */
                                         << split_tb_percentage_str
                                         << '\t'
                                         /* max_split_tb        */
                                         << max_split_tb_str
                                         << '\t'
                                         /* avg_split_tb        */
                                         << avg_split_tb_str
                                         << '\t'
                                         /* max_factor          */
                                         << max_factor_str
                                         << '\t'
                                         /* avg_factor          */
                                         << avg_factor_str << '\n';
    }
}

//!\brief Return the total corrected size of the HIBF in bytes
size_t hibf_statistics::total_hibf_size_in_byte()
{
    if (summaries.empty())
        finalize();

    size_t total_size{};

    // go through each level and collect the memory sizes
    for (auto const & [level, summary] : summaries)
    {
        (void)level;

        total_size += std::reduce(summary.ibf_mem_size.begin(), summary.ibf_mem_size.end());
    }

    size_t const size_in_bits =
        seqan::hibf::build::bin_size_in_bits({.fpr = config.hibf_config.maximum_fpr,
                                              .hash_count = config.hibf_config.number_of_hash_functions,
                                              .elements = total_size});

    return size_in_bits / 8;
}

//!\brief Round bytes to the appropriate unit and convert to string with unit.
[[nodiscard]] std::string hibf_statistics::byte_size_to_formatted_str(size_t const bytes)
{
    size_t iterations{};
    size_t integer{bytes};

    while (integer >> 10u && iterations < 6u)
    {
        integer >>= 10u;
        ++iterations;
    }

    // While this is a bit more involved, we can avoid using floating point numbers.
    auto first_decimal_position = [&]()
    {
        assert(iterations > 0u);
        size_t decimal{bytes};
        decimal -= integer << (iterations * 10u);             // Substract bytes represented by integer, e.g. -5GiB
        decimal >>= (iterations - 1u) * 10u;                  // Shift to next smallest unit, e.g. 800MiB
        decimal = decimal * 1000u / 1024u;                    // Account for using decimal system, i.e. 800MiB != 0.8GiB
        size_t const diff{decimal - (decimal / 100u) * 100u}; // We want to round up to 1 decimal position
        uint32_t const round_up{diff >= 50u};
        decimal += round_up * 100u - diff;
        decimal /= 100u;
        return decimal;
    };

    auto formatted_string = [&]()
    {
        static constexpr int8_t int_to_char_offset{'0'}; // int 0 as char: char{0 + 48} = '0'
        size_t const decimal = iterations ? first_decimal_position() : 0u;
        assert(decimal <= 10u);

        if (!iterations) // No decimals for Bytes
            return std::to_string(integer);
        else if (decimal < 10u) // No need to round integer part
            return std::to_string(integer) + '.' + static_cast<char>(decimal + int_to_char_offset);
        else // Round integer part, e.g., 5.99 MiB should report 6.0 MiB
        {
            ++integer;
            // Check whether rounding results in a change of unit, e.g. 1023.99MiB to 1.0GiB
            if (integer >> 10u)
            {
                ++iterations;
                integer >>= 10u;
            }
            return std::to_string(integer) + ".0";
        }
    };

    std::string result{formatted_string()};
    switch (iterations)
    {
    case 0:
        result += "Bytes";
        break;
    case 1:
        result += "KiB";
        break;
    case 2:
        result += "MiB";
        break;
    case 3:
        result += "GiB";
        break;
    case 4:
        result += "TiB";
        break;
    case 5:
        result += "PiB";
        break;
    default:
        result += "EiB";
        break;
    }

    return result;
}

std::string hibf_statistics::to_formatted_BF_size(size_t const number_of_kmers_to_be_stored) const
{
    size_t const size_in_bits =
        seqan::hibf::build::bin_size_in_bits({.fpr = config.hibf_config.maximum_fpr,
                                              .hash_count = config.hibf_config.number_of_hash_functions,
                                              .elements = number_of_kmers_to_be_stored});
    return byte_size_to_formatted_str(size_in_bits / 8);
}

void hibf_statistics::collect_bins()
{
    std::vector<hibf_statistics::level> ibfs(hibf_layout.max_bins.size() + 1); // 0 = top_level
    robin_hood::unordered_map<std::vector<size_t>, size_t> id_to_pos{};

    // fill id_to_pos map
    id_to_pos[std::vector<size_t>{}] = 0;
    for (size_t i = 0; i < hibf_layout.max_bins.size(); ++i)
        id_to_pos[hibf_layout.max_bins[i].previous_TB_indices] = i + 1;

    for (auto const & user_bin_info : hibf_layout.user_bins)
    {
        std::vector<size_t> prev{};

        // add user bin index to previous merged bins
        for (size_t i = 0; i < user_bin_info.previous_TB_indices.size(); ++i)
        {
            auto & ibf = ibfs[id_to_pos.at(prev)];
            auto const target_tb_index = user_bin_info.previous_TB_indices[i];

            bool found_merged_bin{false};
            for (auto & previous_bins_to_check : ibf.bins)
            {
                if (previous_bins_to_check.tb_index == target_tb_index)
                {
                    found_merged_bin = true;
                    previous_bins_to_check.user_bin_indices.push_back(user_bin_info.idx);
                    ++previous_bins_to_check.num_contained_ubs;
                }
            }

            if (!found_merged_bin)
            {
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Warray-bounds="
#    pragma GCC diagnostic ignored "-Wstringop-overflow="
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
                // ibf.bins.emplace_back(hibf_statistics::bin_kind::merged, 1, std::vector<size_t>{user_bin_info.idx});
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
#    pragma GCC diagnostic pop
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
                ibf.bins.back().tb_index = target_tb_index;
                auto next = prev;
                next.push_back(target_tb_index);
                ibf.bins.back().child_level_idx = id_to_pos.at(next);
            }
            prev.push_back(target_tb_index);
        }

        // emplace a split bin at last since every user bin is on its lowest level single or split
        auto & ibf = ibfs[id_to_pos.at(prev)];
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Warray-bounds="
#    pragma GCC diagnostic ignored "-Wstringop-overflow="
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
        ibf.bins.emplace_back(hibf_statistics::bin_kind::split,
                              user_bin_info.number_of_technical_bins,
                              std::vector<size_t>{user_bin_info.idx});
#if CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
#    pragma GCC diagnostic pop
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_MEMMOV
        ibf.bins.back().tb_index = user_bin_info.storage_TB_id;
    }

    for (auto & ibf : ibfs)
        for (auto & bin : ibf.bins)
            if (bin.kind == hibf_statistics::bin_kind::merged)
                bin.child_level = ibfs[bin.child_level_idx];

    top_level_ibf = std::move(ibfs[0]);
}

void hibf_statistics::compute_cardinalities(level & curr_level)
{
    for (bin & current_bin : curr_level.bins)
    {
        if (current_bin.kind == bin_kind::merged)
        {
            if (config.hibf_config.disable_estimate_union)
            {
                size_t sum{};
                for (size_t i = 0; i < current_bin.user_bin_indices.size(); ++i)
                    sum += counts[current_bin.user_bin_indices[i]]; // TODO should be kmer_counts
                current_bin.cardinality = sum;
            }
            else
            {
                assert(!current_bin.user_bin_indices.empty());
                seqan::hibf::sketch::hyperloglog hll = sketches[current_bin.user_bin_indices[0]];

                for (size_t i = 1; i < current_bin.user_bin_indices.size(); ++i)
                    hll.merge(sketches[current_bin.user_bin_indices[i]]);

                current_bin.cardinality = hll.estimate();
            }

            compute_cardinalities(current_bin.child_level);
        }
        else if (current_bin.kind == bin_kind::split) // bin_kind::split
        {
            assert(current_bin.user_bin_indices.size() == 1);
            current_bin.cardinality = counts[current_bin.user_bin_indices[0]];
        }
    }
}

void hibf_statistics::compute_total_query_cost(level & curr_level)
{
    // Compute number of technical bins in current level (<= tmax)
    size_t number_of_tbs{0};
    size_t level_kmer_count{0};
    size_t index{0};
    std::vector<size_t> merged_bin_indices{};
    std::vector<seqan::hibf::sketch::hyperloglog> merged_bin_sketches{};

    for (bin const & current_bin : curr_level.bins)
    {
        if (current_bin.kind == bin_kind::merged)
        {
            ++number_of_tbs;
            merged_bin_indices.push_back(index);

            if (!config.hibf_config.disable_estimate_union)
            {
                // compute merged_bin_sketch
                assert(!current_bin.user_bin_indices.empty());
                seqan::hibf::sketch::hyperloglog hll = sketches[current_bin.user_bin_indices[0]];

                for (size_t i = 1; i < current_bin.user_bin_indices.size(); ++i)
                    hll.merge(sketches[current_bin.user_bin_indices[i]]);

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
    assert(number_of_tbs <= config.hibf_config.tmax);

    // Add cost of querying the current IBF
    // (how costly is querying number_of_tbs (e.g. 128 tbs) compared to 64 tbs given the current FPR)
    curr_level.current_query_cost += ibf_query_cost::interpolated(number_of_tbs, config.hibf_config.maximum_fpr);

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
        if (!config.hibf_config.disable_estimate_union)
        {
            double const current_estimate = merged_bin_sketches[i].estimate();

            for (size_t j = i + 1; j < merged_bin_indices.size(); ++j)
            {
                seqan::hibf::sketch::hyperloglog tmp =
                    merged_bin_sketches[i]; // copy needed, s.t. current is not modified
                double union_estimate = tmp.merge_and_estimate(merged_bin_sketches[j]);
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

void hibf_statistics::gather_statistics(level const & curr_level, size_t const level_summary_index)
{
    level_summary & summary = summaries[level_summary_index];
    summary.num_ibfs += 1;

    size_t max_cardinality{}, max_cardinality_no_corr{}, num_tbs{}, num_ubs{}, num_split_tbs{}, num_merged_tbs{},
        num_split_ubs{}, num_merged_ubs{}, max_split_tb_span{}, split_tb_kmers{}, max_ubs_in_merged{},
        split_tb_corr_kmers{};

    for (bin const & current_bin : curr_level.bins)
    {
        size_t uncorrected_cardinality{};
        size_t corrected_cardinality{};

        if (current_bin.kind == bin_kind::split)
        {
            uncorrected_cardinality =
                (current_bin.cardinality + current_bin.num_spanning_tbs - 1) / current_bin.num_spanning_tbs; // round up
            corrected_cardinality = std::ceil(uncorrected_cardinality * (fp_correction)[current_bin.num_spanning_tbs]);
        }
        else // current_bin.kind == bin_kind::merged
        {
            uncorrected_cardinality = current_bin.cardinality;
            corrected_cardinality = std::ceil(uncorrected_cardinality * merged_fpr_correction_factor);
        }

        max_cardinality = std::max(max_cardinality, corrected_cardinality);
        max_cardinality_no_corr = std::max(max_cardinality_no_corr, uncorrected_cardinality);

        num_tbs += current_bin.num_spanning_tbs;
        num_ubs += current_bin.num_contained_ubs;

        if (current_bin.kind == bin_kind::split)
        {
            num_split_tbs += current_bin.num_spanning_tbs;
            num_split_ubs += 1;
            split_tb_corr_kmers += corrected_cardinality * current_bin.num_spanning_tbs;
            split_tb_kmers += uncorrected_cardinality * current_bin.num_spanning_tbs;
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

} // namespace chopper::layout
