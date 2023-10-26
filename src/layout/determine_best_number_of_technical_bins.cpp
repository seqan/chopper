// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/next_multiple_of_64.hpp>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

std::pair<seqan::hibf::layout::layout, std::vector<seqan::hibf::sketch::hyperloglog>>
determine_best_number_of_technical_bins(chopper::configuration & config)
{
    seqan::hibf::layout::layout best_layout;

    std::set<size_t> potential_t_max = [&]()
    {
        std::set<size_t> result;

        for (size_t t_max = 64; t_max <= config.hibf_config.tmax; t_max *= 2)
            result.insert(t_max);

        // Additionally, add the t_max that is closest to the sqrt() of the number of
        // user bins, as it is expected to evenly spread bins and may perform well.
        size_t const user_bin_count{config.hibf_config.number_of_user_bins};
        size_t const sqrt_t_max{next_multiple_of_64(std::ceil(std::sqrt(user_bin_count)))};
        result.insert(sqrt_t_max);

        return result;
    }();

    // with -determine-best-tmax the algorithm is executed multiple times and result with the minimum
    // expected query costs are written to the standard output

    std::ofstream file_out{config.output_filename.string() + ".stats"};

    file_out << "## ### Parameters ###\n"
             << "## number of user bins = " << config.hibf_config.number_of_user_bins << '\n'
             << "## number of hash functions = " << config.hibf_config.number_of_hash_functions << '\n'
             << "## maximum false positive rate = " << config.hibf_config.maximum_fpr << '\n'
             << "## relaxed false positive rate = " << config.hibf_config.relaxed_fpr << '\n';
    hibf_statistics::print_header_to(file_out, config.output_verbose_statistics);

    double best_expected_HIBF_query_cost{std::numeric_limits<double>::infinity()};
    size_t best_t_max{};
    size_t t_max_64_memory{};

    std::vector<size_t> kmer_counts;
    std::vector<seqan::hibf::sketch::hyperloglog> sketches;

    for (size_t const t_max : potential_t_max)
    {
        config.hibf_config.tmax = t_max;

        kmer_counts.clear();
        sketches.clear();

        seqan::hibf::sketch::compute_sketches(config.hibf_config, kmer_counts, sketches);
        seqan::hibf::layout::layout tmp_layout =
            seqan::hibf::layout::compute_layout(config.hibf_config, kmer_counts, sketches);

        chopper::layout::hibf_statistics global_stats{config, sketches, kmer_counts};
        global_stats.hibf_layout = tmp_layout;
        global_stats.finalize();
        global_stats.print_summary_to(t_max_64_memory, file_out, config.output_verbose_statistics);

        // Use result if better than previous one.
        if (global_stats.expected_HIBF_query_cost < best_expected_HIBF_query_cost)
        {
            best_layout = std::move(tmp_layout);
            best_t_max = t_max;
            best_expected_HIBF_query_cost = global_stats.expected_HIBF_query_cost;
        }
        else if (!config.force_all_binnings)
        {
            break;
        }
    }

    file_out << "# Best t_max (regarding expected query runtime): " << best_t_max << '\n';
    config.hibf_config.tmax = best_t_max;

    return {best_layout, sketches};
}

} // namespace chopper::layout
