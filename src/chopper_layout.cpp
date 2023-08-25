// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <iostream>
#include <set>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/ibf_query_cost.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/sketch/output.hpp>

#include <hibf/detail/layout/compute_layout.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

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
             << "## false positive rate = " << config.hibf_config.maximum_false_positive_rate << '\n';
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

int execute(chopper::configuration & config, std::vector<std::string> const & filenames)
{
    assert(config.hibf_config.number_of_user_bins > 0);

    if (config.hibf_config.disable_estimate_union)
        config.hibf_config.disable_rearrangement = true;

    if (config.hibf_config.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (size_t number_samples = config.hibf_config.number_of_user_bins;
            number_samples >= 1ULL << 32) // sqrt is bigger than uint16_t
            throw std::invalid_argument{"Too many samples. Please set a tmax (see help via `-hh`)."}; // GCOVR_EXCL_LINE
        else
            config.hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_samples))));
    }
    else if (config.hibf_config.tmax % 64 != 0)
    {
        config.hibf_config.tmax = chopper::next_multiple_of_64(config.hibf_config.tmax);
        std::cerr << "[CHOPPER LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.hibf_config.tmax << ".\n";
    }

    seqan::hibf::layout::layout hibf_layout;
    std::vector<seqan::hibf::sketch::hyperloglog> sketches;

    if (config.determine_best_tmax)
    {
        std::tie(hibf_layout, sketches) = determine_best_number_of_technical_bins(config);
    }
    else
    {
        std::vector<size_t> kmer_counts;

        hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config, kmer_counts, sketches);

        if (config.output_verbose_statistics)
        {
            size_t dummy{};
            chopper::layout::hibf_statistics global_stats{config, sketches, kmer_counts};
            global_stats.hibf_layout = hibf_layout;
            global_stats.print_header_to(std::cout);
            global_stats.print_summary_to(dummy, std::cout);
        }
    }

    if (!config.disable_sketch_output)
    {
        if (!std::filesystem::exists(config.sketch_directory))
            std::filesystem::create_directory(config.sketch_directory);

        assert(filenames.size() == sketches.size());
        for (size_t i = 0; i < filenames.size(); ++i)
            sketch::write_sketch_file(filenames[i], sketches[i], config);
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    chopper::layout::write_user_bins_to(filenames, fout);
    chopper::layout::write_config_to(config, fout);
    hibf_layout.write_to(fout);

    return 0;
}

} // namespace chopper::layout
