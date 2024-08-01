// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/next_multiple_of_64.hpp>
#include <chopper/sketch/output.hpp>

#include <hibf/misc/iota_vector.hpp>
#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

int execute(chopper::configuration & config, std::vector<std::vector<std::string>> const & filenames)
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

    seqan::hibf::concurrent_timer compute_sketches_timer{};
    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    if (config.determine_best_tmax)
    {
        std::tie(hibf_layout, sketches) = determine_best_number_of_technical_bins(config);
    }
    else
    {
        std::vector<size_t> kmer_counts;

        compute_sketches_timer.start();
        seqan::hibf::sketch::compute_sketches(config.hibf_config, kmer_counts, sketches);
        compute_sketches_timer.stop();
        dp_algorithm_timer.start();
        hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                          kmer_counts,
                                                          sketches,
                                                          seqan::hibf::iota_vector(sketches.size()),
                                                          union_estimation_timer,
                                                          rearrangement_timer);
        dp_algorithm_timer.stop();

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
            sketch::write_sketch_file(filenames[i][0], sketches[i], config);
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    chopper::layout::write_user_bins_to(filenames, fout);
    config.write_to(fout);
    hibf_layout.write_to(fout);

    if (!config.output_timings.empty())
    {
        std::ofstream output_stream{config.output_timings};
        output_stream << std::fixed << std::setprecision(2);
        output_stream << "sketching_in_seconds\t"
                      << "layouting_in_seconds\t"
                      << "union_estimation_in_seconds\t"
                      << "rearrangement_in_seconds\n";
        output_stream << compute_sketches_timer.in_seconds() << '\t';
        output_stream << dp_algorithm_timer.in_seconds() << '\t';
        output_stream << union_estimation_timer.in_seconds() << '\t';
        output_stream << rearrangement_timer.in_seconds() << '\t';
    }

    return 0;
}

} // namespace chopper::layout
