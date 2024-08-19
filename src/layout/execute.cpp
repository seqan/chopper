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

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches)
{
    config.hibf_config.validate_and_set_defaults();

    seqan::hibf::layout::layout hibf_layout;
    std::vector<size_t> kmer_counts;
    seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);

    if (config.determine_best_tmax)
    {
        hibf_layout = determine_best_number_of_technical_bins(config, kmer_counts, sketches);
    }
    else
    {
        config.dp_algorithm_timer.start();
        hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                          kmer_counts,
                                                          sketches,
                                                          seqan::hibf::iota_vector(sketches.size()),
                                                          config.union_estimation_timer,
                                                          config.rearrangement_timer);
        config.dp_algorithm_timer.stop();

        if (config.output_verbose_statistics)
        {
            size_t dummy{};
            chopper::layout::hibf_statistics global_stats{config, sketches, kmer_counts};
            global_stats.hibf_layout = hibf_layout;
            global_stats.print_header_to(std::cout);
            global_stats.print_summary_to(dummy, std::cout);
        }
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    chopper::layout::write_user_bins_to(filenames, fout);
    config.write_to(fout);
    hibf_layout.write_to(fout);

    return 0;
}

} // namespace chopper::layout
