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
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/fast_layout.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/output.hpp>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>

namespace chopper::layout
{

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches)
{
    config.hibf_config.validate_and_set_defaults();

    std::vector<size_t> cardinalities;
    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);

    seqan::hibf::layout::layout hibf_layout;

    if (config.determine_best_tmax)
    {
        hibf_layout = determine_best_number_of_technical_bins(config, cardinalities, sketches);
    }
    else
    {
        config.dp_algorithm_timer.start();
        if (config.fast_layout)
        {
            fast_layout(config,
                        seqan::hibf::iota_vector(sketches.size()),
                        cardinalities,
                        sketches,
                        minHash_sketches,
                        hibf_layout);
            // sort records ascending by the number of bin indices (corresponds to the IBF levels)
            // GCOVR_EXCL_START
            std::ranges::sort(hibf_layout.max_bins,
                                [](auto const & r, auto const & l)
                                {
                                    if (r.previous_TB_indices.size() == l.previous_TB_indices.size())
                                        return std::ranges::lexicographical_compare(r.previous_TB_indices,
                                                                                    l.previous_TB_indices);
                                    else
                                        return r.previous_TB_indices.size() < l.previous_TB_indices.size();
                                });
            // GCOVR_EXCL_STOP
        }
        else
        {
            hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                              cardinalities,
                                                              sketches,
                                                              seqan::hibf::iota_vector(sketches.size()),
                                                              config.union_estimation_timer,
                                                              config.rearrangement_timer);
        }
        config.dp_algorithm_timer.stop();

        if (config.output_verbose_statistics)
        {
            size_t dummy{};
            chopper::layout::hibf_statistics global_stats{config, sketches, cardinalities};
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
