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
#include <chopper/layout/phibf/partition_user_bins.hpp>
#include <chopper/next_multiple_of_64.hpp>

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
    bool const partitioned = config.number_of_partitions >= 2u; // 0 == unset == single HIBF, 1 == single HIBF
    assert(config.number_of_partitions >= 1u); // 0 == fail, 1 == no partitioning, => 2 PHIBF

    if (config.determine_best_tmax && config.fast_layout)
        throw std::invalid_argument{"determine_best_tmax is not supported with fast_layout."};

    if (partitioned && (config.determine_best_tmax))
        throw std::invalid_argument{"The partitioned HIBF (number_of_partitions >= 2) is not supported with "
                                    "determine_best_tmax."};

    config.hibf_config.validate_and_set_defaults();

    std::vector<size_t> cardinalities;
    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);
    std::vector<std::vector<size_t>> positions(config.number_of_partitions);
    std::vector<seqan::hibf::layout::layout> hibf_layouts(config.number_of_partitions);

    if (partitioned)
        phibf::partition_user_bins(config, cardinalities, sketches, minHash_sketches, positions);
    else
        positions[0] = seqan::hibf::iota_vector(sketches.size());

    if (config.determine_best_tmax)
    {
        hibf_layouts[0] = determine_best_number_of_technical_bins(config, cardinalities, sketches);
    }
    else
    {
#pragma omp parallel for schedule(dynamic) num_threads(config.hibf_config.threads)
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            // reset tmax to sqrt(n) where n is the number of user bins in the partition
            auto config_copy = config;

            if (!config.tmax_is_set)
            {
                config_copy.hibf_config.tmax =
                        chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(positions[i].size()))));
            }

            // use timers of config_copy for thread safe start() and stop() and add result to config afterwards
            config_copy.dp_algorithm_timer.start();
            if (config.fast_layout)
            {
                fast_layout(config_copy,
                            std::move(positions[i]),
                            cardinalities,
                            sketches,
                            minHash_sketches,
                            hibf_layouts[i]);
            }
            else
            {
                hibf_layouts[i] = seqan::hibf::layout::compute_layout(config_copy.hibf_config,
                                                                      cardinalities,
                                                                      sketches,
                                                                      std::move(positions[i]),
                                                                      config_copy.union_estimation_timer,
                                                                      config_copy.rearrangement_timer);
                config.union_estimation_timer += config_copy.union_estimation_timer;
                config.rearrangement_timer += config_copy.rearrangement_timer;
            }
            config_copy.dp_algorithm_timer.stop();

            config.dp_algorithm_timer += config_copy.dp_algorithm_timer;
        }
    }

    // this fails because each layout has only part of the user bins
    for (size_t i = 0; i < config.number_of_partitions; ++i)
        hibf_layouts[i].validate(config.hibf_config);

    // With determine_best_tmax, determine_best_number_of_technical_bins already wrote the statistics of every tmax to
    // `<output_filename>.stats`.
    if (!config.determine_best_tmax && config.output_verbose_statistics)
    {
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            size_t dummy{};
            chopper::layout::hibf_statistics global_stats{config, sketches, cardinalities};
            global_stats.hibf_layout = hibf_layouts[i];
            global_stats.print_header_to(std::cout);
            global_stats.print_summary_to(dummy, std::cout);
        }
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    chopper::layout::write_user_bins_to(filenames, fout);
    config.write_to(fout);

    for (size_t i = 0; i < config.number_of_partitions; ++i)
        hibf_layouts[i].write_to(fout);

    return 0;
}

} // namespace chopper::layout
