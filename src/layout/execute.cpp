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
#include <numeric>
#include <random>
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

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/toolbox.hpp>

namespace chopper::layout
{

void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions)
{
    // all approaches need sorted positions
    std::vector<size_t> const sorted_positions = [&cardinalities]()
    {
        std::vector<size_t> ps;
        ps.resize(cardinalities.size());
        std::iota(ps.begin(), ps.end(), 0);
        seqan::hibf::sketch::toolbox::sort_by_cardinalities(cardinalities, ps);
        return ps;
    }();

    if (config.partitioning_approach == partitioning_scheme::blocked)
    {
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);
        size_t const block_size =
            std::min(u_bins_per_part,
                     chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));

        size_t current_part{0u};
        size_t current_block_count{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            positions[current_part].push_back(current_user_bin_id);
            ++current_block_count;

            if (current_block_count >= block_size)
            {
                current_block_count = 0;
                ++current_part;
                if (current_part == config.number_of_partitions) // we need to circle back to the first partition
                    current_part = 0;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::sorted)
    {
        size_t const sum_of_cardinalities = std::accumulate(cardinalities.begin(), cardinalities.end(), size_t{});
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        size_t current_cardinality{0u};
        size_t current_part{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            positions[current_part].push_back(current_user_bin_id);
            current_cardinality += cardinalities[current_user_bin_id];

            if (current_cardinality >= cardinality_per_part)
            {
                current_cardinality = 0;
                ++current_part;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::folded)
    {
        size_t const sum_of_cardinalities = std::accumulate(cardinalities.begin(), cardinalities.end(), size_t{});
        size_t const cardinality_per_part_halved =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions * 2);

        size_t current_cardinality{0u};
        std::vector<size_t> const parts = [&config]()
        {
            size_t const len{config.number_of_partitions};
            std::vector<size_t> result(len * 2);
            std::iota(result.begin(), result.begin() + len, 0);
            std::copy(result.rbegin() + len, result.rend(), result.begin() + len);
            return result;
        }();
        size_t current_part{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            positions[parts[current_part]].push_back(current_user_bin_id);
            current_cardinality += cardinalities[current_user_bin_id];

            if (current_cardinality >= cardinality_per_part_halved)
            {
                current_cardinality = 0;
                ++current_part;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::weighted_fold)
    {
        size_t const sum_of_cardinalities = std::accumulate(cardinalities.begin(), cardinalities.end(), size_t{});
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);

        size_t current_big_pos{0};                          // the next largest user bin to assign to a partition
        size_t current_small_pos{cardinalities.size() - 1}; // the next small user bin

        for (size_t current_part = 0; current_part + 1 < config.number_of_partitions; ++current_part)
        {
            size_t current_cardinality{0};
            std::vector<size_t> small_bins;
            size_t new_small_bin_addition{0};

            auto compute_score = [&]()
            {
                double const weight = static_cast<double>(current_cardinality) / cardinality_per_part;
                double const amount =
                    static_cast<double>(positions[current_part].size() + small_bins.size() + new_small_bin_addition)
                    / u_bins_per_part;
                return std::abs(1.0 - weight) + std::abs(1.0 - amount);
            };

            // first add all large bins that fit
            while (current_cardinality < cardinality_per_part)
            {
                positions[current_part].push_back(sorted_positions[current_big_pos]);
                current_cardinality += cardinalities[sorted_positions[current_big_pos]];
                ++current_big_pos;
            }

            double local_optimum = compute_score();

            // then remove big bins and add small bins until a local optima is reached
            while (true)
            {
                size_t const cache_last_small_pos{current_small_pos};
                // remove a big user bin and fill the partition with small user bins
                current_cardinality -= cardinalities[sorted_positions[current_big_pos]];

                // can we further improve the ratio by adding more small bins?
                double improved_score{};
                do
                {
                    improved_score = compute_score();
                    current_cardinality += cardinalities[sorted_positions[current_small_pos]];
                    --current_small_pos;
                    ++new_small_bin_addition;
                }
                while (compute_score() < improved_score); // smaller is better
                // remove overstep
                ++current_small_pos;
                current_cardinality -= cardinalities[sorted_positions[current_small_pos]];
                --new_small_bin_addition;

                if (local_optimum < compute_score()) // score would increase. Stop
                {
                    current_small_pos = cache_last_small_pos;
                    break;
                }
                else // update
                {
                    positions[current_part].pop_back();
                    --current_big_pos;
                    for (size_t pos = cache_last_small_pos; pos > current_small_pos; --pos)
                        small_bins.push_back(sorted_positions[pos]);
                }
            }
            positions[current_part].insert(positions[current_part].end(), small_bins.begin(), small_bins.end());
        }

        // remaining user bins go to last partition
        while (current_big_pos <= current_small_pos)
        {
            positions[config.number_of_partitions - 1].push_back(sorted_positions[current_big_pos]);
            ++current_big_pos;
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::similarity)
    {
        uint8_t sketch_bits{10};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        size_t const sum_of_cardinalities = std::accumulate(cardinalities.begin(), cardinalities.end(), size_t{});
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);
        size_t const block_size =
            std::min(u_bins_per_part,
                     chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));
        size_t const number_of_blocks = seqan::hibf::divide_and_ceil(cardinalities.size(), block_size);

        // initialise partitions with the first config.number_of_partitions blocks
        assert(number_of_blocks >= config.number_of_partitions);
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            for (size_t x = 0; x < block_size; ++x)
            {
                partition_sketches[i].merge(sketches[sorted_positions[i + x]]);
                positions[i].push_back(sorted_positions[i + x]);
            }
        }

        // assign the rest by similarity
        // but don't move from largest to smallest but pick the next block to process randomly.
        // this probably leads to more evenly distributed partitions (evenly in terms of number of user bins)
        std::vector<size_t> indices(number_of_blocks - config.number_of_partitions);
        std::iota(indices.begin(), indices.end(), config.number_of_partitions);
        std::random_device shuffle_random_device;
        std::mt19937 shuffle_engine(shuffle_random_device());
        std::shuffle(indices.begin(), indices.end(), shuffle_engine);

        for (size_t const i : indices)
        {
            seqan::hibf::sketch::hyperloglog current_sketch(sketch_bits);

            // initialise sketch of the current block of indices
            for (size_t x = 0; x < block_size && (i + x) < sorted_positions.size(); ++x)
                current_sketch.merge(sketches[sorted_positions[i + x]]);

            // search best partition fit by similarity
            // similarity here is defined as:
            // "whose (<-partition) effective text size will increase the least when current_block is added to it"
            size_t smallest_change{std::numeric_limits<size_t>::max()};
            size_t best_p{0};
            for (size_t p = 0; p < config.number_of_partitions; ++p)
            {
                seqan::hibf::sketch::hyperloglog tmp = current_sketch;
                tmp.merge(partition_sketches[p]);

                if (tmp.estimate() < smallest_change && partition_sketches[p].estimate() < cardinality_per_part)
                {
                    smallest_change = tmp.estimate();
                    best_p = p;
                }
            }

            // now that we know which partition fits best (`best_p`), add those indices to it
            for (size_t x = 0; x < block_size && (i + x) < sorted_positions.size(); ++x)
                positions[best_p].push_back(i + x);
            partition_sketches[best_p].merge(current_sketch);
        }
    }
}

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

    if (config.number_of_partitions < 2) // 0 == unset == single HIBF, 1 == single HIBF
    {
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

            std::vector<size_t> positions = [&kmer_counts]()
            {
                std::vector<size_t> ps;
                ps.resize(kmer_counts.size());
                std::iota(ps.begin(), ps.end(), 0);
                return ps;
            }(); // GCOVR_EXCL_LINE

            dp_algorithm_timer.start();
            hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                              kmer_counts,
                                                              sketches,
                                                              std::move(positions),
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
    }
    else
    {
        std::vector<size_t> cardinalities;
        std::vector<seqan::hibf::sketch::hyperloglog> sketches;
        std::vector<std::vector<size_t>> positions(config.number_of_partitions); // asign positions for each partition

        // compute sketches of all user bins
        seqan::hibf::concurrent_timer compute_sketches_timer{};
        compute_sketches_timer.start();
        seqan::hibf::sketch::compute_sketches(config.hibf_config, cardinalities, sketches);
        compute_sketches_timer.stop();

        partition_user_bins(config, cardinalities, sketches, positions);

        std::vector<seqan::hibf::layout::layout> hibf_layouts(config.number_of_partitions); // multiple layouts

#pragma omp parallel for schedule(dynamic) num_threads(config.hibf_config.threads)
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            seqan::hibf::concurrent_timer union_estimation_timer{};
            seqan::hibf::concurrent_timer rearrangement_timer{};
            seqan::hibf::concurrent_timer dp_algorithm_timer{};

            // reset tmax to fit number of user bins in layout
            auto local_hibf_config = config.hibf_config; // every thread needs to set individual tmax
            local_hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(positions[i].size()))));

            dp_algorithm_timer.start();
            hibf_layouts[i] = seqan::hibf::layout::compute_layout(local_hibf_config,
                                                                  cardinalities,
                                                                  sketches,
                                                                  std::move(positions[i]),
                                                                  union_estimation_timer,
                                                                  rearrangement_timer);
            dp_algorithm_timer.stop();

            if (!config.output_timings.empty())
            {
                std::ofstream output_stream{config.output_timings, std::ios_base::app};
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
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);

        for (size_t i = 0; i < config.number_of_partitions; ++i)
            hibf_layouts[i].write_to(fout);
    }

    return 0;
}

} // namespace chopper::layout
