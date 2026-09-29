// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <chopper/layout/fast_layout.hpp>
#include <chopper/layout/partition_user_bins.hpp>

#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp>

namespace chopper::layout
{

/*!\brief Computes the layout of a subset of user bins with the regular (DP-based) HIBF layout algorithm.
 * \param[in] config        The configuration (uses `hibf_config`).
 * \param[in] positions     The global indices of the user bins to lay out.
 * \param[in] cardinalities The cardinality of each user bin, indexed by global user bin index.
 * \param[in] sketches      The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \returns A layout whose top level is the IBF of the given subset. `user_bins[i].idx` are global user bin indices;
 *          `previous_TB_indices` and `max_bins` are relative to this subset's top-level IBF.
 *
 * The union estimation and rearrangement timers are local and discarded.
 */
seqan::hibf::layout::layout general_layout(chopper::configuration const & config,
                                           std::vector<size_t> positions,
                                           std::vector<size_t> const & cardinalities,
                                           std::vector<seqan::hibf::sketch::hyperloglog> const & sketches)
{
    seqan::hibf::layout::layout hibf_layout;

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    dp_algorithm_timer.start();
    hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                      cardinalities,
                                                      sketches,
                                                      std::move(positions),
                                                      union_estimation_timer,
                                                      rearrangement_timer);
    dp_algorithm_timer.stop();

    return hibf_layout;
}

/*!\brief Decides whether the lower-level IBF for a merged bin is laid out by the fast layout or by general_layout.
 * \param[in] config        The configuration (uses `hibf_config.tmax`).
 * \param[in] positions     The global indices of the user bins in the merged bin.
 * \param[in] cardinalities The cardinality of each user bin, indexed by global user bin index.
 * \returns `true` if the fast layout should be used:
 *          - `false` if there are fewer than `64 * tmax` user bins. With few user bins per technical bin, the greedy
 *            fast layout merges only a handful of bins at a time, and the resulting heavy splitting on lower levels
 *            raises the FPR correction.
 *          - `true` if there are more than 500'000 user bins. The DP layout would take more than about half a day.
 *          - otherwise `true` only if no user bin is larger than `sum_of_cardinalities / tmax`, i.e., no user bin
 *            needs splitting and a merge-only layout suffices.
 */
bool do_I_need_a_fast_layout(chopper::configuration const & config,
                             std::vector<size_t> const & positions,
                             std::vector<size_t> const & cardinalities)
{
    // the fast layout heuristic would greedily merge even if merging only 2 bins at a time
    // merging only little number of bins is highly disadvantegous for lower levels because few bins
    // will be heavily split and this will raise the fpr correction for split bins
    // Thus, if the average number of user bins per technical bin is less then 64, we should not fast layout
    if (positions.size() < (64 * config.hibf_config.tmax))
        return false;

    if (positions.size() > 500'000) // layout takes more than half a day (should this be a user option?)
        return true;

    size_t largest_size{0};
    size_t sum_of_cardinalities{0};

    for (size_t const i : positions)
    {
        sum_of_cardinalities += cardinalities[i];
        largest_size = std::max(largest_size, cardinalities[i]);
    }

    size_t const cardinality_per_tb = sum_of_cardinalities / config.hibf_config.tmax;

    bool const largest_user_bin_might_be_split = largest_size > cardinality_per_tb;

    // if no splitting is needed, its worth it to use a fast-merge-only algorithm
    if (!largest_user_bin_might_be_split)
        return true;

    return false;
}

/*!\brief Records a lower-level IBF, computed by partition_user_bins, in the layout.
 * \param[in]     config      The configuration (uses the FPR settings of `hibf_config`).
 * \param[in,out] hibf_layout The global layout. Its `user_bins` must be indexed by global user bin index
 *                            (`user_bins[i].idx == i`), as initialised by fast_layout.
 * \param[in]     partitions  The technical bins of the new IBF: `partitions[t]` holds the global user bin indices
 *                            in technical bin `t`.
 * \param[in]     sketches    The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in]     previous    The merged-bin path identifying the new IBF. Every user bin in `partitions` must
 *                            currently have exactly this path as `previous_TB_indices`.
 *
 * - **Merged bin** (more than one user bin): `t` is appended to each user bin's `previous_TB_indices`. Their
 *   `storage_TB_id` is set later, when the next level down is recorded.
 * - **Split or single bin** (one user bin): `storage_TB_id` is set to `t`, and `number_of_technical_bins` to the
 *   number of consecutive technical bins holding the same user bin. partition_user_bins places split bins
 *   contiguously.
 * - **Empty technical bins** are skipped.
 *
 * The technical bin with the largest FPR-corrected size (relaxed correction for merged bins, split correction for
 * split bins) is appended to `hibf_layout.max_bins` as `(previous, max_bin_id)`.
 *
 * Not thread-safe; callers serialise it with `omp critical`.
 */
void add_level_to_layout(chopper::configuration const & config,
                         seqan::hibf::layout::layout & hibf_layout,
                         std::vector<std::vector<size_t>> const & partitions,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<size_t> const & previous)
{
    size_t max_bin_id{0};
    size_t max_size{0};

    auto const split_fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = config.hibf_config.maximum_fpr, //
                                                     .hash_count = config.hibf_config.number_of_hash_functions,
                                                     .t_max = partitions.size()});

    double const relaxed_fpr_correction = seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config.hibf_config.maximum_fpr, //
         .relaxed_fpr = config.hibf_config.relaxed_fpr,
         .hash_count = config.hibf_config.number_of_hash_functions});

    // we assume here that the user bins have been sorted by user bin id such that pos = idx
    for (size_t partition_idx{0}; partition_idx < partitions.size(); ++partition_idx)
    {
        auto const & partition = partitions[partition_idx];

        if (partition.size() > 1) // merged bin
        {
            seqan::hibf::sketch::hyperloglog current_sketch{sketches[0]}; // ensure same bit size
            current_sketch.reset();

            for (size_t const user_bin_id : partition)
            {
                assert(hibf_layout.user_bins[user_bin_id].idx == user_bin_id);
                auto & current_user_bin = hibf_layout.user_bins[user_bin_id];

                // update
                assert(previous == current_user_bin.previous_TB_indices);
                current_user_bin.previous_TB_indices.push_back(partition_idx);
                current_sketch.merge(sketches[user_bin_id]);
            }

            // update max_bin_id, max_size
            size_t const current_size = current_sketch.estimate() * relaxed_fpr_correction;
            if (current_size > max_size)
            {
                max_bin_id = partition_idx;
                max_size = current_size;
            }
        }
        else if (partition.size() == 0) // should not happen.. dge case?
        {
            continue;
        }
        else // single or split bin (partition.size() == 1)
        {
            auto & current_user_bin = hibf_layout.user_bins[partitions[partition_idx][0]];
            assert(current_user_bin.idx == partitions[partition_idx][0]);
            current_user_bin.storage_TB_id = partition_idx;
            current_user_bin.number_of_technical_bins = 1; // initialise to 1

            while (partition_idx + 1 < partitions.size() && partitions[partition_idx].size() == 1
                   && partitions[partition_idx + 1].size() == 1
                   && partitions[partition_idx][0] == partitions[partition_idx + 1][0])
            {
                ++current_user_bin.number_of_technical_bins;
                ++partition_idx;
            }

            // update max_bin_id, max_size
            size_t const current_size = sketches[current_user_bin.idx].estimate()
                                      * split_fpr_correction[current_user_bin.number_of_technical_bins];
            if (current_size > max_size)
            {
                max_bin_id = current_user_bin.storage_TB_id;
                max_size = current_size;
            }
        }
    }

    hibf_layout.max_bins.emplace_back(previous, max_bin_id); // add lower level meta information
}

/*!\brief Grafts a layout computed for a merged bin (see general_layout) into the global layout.
 * \param[in,out] child_layout The layout of the merged bin's subtree. Its `max_bins` are modified (prefixed with
 *                             `new_previous`) and copied.
 * \param[in,out] hibf_layout  The global layout. Its `user_bins` must be indexed by global user bin index.
 * \param[in]     new_previous The path of the merged bin that the child layout's top level is attached to.
 *
 * - The child's top-level IBF is added to `hibf_layout.max_bins` as `(new_previous, child.top_level_max_bin_id)`.
 * - Every other child max bin is added with `new_previous` prepended to its path.
 * - For every child user bin, the child's relative path is appended to the global user bin's `previous_TB_indices`
 *   (which already equals `new_previous`), and `storage_TB_id` and `number_of_technical_bins` are copied.
 *
 * Not thread-safe; callers serialise it with `omp critical`.
 */
void update_layout_from_child_layout(seqan::hibf::layout::layout & child_layout,
                                     seqan::hibf::layout::layout & hibf_layout,
                                     std::vector<size_t> const & new_previous)
{
    hibf_layout.max_bins.emplace_back(new_previous, child_layout.top_level_max_bin_id);

    for (auto & max_bin : child_layout.max_bins)
    {
        max_bin.previous_TB_indices.insert(max_bin.previous_TB_indices.begin(),
                                           new_previous.begin(),
                                           new_previous.end());
        hibf_layout.max_bins.push_back(max_bin);
    }

    for (auto const & user_bin : child_layout.user_bins)
    {
        auto & actual_user_bin = hibf_layout.user_bins[user_bin.idx];

        actual_user_bin.previous_TB_indices.insert(actual_user_bin.previous_TB_indices.end(),
                                                   user_bin.previous_TB_indices.begin(),
                                                   user_bin.previous_TB_indices.end());
        actual_user_bin.number_of_technical_bins = user_bin.number_of_technical_bins;
        actual_user_bin.storage_TB_id = user_bin.storage_TB_id;
    }
}

/*!\brief Lays out the lower-level IBF of one merged bin with the fast layout, and recursively lays out its children.
 * \param[in]     config           The configuration.
 * \param[in]     positions        The global indices of the user bins in the merged bin.
 * \param[in]     cardinalities    The cardinality of each user bin, indexed by global user bin index.
 * \param[in]     sketches         The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in]     minHash_sketches The MinHash tables of each user bin, indexed by global user bin index.
 * \param[in,out] hibf_layout      The global layout, updated in place.
 * \param[in]     previous         The merged-bin path of this IBF.
 *
 * Partitions `positions` into `tmax` technical bins (partition_user_bins) and records them (add_level_to_layout).
 * Each resulting merged bin is then handled like in fast_layout: another fast-layout recursion or a general_layout,
 * decided by do_I_need_a_fast_layout. The recursion runs sequentially within the calling OpenMP task. Only the
 * updates to `hibf_layout` are in critical sections.
 */
void fast_layout_recursion(chopper::configuration const & config,
                           std::vector<size_t> const & positions,
                           std::vector<size_t> const & cardinalities,
                           std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                           std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                           seqan::hibf::layout::layout & hibf_layout,
                           std::vector<size_t> const & previous)
{
    std::vector<std::vector<size_t>> tmax_partitions(config.hibf_config.tmax);

    // here we assume that we want to start with a fast layout
    partition_user_bins(config, positions, cardinalities, sketches, minHash_sketches, tmax_partitions);

#pragma omp critical
    {
        add_level_to_layout(config, hibf_layout, tmax_partitions, sketches, previous);
    }

    for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
    {
        auto const & partition = tmax_partitions[partition_idx];
        auto const new_previous = [&]()
        {
            auto cpy{previous};
            cpy.push_back(partition_idx);
            return cpy;
        }();

        if (partition.empty() || partition.size() == 1) // nothing to merge
            continue;

        if (do_I_need_a_fast_layout(config, partition, cardinalities))
        {
            fast_layout_recursion(config,
                                  partition,
                                  cardinalities,
                                  sketches,
                                  minHash_sketches,
                                  hibf_layout,
                                  new_previous); // recurse fast_layout
        }
        else
        {
            auto child_layout = general_layout(config, partition, cardinalities, sketches);

#pragma omp critical
            {
                update_layout_from_child_layout(child_layout, hibf_layout, new_previous);
            }
        }
    }
}

void fast_layout(chopper::configuration const & config,
                 std::vector<size_t> const & positions,
                 std::vector<size_t> const & cardinalities,
                 std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                 std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                 seqan::hibf::layout::layout & hibf_layout)
{
    std::vector<std::vector<size_t>> tmax_partitions(config.hibf_config.tmax);

    // here we assume that we want to start with a fast layout
    config.intital_partition_timer.start();
    partition_user_bins(config, positions, cardinalities, sketches, minHash_sketches, tmax_partitions);
    config.intital_partition_timer.stop();

    auto const split_fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = config.hibf_config.maximum_fpr, //
                                                     .hash_count = config.hibf_config.number_of_hash_functions,
                                                     .t_max = config.hibf_config.tmax});

    double const relaxed_fpr_correction = seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config.hibf_config.maximum_fpr, //
         .relaxed_fpr = config.hibf_config.relaxed_fpr,
         .hash_count = config.hibf_config.number_of_hash_functions});

    size_t max_bin_id{0};
    size_t max_size{0};
    hibf_layout.user_bins.resize(config.hibf_config.number_of_user_bins);

    // initialise user bins in layout
    for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
    {
        if (tmax_partitions[partition_idx].size() > 1) // merged bin
        {
            seqan::hibf::sketch::hyperloglog current_sketch{sketches[0]}; // ensure same bit size
            current_sketch.reset();

            for (size_t const user_bin_id : tmax_partitions[partition_idx])
            {
                hibf_layout.user_bins[user_bin_id] = {.previous_TB_indices = {partition_idx},
                                                      .storage_TB_id = 0 /*not determiend yet*/,
                                                      .number_of_technical_bins = 1 /*not determiend yet*/,
                                                      .idx = user_bin_id};
                current_sketch.merge(sketches[user_bin_id]);
            }

            // update max_bin_id, max_size
            size_t const current_size = current_sketch.estimate() * relaxed_fpr_correction;
            if (current_size > max_size)
            {
                max_bin_id = partition_idx;
                max_size = current_size;
            }
        }
        else if (tmax_partitions[partition_idx].size() == 0) // should not happen. Edge case?
        {
            continue;
        }
        else // single or split bin (tmax_partitions[partition_idx].size() == 1)
        {
            assert(tmax_partitions[partition_idx].size() == 1);
            size_t const user_bin_id = tmax_partitions[partition_idx][0];
            hibf_layout.user_bins[user_bin_id] = {.previous_TB_indices = {},
                                                  .storage_TB_id = partition_idx,
                                                  .number_of_technical_bins = 1 /*determiend below*/,
                                                  .idx = user_bin_id};

            while (partition_idx + 1 < tmax_partitions.size() && tmax_partitions[partition_idx].size() == 1
                   && tmax_partitions[partition_idx + 1].size() == 1
                   && tmax_partitions[partition_idx][0] == tmax_partitions[partition_idx + 1][0])
            {
                ++hibf_layout.user_bins[user_bin_id].number_of_technical_bins;
                ++partition_idx;
            }

            // update max_bin_id, max_size
            size_t const current_size =
                sketches[user_bin_id].estimate()
                * split_fpr_correction[hibf_layout.user_bins[user_bin_id].number_of_technical_bins];
            if (current_size > max_size)
            {
                max_bin_id = hibf_layout.user_bins[user_bin_id].storage_TB_id;
                max_size = current_size;
            }
        }
    }

    hibf_layout.top_level_max_bin_id = max_bin_id;

    config.small_layouts_timer.start();
#pragma omp parallel num_threads(config.hibf_config.threads)
#pragma omp single
    {
#pragma omp taskloop
        for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
        {
            auto const & partition = tmax_partitions[partition_idx];

            if (partition.empty() || partition.size() == 1) // nothing to merge
                continue;

            if (do_I_need_a_fast_layout(config, partition, cardinalities))
            {
                fast_layout_recursion(config,
                                      partition,
                                      cardinalities,
                                      sketches,
                                      minHash_sketches,
                                      hibf_layout,
                                      {partition_idx}); // recurse fast_layout
            }
            else
            {
                auto small_layout = general_layout(config, partition, cardinalities, sketches);

#pragma omp critical
                {
                    update_layout_from_child_layout(small_layout, hibf_layout, std::vector<size_t>{partition_idx});
                }
            }
        }
    }
    config.small_layouts_timer.stop();

    // sort records ascending by the number of bin indices (corresponds to the IBF levels)
    // GCOVR_EXCL_START
    std::ranges::sort(hibf_layout.max_bins,
                      std::ranges::less{},
                      [](auto const & mb)
                      {
                          // std::cref: compare the vector by reference instead of copying it
                          return std::make_tuple(mb.previous_TB_indices.size(), std::cref(mb.previous_TB_indices));
                      });
    // GCOVR_EXCL_STOP

#ifndef NDEBUG
    // sanity check in debug
    std::vector<size_t> layout_user_bins{};
    for (auto & user_bin : hibf_layout.user_bins)
        layout_user_bins.push_back(user_bin.idx);
    if (!std::ranges::is_permutation(layout_user_bins, positions))
        throw std::logic_error{"Not all/Wrong user bins have been assigned to the layout!"};
#endif
}

} // namespace chopper::layout