// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout
{

/*!\brief Computes an HIBF layout with the fast (LSH/similarity-based) layout algorithm.
 * \param[in]  config           The configuration. Uses `hibf_config` (`tmax`, `number_of_user_bins`, FPR settings)
 *                              and the fast-layout timers.
 * \param[in]  positions        The global indices of all user bins. Must be a permutation of
 *                              `[0, config.hibf_config.number_of_user_bins)`.
 * \param[in]  cardinalities    The cardinality of each user bin, indexed by global user bin index.
 * \param[in]  sketches         The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in]  minHash_sketches The MinHash tables of each user bin, indexed by global user bin index.
 * \param[out] hibf_layout      The resulting layout. Expected to be empty on entry. `user_bins` is resized to
 *                              `number_of_user_bins`, so `user_bins[i].idx == i`.
 *
 * 1. **Top level:** lsh_distributed_ibf_layout distributes all user bins onto `tmax` technical bins. The user bins are
 *    initialised in `hibf_layout`: merged bins get `previous_TB_indices = {t}`; split and single bins get
 *    `storage_TB_id` and their number of consecutive technical bins. `top_level_max_bin_id` is set to the technical
 *    bin with the largest FPR-corrected size.
 * 2. **Lower levels:** Merged bins are processed in parallel (OpenMP `taskloop`). Depending on
 *    do_I_need_a_fast_layout, each is laid out recursively with the fast layout or with the regular DP layout,
 *    and the result is grafted into `hibf_layout`.
 * 3. `max_bins` is sorted by level (path length), then lexicographically by path.
 *
 * The layout is not validated; execute validates it before writing it.
 */
void fast_layout(chopper::configuration const & config,
                 std::vector<size_t> const & positions,
                 std::vector<size_t> const & cardinalities,
                 std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                 std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                 seqan::hibf::layout::layout & hibf_layout);

} // namespace chopper::layout
