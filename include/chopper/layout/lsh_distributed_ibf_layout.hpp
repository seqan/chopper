// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::lsh_distributed_ibf_layout.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout
{

/*!\brief Distributes user bins onto `tmax` technical bins of a single IBF.
 * \param[in]  config            The configuration. Uses `tmax`, `sketch_bits`, `maximum_fpr`, `relaxed_fpr` and
 *                               `number_of_hash_functions` of `config.hibf_config`, and the layout timers.
 * \param[in]  positions         The indices of the user bins to distribute (need not be sorted).
 * \param[in]  cardinalities     The cardinality of each user bin, indexed by the values in `positions`.
 * \param[in]  sketches          The HyperLogLog sketch of each user bin, indexed by the values in `positions`.
 * \param[in]  minHash_sketches  The MinHash tables of each user bin, indexed by the values in `positions`.
 *                               Each table needs at least 3 sketches of at least 5 hashes each.
 * \param[out] technical_bins    Must have size `tmax` on entry. On return, `technical_bins[i]` holds the user bin
 *                               indices assigned to technical bin `i`.
 *
 * User bins are sorted by descending cardinality and handled in two groups:
 *
 * 1. **Split bins**: The largest user bins (cardinality above a threshold, see find_bins_to_be_split) are
 *    distributed over the technical bins at the *back* of `technical_bins` by determine_split_bins. At least one
 *    technical bin is left for merged bins. If fewer user bins remain than technical bins would be left, more
 *    technical bins are given to the split bins.
 * 2. **Merged bins**: The remaining user bins go into the technical bins at the *front* of `technical_bins`. They are
 *    clustered by similarity with MinHash LSH and assigned by HyperLogLog union estimates (see lsh_sim_approach).
 *    The target size per merged technical bin is `max(max_split_size / relaxed_fpr_correction, split_threshold)`.
 *
 * The first `split_threshold` is `ceil(joint_estimate / tmax)`, a lower bound that holds only if all user bins are
 * identical. If there are merged bins, the threshold is recalibrated **once**, based on the ratio of
 * `max_merged_size * relaxed_fpr_correction` to `max_split_size` (averaged with `max_merged_size` if nothing was
 * split). `technical_bins` is then cleared and both steps run again with the new threshold.
 */
void lsh_distributed_ibf_layout(chopper::configuration const & config,
                                std::vector<size_t> const & positions,
                                std::vector<size_t> const & cardinalities,
                                std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                std::vector<std::vector<size_t>> & technical_bins);

} // namespace chopper::layout
