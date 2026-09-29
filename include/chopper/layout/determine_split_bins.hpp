// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::determine_split_bins.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>

namespace chopper::layout
{

/*!\brief Distributes the largest user bins as split bins over technical bins at the back of `partitions`.
 * \param[in]     config             The configuration (uses `maximum_fpr` and `number_of_hash_functions` of
 *                                   `config.hibf_config`).
 * \param[in]     positions          User bin indices, sorted by descending cardinality. The first `num_user_bins`
 *                                   are split.
 * \param[in]     cardinalities      The cardinality of each user bin, indexed by the values in `positions`.
 * \param[in]     num_technical_bins The maximum number of technical bins for the split user bins. Must be at least
 *                                   `num_user_bins`.
 * \param[in]     num_user_bins      The number of user bins to split.
 * \param[in,out] partitions         The technical bins. Split user bin indices are appended to the last technical
 *                                   bins, one index per technical bin. `positions[num_user_bins - 1]` gets the last
 *                                   technical bins, `positions[0]` the first of them. The technical bins of each user
 *                                   bin are consecutive.
 * \returns A pair of
 *          1. the number of technical bins used, which may be less than `num_technical_bins`, and
 *          2. the largest FPR-corrected cardinality per technical bin.
 *          `{0, 0}` if `num_user_bins == 0`.
 *
 * A dynamic programming algorithm assigns each user bin a number of technical bins, such that the largest
 * FPR-corrected cardinality per technical bin, `ceil(cardinality * fpr_correction[k] / k)` for a user bin in `k`
 * technical bins, is minimal. If several numbers of technical bins give the same minimum, the smallest is used.
 */
std::pair<size_t, size_t> determine_split_bins(chopper::configuration const & config,
                                               std::vector<size_t> const & positions,
                                               std::vector<size_t> const & cardinalities,
                                               size_t const num_technical_bins,
                                               size_t const num_user_bins,
                                               std::vector<std::vector<size_t>> & partitions);

} // namespace chopper::layout
