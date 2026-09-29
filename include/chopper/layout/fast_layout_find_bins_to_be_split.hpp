// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::find_bins_to_be_split.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace chopper::layout
{

/*!\brief Determines which user bins are split and how many technical bins the split user bins occupy.
 * \param[in] sorted_positions User bin indices, sorted by descending cardinality
 *                             (see seqan::hibf::sketch::toolbox::sort_by_cardinalities).
 * \param[in] cardinalities    The cardinality of each user bin, indexed by the values in `sorted_positions`.
 * \param[in] threshold        The initial cardinality threshold. User bins with a cardinality greater than
 *                             `threshold` are split. Must be greater than 0.
 * \param[in] max_bins         The maximum number of technical bins that split user bins may occupy.
 * \returns A pair of
 *          1. the number of user bins to split, i.e., the length of the prefix of `sorted_positions`, and
 *          2. the number of technical bins for these split user bins, which is at most `max_bins`.
 *
 * The user bins with a cardinality greater than `threshold` form a prefix of `sorted_positions`. The number of
 * technical bins they need is `sum / threshold` (rounded down), where `sum` is the sum of their cardinalities.
 * If this exceeds `max_bins`, `threshold` is set to `max(threshold + 1, sum / max_bins)` (rounded down) and the
 * prefix is computed again, until the split user bins fit into `max_bins` technical bins.
 *
 * The number of split user bins never exceeds the number of technical bins, so each split user bin gets at least one
 * technical bin: Each split user bin has a cardinality greater than `threshold`, hence `sum / threshold` is at least
 * the number of split user bins. This is asserted, not enforced.
 */
inline std::pair<size_t, size_t> find_bins_to_be_split(std::vector<size_t> const & sorted_positions,
                                                       std::vector<size_t> const & cardinalities,
                                                       size_t threshold,
                                                       size_t const max_bins)
{
    assert(threshold > 0);
    assert(max_bins > 0);

    size_t idx{0};
    size_t sum{0};
    // update idx and sum
    auto find_idx_and_sum = [&]()
    {
        while (idx < sorted_positions.size() && cardinalities[sorted_positions[idx]] > threshold)
        {
            sum += cardinalities[sorted_positions[idx]];
            ++idx;
        }
    };
    find_idx_and_sum();

    // SInce the threshold is the expected size of each technical bin, the number of split bins needed i:
    size_t number_of_split_bins = sum / threshold;

    // If number_of_split_bins is more than the available technical bins, the threshold must be adjusted
    while (number_of_split_bins > max_bins)
    {
        // Since number_of_split_bins is too large iwth the given threshold
        // we need to increase the threshold such that less user bins are targetted for splitting.
        // The new threshold is therefore adjusted to the new expected average technical bin size of
        // distributing the split bin content evenly to `max_bin` technical bins.
        // max(threshold, ...) to avoid an endless loop.
        threshold = std::max<size_t>(threshold + 1, static_cast<double>(sum) / max_bins);

        // update idx, sum and number_of_split_bins
        idx = 0;
        sum = 0;
        find_idx_and_sum();
        number_of_split_bins = sum / threshold;
    }

    assert(idx <= number_of_split_bins); // there should never be more user bins than available split bins

    return {idx, number_of_split_bins};
}

} // namespace chopper::layout
