// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the MinHash LSH primitives shared by the fast layout and the partitioned HIBF.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <ranges>
#include <vector>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout
{

/*!\brief Combines the first `number_of_hashes_to_consider` MinHash values of a sketch into a single LSH key.
 * \param[in] sketch                       A single MinHash sketch (one row of seqan::hibf::sketch::minhashes::table).
 * \param[in] number_of_hashes_to_consider The number of leading hashes to combine (LSH parameter r).
 *                                         Must be `<= sketch.size()`.
 * \returns The sum of the first `number_of_hashes_to_consider` hashes, with unsigned wrap-around.
 *
 * This is the AND step of the LSH AND-OR scheme: two user bins get the same key only if all r hashes agree,
 * apart from sum collisions.
 */
inline uint64_t lsh_hash_the_sketch(std::vector<uint64_t> const & sketch, size_t const number_of_hashes_to_consider)
{
    assert(number_of_hashes_to_consider <= sketch.size());
    return std::reduce(sketch.begin(), sketch.begin() + number_of_hashes_to_consider);
}

/*!\brief Builds the LSH collision table of the current clusters for one LSH band.
 * \tparam    cluster_t                       The cluster type. Either chopper::layout::Cluster, whose
 *                                            `contained_user_bins()` is a range of user bins, or a type whose
 *                                            `contained_user_bins()` is a range of ranges of user bins (e.g., the
 *                                            MultiCluster of the partitioned HIBF).
 * \param[in] clusters                        The current clusters. Clusters that were moved are skipped.
 * \param[in] minHash_sketches                The MinHash tables of all user bins, indexed by global user bin index.
 * \param[in] current_sketch_index            The sketch (row of the MinHash table) to use (LSH band index).
 * \param[in] current_number_of_sketch_hashes The number of hashes combined per key (LSH parameter r).
 * \returns A map from LSH key to the sorted, unique ids of the representative clusters that produced the key.
 *
 * Each user bin in a valid cluster adds its key, and the cluster's id is stored under that key. A multi-member
 * cluster can therefore appear under several keys, so clusters that share a key with *any* member collide.
 */
template <typename cluster_t>
auto LSH_fill_hashtable(std::vector<cluster_t> const & clusters,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        size_t const current_sketch_index,
                        size_t const current_number_of_sketch_hashes)
{
    robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

    [[maybe_unused]] size_t processed_user_bins{0}; // only for sanity check

    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        auto const & current = clusters[pos];
        assert(current.is_valid(pos));

        if (current.has_been_moved()) // cluster has been moved somewhere else, don't process
            continue;

        auto insert = [&](size_t const user_bin_idx)
        {
            ++processed_user_bins;
            uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx].table[current_sketch_index],
                                                     current_number_of_sketch_hashes);
            table[key].push_back(current.id()); // insert representative for all user bins
        };

        using user_bins_t = std::remove_cvref_t<decltype(current.contained_user_bins())>;

        if constexpr (std::ranges::range<std::ranges::range_value_t<user_bins_t>>) // e.g. MultiCluster
        {
            for (auto const & similarity_cluster : current.contained_user_bins())
                for (size_t const user_bin_idx : similarity_cluster)
                    insert(user_bin_idx);
        }
        else
        {
            for (size_t const user_bin_idx : current.contained_user_bins())
                insert(user_bin_idx);
        }
    }
    assert(processed_user_bins == clusters.size()); // all user bins should've been processed by one of the clusters

    // uniquify list. Since I am inserting representative_idx's into the table, the same number can
    // be inserted into multiple splots, and multiple times in the same slot.
    for (auto & [key, list] : table)
    {
        std::ranges::sort(list);
        auto const [first, last] = std::ranges::unique(list);
        list.erase(first, last);
    }

    return table;
}

} // namespace chopper::layout
