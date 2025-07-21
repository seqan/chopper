// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <chopper/layout/determine_split_bins.hpp>
#include <chopper/layout/fast_layout_cluster.hpp>
#include <chopper/layout/fast_layout_find_bins_to_be_split.hpp>
#include <chopper/layout/partition_user_bins.hpp>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/sketch/toolbox.hpp>

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
uint64_t lsh_hash_the_sketch(std::vector<uint64_t> const & sketch, size_t const number_of_hashes_to_consider)
{
    assert(number_of_hashes_to_consider <= sketch.size());
    return std::reduce(sketch.begin(), sketch.begin() + number_of_hashes_to_consider);
}

/*!\brief Builds the LSH collision table of the current clusters for one LSH band.
 * \param[in] clusters                        The current clusters. Clusters that were moved are skipped.
 * \param[in] minHash_sketches                The MinHash tables of all user bins, indexed by global user bin index.
 * \param[in] current_sketch_index            The sketch (row of the MinHash table) to use (LSH band index).
 * \param[in] current_number_of_sketch_hashes The number of hashes combined per key (LSH parameter r).
 * \returns A map from LSH key to the sorted, unique ids of the representative clusters that produced the key.
 *
 * Each user bin in a valid cluster adds its key, and the cluster's id is stored under that key. A multi-member
 * cluster can therefore appear under several keys, so clusters that share a key with *any* member collide.
 */
auto LSH_fill_hashtable(std::vector<Cluster> const & clusters,
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

        for (size_t const user_bin_idx : current.contained_user_bins())
        {
            ++processed_user_bins;
            uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx].table[current_sketch_index],
                                                     current_number_of_sketch_hashes);
            table[key].push_back(current.id()); // insert representative for all user bins
        }
    }
    assert(processed_user_bins == clusters.size()); // all user bins should've been processed by one of the clusters

    // uniquify list. Since I am inserting representative_idx's into the table, the same number can
    // be inserted into multiple splots, and multiple times in the same slot.
    for (auto & [key, list] : table)
    {
        std::ranges::sort(list);
        auto const ret = std::ranges::unique(list);
        list.erase(ret.begin(), ret.end());
    }

    return table;
}

/*!\brief Clusters very similar user bins by iterative MinHash LSH.
 * \param[in] minHash_sketches            The MinHash tables of all user bins, indexed by global user bin index.
 * \param[in] positions                   The global indices of the user bins to cluster.
 * \param[in] cardinalities               The cardinality of each user bin, indexed by global user bin index.
 * \param[in] sketches                    The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in] average_technical_bin_size  The threshold. Stop clustering as soon as a cluster's estimated cardinality
 *                                        reaches this.
 * \param[in] config                      The configuration (uses `hibf_config.sketch_bits`).
 * \returns One Cluster per entry in `positions`. Cluster `i` is created with id `i` (a local index) and holds user bin
 *          `positions[i]` (a global index). After merging, a cluster is either valid (`id() == i`, contains >= 1
 *          user bins) or moved (empty, and `moved_to_cluster_id()` points to the cluster it was merged into,
 *          possibly through a chain of moves; resolve with LSH_find_representative_cluster).
 *
 * In each round, one LSH band (`current_sketch_index`) is used to build a collision table (LSH_fill_hashtable), and
 * all clusters in a bucket are merged into the representative of the first one (OR step). The merged cluster's
 * cardinality is re-estimated from the union of its HyperLogLog sketches. Rounds stop after
 * `number_of_max_minHash_sketches` (b = 3) bands or when the largest cluster cardinality reaches
 * `average_technical_bin_size`. The number of hashes per key stays at `minHash_sketch_size` (r = 5) in all rounds.
 * b and r were chosen by experiment.
 */
std::vector<Cluster> very_similar_LSH_clustering(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                                 std::vector<size_t> const & positions,
                                                 std::vector<size_t> const & cardinalities,
                                                 std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                                 size_t const average_technical_bin_size,
                                                 chopper::configuration const & config)
{
    // The following two parameters are experimentally derived:
    size_t const number_of_max_minHash_sketches{3}; // LSH ADD+OR parameter b
    size_t const minHash_sketch_size{5};            // LSH ADD+OR parameter r
    size_t const number_of_user_bins{positions.size()};
    seqan::hibf::sketch::hyperloglog const empty_sketch{config.hibf_config.sketch_bits};

    assert(!minHash_sketches.empty());
    assert(minHash_sketches[0].table.size() >= number_of_max_minHash_sketches);
    assert(minHash_sketches[0].table[0].size() >= minHash_sketch_size);
    assert(number_of_user_bins <= minHash_sketches.size());

    // initialise clusters with a signle user bin per cluster.
    // clusters are either
    // 1) of size 1; containing an id != position where the id points to the cluster it has been moved to
    //    e.g. cluster[Y] = {Z} (Y has been moved into Z, so Z could look likes this cluster[Z] = {Z, Y})
    // 2) of size >= 1; with the first entry beging id == position (a valid cluster)
    //    e.g. cluster[X] = {X}       // valid singleton
    //    e.g. cluster[X] = {X, a, b, c, ...}   // valid cluster with more joined entries
    // The clusters could me moved recursively, s.t.
    // cluster[A] = {B}
    // cluster[B] = {C}
    // cluster[C] = {C, A, B} // is valid cluster since cluster[C][0] == C; contains A and B
    std::vector<Cluster> clusters;
    clusters.reserve(number_of_user_bins);

    std::vector<size_t> current_cluster_cardinality(number_of_user_bins);
    std::vector<seqan::hibf::sketch::hyperloglog> current_cluster_sketches(number_of_user_bins, empty_sketch);
    size_t current_max_cluster_size{0};
    size_t current_sketch_index{0};

    for (size_t pos = 0; pos < number_of_user_bins; ++pos)
    {
        clusters.emplace_back(pos, positions[pos]);
        current_cluster_cardinality[pos] = cardinalities[positions[pos]];
        current_cluster_sketches[pos] = sketches[positions[pos]];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[positions[pos]]);
    }

    auto get_cluster = [&clusters](size_t const current) -> Cluster &
    {
        size_t const cluster_id = LSH_find_representative_cluster(clusters, current);
        return clusters[cluster_id];
    };

    // refine clusters
    while (current_max_cluster_size < average_technical_bin_size
           && current_sketch_index < number_of_max_minHash_sketches)
    {
        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, minHash_sketch_size);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            if (list.size() <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            auto & representative_cluster = get_cluster(list[0]);
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            auto & representative_cluster_sketch = current_cluster_sketches[representative_cluster.id()];

            for (size_t const current : std::views::drop(list, 1))
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                auto & next_cluster = get_cluster(current);

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.empty());
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now

                representative_cluster_sketch.merge(current_cluster_sketches[next_cluster.id()]);
            }

            current_cluster_cardinality[representative_cluster.id()] = representative_cluster_sketch.estimate();
        }
        current_max_cluster_size = std::ranges::max(current_cluster_cardinality);

        ++current_sketch_index;
    }

    return clusters;
}

/*!\brief Orders the clusters so that lsh_sim_approach can take the leading ones as partition seeds.
 * \param[in,out] clusters      The clusters returned by very_similar_LSH_clustering. Reordered in place.
 * \param[in]     cardinalities The cardinality of each user bin, indexed by global user bin index.
 * \param[in]     config        The configuration (uses `hibf_config.tmax`).
 * \throws std::runtime_error if an empty cluster ends up before a non-empty one (sanity check).
 *
 * 1. The user bins inside each cluster are sorted by descending cardinality, so `contained_user_bins().front()` is
 *    the largest.
 * 2. The first `tmax` positions receive the largest clusters by number of user bins, in descending order. Ties are
 *    broken by the cardinality of the largest user bin.
 * 3. The remaining clusters are sorted by the cardinality of their largest user bin, in descending order. Empty
 *    (moved) clusters go last.
 *
 * After this, a cluster's position no longer matches its id(), so moved-to links (and
 * LSH_find_representative_cluster) must not be used anymore.
 */
void post_process_clusters(std::vector<Cluster> & clusters,
                           std::vector<size_t> const & cardinalities,
                           chopper::configuration const & config)
{
    // clusters are done. Start post processing
    // since post processing involves re-ordering the clusters, the moved_to_cluster_id value of a cluster will not
    // refer to the position of the cluster in the `clusters` vecto anymore but the cluster with the resprive id()
    // would neet to be found
    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        assert(clusters[pos].is_valid(pos));
        clusters[pos].sort_by_cardinality(cardinalities);
    }

    // push largest p clusters to the front
    std::ranges::partial_sort(clusters,
                              std::ranges::next(clusters.begin(), config.hibf_config.tmax, clusters.end()),
                              [&cardinalities](auto const & v1, auto const & v2)
                              {
                                  // Note: If v2 is empty, so is v1.
                                  if (v1.size() == v2.size() && !v2.empty())
                                      return cardinalities[v1.contained_user_bins().front()]
                                           > cardinalities[v2.contained_user_bins().front()];

                                  return v1.size() > v2.size();
                              });

    // after filling up the partitions with the biggest clusters, sort the clusters by cardinality of the biggest ub
    // s.t. that euqally sizes ub are assigned after each other and the small stuff is added at last.
    // the largest ub is already at the start because of former sorting.
    std::ranges::sort(std::ranges::next(clusters.begin(), config.hibf_config.tmax, clusters.end()),
                      clusters.end(),
                      [&cardinalities](auto const & v1, auto const & v2)
                      {
                          if (v1.empty())
                              return false; // v1 can never be larger than v2 then

                          if (v2.empty()) // and v1 is not, since the first if would catch
                              return true;

                          return cardinalities[v1.contained_user_bins().front()]
                               > cardinalities[v2.contained_user_bins().front()];
                      });

    assert(clusters.size() < 2 || clusters[0].size() >= clusters[1].size()); // sanity check

#ifndef NDEBUG
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        // once empty - always empty; all empty clusters should be at the end
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty())
            throw std::runtime_error{"sorting did not work"};
    }
#endif
}

/*!\brief Assigns a whole cluster of user bins to the partition where adding it costs the least.
 * \param[in]     config                      The configuration (uses `hibf_config.tmax` and `sketch_bits`).
 * \param[in]     number_of_partitions        Only partitions `[0, number_of_partitions)` are considered.
 * \param[in,out] corrected_estimate_per_part The current target cardinality per partition. Raised to the chosen
 *                                            partition's new estimate if that is larger.
 * \param[in]     cluster                     The global user bin indices to assign together.
 * \param[in]     cardinalities               The cardinality of each user bin, indexed by global user bin index.
 * \param[in]     sketches                    The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in,out] positions                   The user bins per partition. `cluster` is appended to the chosen one.
 * \param[in,out] partition_sketches          The union sketch per partition. Updated for the chosen partition.
 * \param[in,out] max_partition_cardinality   The largest user bin cardinality per partition. Updated.
 * \param[in,out] min_partition_cardinality   The smallest user bin cardinality per partition. Updated.
 * \returns `true`. Throws std::runtime_error if no partition was selected, i.e., if `number_of_partitions == 0`.
 *
 * For each partition p, the cost of adding the cluster is the sum of
 * - `union - |p|`: the new k-mers the cluster adds to p (HyperLogLog estimates),
 * - `tmax * max(0, union - corrected_estimate_per_part)`: growth of the IBF bin size beyond the current target, and
 * - a lower-level penalty depending on `max_card`, the largest user bin cardinality in `cluster`:
 *   - p already holds more than `tmax` user bins (lower level exists): `max_card * log_tmax(#UBs after adding)`,
 *     an estimate of how often the content is stored again on lower levels;
 *   - adding the cluster pushes p above `tmax` user bins (new lower level): `min(min_p, max_card) * tmax`;
 *   - otherwise: `max_p - max_card` if the cluster is smaller than every user bin in p (wasted space), or
 *     `(max_card - max_p) * tmax` if it is larger than every user bin in p (IBF grows), else 0.
 *
 * The partition with the smallest cost is chosen. A partition with zero cost always replaces the current best, so if
 * several have zero cost, the last one wins.
 */
bool find_best_partition(chopper::configuration const & config,
                         size_t const number_of_partitions,
                         size_t & corrected_estimate_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                         std::vector<size_t> & max_partition_cardinality,
                         std::vector<size_t> & min_partition_cardinality)
{
    seqan::hibf::sketch::hyperloglog const current_sketch = [&sketches, &cluster, &config]()
    {
        seqan::hibf::sketch::hyperloglog result{config.hibf_config.sketch_bits};

        for (size_t const user_bin_idx : cluster)
            result.merge(sketches[user_bin_idx]);

        return result;
    }();

    size_t const max_card = [&cardinalities, &cluster]()
    {
        size_t max{0};

        for (size_t const user_bin_idx : cluster)
            max = std::max(max, cardinalities[user_bin_idx]);

        return max;
    }();

    // Search best partition fit by similarity. Similarity here is defined as:
    // "whose (<-partition) effective text size is subsumed most by the current user bin". Or in other words:
    // "which partition has the largest intersection with user bin b compared to its own (partition) size."
    size_t smallest_change{std::numeric_limits<size_t>::max()};
    size_t best_p{0};
    bool best_p_found{false};

    auto penalty_lower_level = [&](size_t const additional_number_of_user_bins, size_t const p) -> size_t
    {
        assert(positions[p].size() != 0); // partitions should be initialised beforehand
        size_t const min = min_partition_cardinality[p];
        size_t const max = max_partition_cardinality[p];

        if (positions[p].size() > config.hibf_config.tmax) // already a third level
        {
            // if there must already be another lower level because the current merged bin contains more than tmax
            // user bins, then the current user bin is very likely stored multiple times. Therefore, the penalty is set
            // to the cardinality of the current user bin times the number of levels, e.g. the number of times this user
            // bin needs to be stored additionally
            size_t const num_ubs_in_merged_bin{positions[p].size() + additional_number_of_user_bins};
            double const levels = std::log(num_ubs_in_merged_bin) / std::log(config.hibf_config.tmax);
            return static_cast<size_t>(max_card * levels);
        }
        else if (positions[p].size() + additional_number_of_user_bins > config.hibf_config.tmax) // now a third level
        {
            // if the current merged bin contains exactly tmax UBS, adding otherone must
            // result in another lower level. Most likely, the smallest user bin will end up on the lower level
            // therefore the penalty is set to 'min * tmax'
            // of course, there could also be a third level with a lower number of user bins, but this is hard to
            // estimate.
            size_t const penalty = std::min(min, max_card) * config.hibf_config.tmax;
            return penalty;
        }
        else // positions[p].size() + additional_number_of_user_bins <= tmax
        {
            // if the new user bin is smaller than all other already contained user bins
            // the waste of space is high if stored in a single technical bin
            if (max_card < min)
                return (max - max_card);
            // if the new user bin is bigger than all other already contained user bins, the IBF size increases
            else if (max_card > max)
                return (max_card - max) * config.hibf_config.tmax;
            // else, end if-else-block and zero is returned
        }

        return 0u;
    };

    for (size_t p = 0; p < number_of_partitions; ++p)
    {
        seqan::hibf::sketch::hyperloglog union_sketch = current_sketch;
        union_sketch.merge(partition_sketches[p]);
        size_t const union_estimate = union_sketch.estimate();
        size_t const current_partition_size = partition_sketches[p].estimate();

        assert(union_estimate >= current_partition_size);
        size_t const penalty_current_bin = union_estimate - current_partition_size;
        size_t const penalty_current_ibf =
            config.hibf_config.tmax
            * ((union_estimate <= corrected_estimate_per_part) ? 0u : union_estimate - corrected_estimate_per_part);
        size_t const change = penalty_current_bin + penalty_current_ibf + penalty_lower_level(cluster.size(), p);

        if (change == 0 || /* If there is no penalty at all, this is a best fit even if the partition is "full"*/
            (smallest_change > change))
        {
            smallest_change = change;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        throw std::runtime_error{"currently there are no safety measures if a partition is not found"};

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        max_partition_cardinality[best_p] = std::max(max_partition_cardinality[best_p], cardinalities[user_bin_idx]);
        min_partition_cardinality[best_p] = std::min(min_partition_cardinality[best_p], cardinalities[user_bin_idx]);
    }
    partition_sketches[best_p].merge(current_sketch);
    corrected_estimate_per_part = std::max<size_t>(corrected_estimate_per_part, partition_sketches[best_p].estimate());

    return true;
}

/*!\brief Distributes the merged-bin candidates onto `number_of_remaining_tbs` partitions by LSH clustering and
 *        similarity-based assignment.
 * \param[in]     config                      The configuration (uses `hibf_config` and the LSH/search timers).
 * \param[in]     sorted_positions2           The global indices of the user bins to distribute, sorted by descending
 *                                            cardinality (the remainder after split bins were removed).
 * \param[in]     cardinalities               The cardinality of each user bin, indexed by global user bin index.
 * \param[in]     sketches                    The HyperLogLog sketch of each user bin, indexed by global user bin index.
 * \param[in]     minHash_sketches            The MinHash tables of each user bin, indexed by global user bin index.
 * \param[in,out] partitions                  Receives the assignment in `partitions[0, number_of_remaining_tbs)`.
 *                                            Must have at least `number_of_remaining_tbs` entries.
 * \param[in]     number_of_remaining_tbs     The number of technical bins available for merged bins.
 * \param[in]     technical_bin_size_threshold The target cardinality per technical bin. Stops LSH clustering and
 *                                            triggers spill-over while seeding partitions.
 * \param[in]     sum_of_cardinalities        The sum of the cardinalities of *all* user bins of this IBF.
 * \returns The largest estimated cardinality of any of the `number_of_remaining_tbs` partitions.
 *
 * 1. **Cluster:** very_similar_LSH_clustering followed by post_process_clusters.
 * 2. **Ensure enough clusters:** If there are fewer non-empty clusters than partitions, user bins are moved out of
 *    the last cluster with more than one user bin, into their own clusters, until there are enough clusters.
 * 3. **Seed:** Partition p receives the next cluster in order. A cluster is spread over consecutive partitions
 *    whenever `tmax` user bins have been placed or the partition's estimate exceeds `technical_bin_size_threshold`.
 *    A cluster with more than `tmax` user bins whose total cardinality exceeds `0.05 * sum_of_cardinalities / tmax`
 *    is placed in full. Any other cluster places at most `tmax` user bins. User bins not placed, or left over
 *    because the partitions ran out, go to a list of remaining clusters.
 * 4. **Assign the rest:** Each remaining cluster, plus all clusters that were not used as seeds, is assigned as a
 *    whole by find_best_partition. The per-partition target starts at `technical_bin_size_threshold` and only grows.
 */
size_t lsh_sim_approach(chopper::configuration const & config,
                        std::vector<size_t> const & sorted_positions2,
                        std::vector<size_t> const & cardinalities,
                        std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        std::vector<std::vector<size_t>> & partitions,
                        size_t const number_of_remaining_tbs,
                        size_t const technical_bin_size_threshold,
                        size_t const sum_of_cardinalities)
{
    uint8_t const sketch_bits{config.hibf_config.sketch_bits};
    std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(number_of_remaining_tbs,
                                                                     seqan::hibf::sketch::hyperloglog(sketch_bits));

    std::vector<size_t> max_partition_cardinality(number_of_remaining_tbs, 0u);
    std::vector<size_t> min_partition_cardinality(number_of_remaining_tbs, std::numeric_limits<size_t>::max());

    // initial partitioning using locality sensitive hashing (LSH)
    config.lsh_algorithm_timer.start();
    std::vector<Cluster> clusters = very_similar_LSH_clustering(minHash_sketches,
                                                                sorted_positions2,
                                                                cardinalities,
                                                                sketches,
                                                                technical_bin_size_threshold,
                                                                config);
    post_process_clusters(clusters, cardinalities, config);
    config.lsh_algorithm_timer.stop();

    // There must be more non-empty clusters than technical bins
    size_t number_of_remaining_clusters = std::ranges::count_if(clusters,
                                                                [](auto const & c)
                                                                {
                                                                    return !c.empty();
                                                                });
    clusters.resize(std::max<size_t>(number_of_remaining_tbs, clusters.size()));
    while (number_of_remaining_tbs > number_of_remaining_clusters)
    {
        // If there are not enough remaining clusters, we need to split up clusters. This is rather an edge case.
        auto empty_cluster_it = clusters.begin() + number_of_remaining_clusters;
        // We can only split up clusters of size > 1. Start with smaller ones to the back of `clusters`
        auto cluster_it = std::ranges::find_if(clusters.rbegin(),
                                               clusters.rend(),
                                               [](auto const & c)
                                               {
                                                   return c.size() > 1;
                                               });
        assert(cluster_it != clusters.rend());
        for (; cluster_it->size() != 1 && number_of_remaining_tbs > number_of_remaining_clusters;)
        {
            assert(empty_cluster_it->empty());
            empty_cluster_it->add_user_bin(cluster_it->pop_back()); // not super efficient but fine
            ++number_of_remaining_clusters;
            ++empty_cluster_it;
        }
    }
    assert(number_of_remaining_tbs <= static_cast<size_t>(std::ranges::count_if(clusters,
                                                                                [](auto const & c)
                                                                                {
                                                                                    return !c.empty();
                                                                                })));

    std::vector<std::vector<size_t>> remaining_clusters{};

    // initialise partitions with the first p largest clusters (post_processing sorts by size)
    size_t cidx{0}; // current cluster index
    for (size_t p = 0; p < number_of_remaining_tbs; ++p)
    {
        assert(!clusters[cidx].empty());
        auto const & cluster = clusters[cidx].contained_user_bins();
        bool split_cluster = false;

        if (cluster.size() > config.hibf_config.tmax)
        {
            size_t card{0};
            for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
                card += cardinalities[cluster[uidx]];

            if (card > 0.05 * sum_of_cardinalities / config.hibf_config.tmax)
                split_cluster = true;
        }

        size_t end = (split_cluster) ? cluster.size() : std::min(cluster.size(), config.hibf_config.tmax);
        for (size_t uidx = 0; uidx < end; ++uidx)
        {
            size_t const user_bin_idx = cluster[uidx];
            // if a single cluster already exceeds the cardinality_per_part,
            // then the remaining user bins of the cluster must spill over into the next partition
            if ((uidx != 0 && (uidx % config.hibf_config.tmax == 0))
                || partition_sketches[p].estimate() > technical_bin_size_threshold)
            {
                ++p;

                if (p >= number_of_remaining_tbs)
                {
                    split_cluster = true;
                    end = uidx;
                    break;
                }
            }

            partition_sketches[p].merge(sketches[user_bin_idx]);
            partitions[p].push_back(user_bin_idx);
            max_partition_cardinality[p] = std::max(max_partition_cardinality[p], cardinalities[user_bin_idx]);
            min_partition_cardinality[p] = std::min(min_partition_cardinality[p], cardinalities[user_bin_idx]);
        }

        if (split_cluster)
        {
            std::vector<size_t> remainder(cluster.begin() + end, cluster.end());
            remaining_clusters.insert(remaining_clusters.end(), remainder);
        }

        ++cidx;
    }

    for (size_t i = cidx; i < clusters.size(); ++i)
    {
        if (clusters[i].empty())
            break;

        remaining_clusters.insert(remaining_clusters.end(), clusters[i].contained_user_bins());
    }

    // assign the rest by similarity
    size_t merged_threshold{technical_bin_size_threshold};
    for (size_t ridx = 0; ridx < remaining_clusters.size(); ++ridx)
    {
        auto const & cluster = remaining_clusters[ridx];

        config.search_partition_algorithm_timer.start();
        find_best_partition(config,
                            number_of_remaining_tbs,
                            merged_threshold,
                            cluster,
                            cardinalities,
                            sketches,
                            partitions,
                            partition_sketches,
                            max_partition_cardinality,
                            min_partition_cardinality);
        config.search_partition_algorithm_timer.stop();
    }

    // compute actual max size
    size_t max_size{0};
    for (auto const & sketch : partition_sketches)
        max_size = std::max(max_size, (size_t)sketch.estimate());

    return max_size;
}

void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & positions,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<std::vector<size_t>> & partitions)
{
    // all approaches need sorted positions
    std::vector<size_t> const sorted_positions = [&positions, &cardinalities]()
    {
        std::vector<size_t> ps(positions.begin(), positions.end());
        seqan::hibf::sketch::toolbox::sort_by_cardinalities(cardinalities, ps);
        return ps;
    }();

    auto const [sum_of_cardinalities, joint_estimate] = [&]()
    {
        size_t sum{0};
        seqan::hibf::sketch::hyperloglog sketch{config.hibf_config.sketch_bits};

        for (size_t const pos : positions)
        {
            sum += cardinalities[pos];
            sketch.merge(sketches[pos]);
        }

        return std::tuple<size_t, size_t>{sum, sketch.estimate()};
    }();

    double const relaxed_fpr_correction = seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config.hibf_config.maximum_fpr, //
         .relaxed_fpr = config.hibf_config.relaxed_fpr,
         .hash_count = config.hibf_config.number_of_hash_functions});

    size_t idx{0}; // start in sorted positions
    size_t number_of_split_tbs{0};
    size_t number_of_merged_tbs{config.hibf_config.tmax};
    size_t max_split_size{0};
    size_t max_merged_size{0};

    // Start with the lower bound on the split threshold: All user bin content is the same and when split/merged
    // one technical bin is never more than joint_estimate/tmax per bin.
    // (unrealistic but a good starting point. threshold will be revised in a second iteration)
    size_t split_threshold = seqan::hibf::divide_and_ceil(joint_estimate, config.hibf_config.tmax);

    auto partition_split_bins = [&]()
    {
        size_t number_of_potential_split_bins{0};    // determined in find_bins_to_be_split
        size_t max_tbs{config.hibf_config.tmax - 1}; // leave one bin for merging

        std::tie(idx, number_of_potential_split_bins) =
            find_bins_to_be_split(sorted_positions, cardinalities, split_threshold, max_tbs);

        // if (remaining UBs < remaining TBs for merging)
        // then assign more split bins, such that there are at least as manu UBs as TBs for merging
        if (sorted_positions.size() - idx < config.hibf_config.tmax - number_of_potential_split_bins)
            number_of_potential_split_bins +=
                (config.hibf_config.tmax - number_of_potential_split_bins) - (sorted_positions.size() - idx);

        std::tie(number_of_split_tbs, max_split_size) =
            chopper::layout::determine_split_bins(config,
                                                  sorted_positions,
                                                  cardinalities,
                                                  number_of_potential_split_bins,
                                                  idx,
                                                  partitions);
        number_of_merged_tbs = config.hibf_config.tmax - number_of_split_tbs;
    };

    auto partition_merged_bins = [&]()
    {
        // determine number of split bins
        std::vector<size_t> const sorted_positions2(sorted_positions.begin() + idx, sorted_positions.end());

        // distribute the rest to merged bins
        size_t const corrected_max_split_size = max_split_size / relaxed_fpr_correction;
        size_t const merged_threshold = std::max(corrected_max_split_size, split_threshold);

        max_merged_size = lsh_sim_approach(config,
                                           sorted_positions2,
                                           cardinalities,
                                           sketches,
                                           minHash_sketches,
                                           partitions,
                                           number_of_merged_tbs,
                                           merged_threshold,
                                           sum_of_cardinalities);
    };

    partition_split_bins();

    // All user bins can be assigned as split bins (idx == sorted_positions.size()).
    // In that case there are no remaining bins to distribute via partition_merged_bins.
    // And no reconfiguration of the threshold needs to be done since splitting is done with an "optimal" DP
    if (idx < sorted_positions.size())
    {
        partition_merged_bins();

        int64_t const difference =
            static_cast<int64_t>(max_merged_size * relaxed_fpr_correction) - static_cast<int64_t>(max_split_size);

        if (number_of_split_tbs == 0)
            split_threshold = (split_threshold + max_merged_size) / 2; // increase threshold
        else if (difference > 0)                                       // need more merged bins -> increase threshold
            split_threshold = static_cast<double>(split_threshold)
                            * ((static_cast<double>(max_merged_size) * relaxed_fpr_correction)
                               / static_cast<double>(max_split_size));
        else // need more split bins -> decrease threshold
            split_threshold = std::max<double>(1.0,
                                               static_cast<double>(split_threshold)
                                                   * ((static_cast<double>(max_merged_size) * relaxed_fpr_correction)
                                                      / static_cast<double>(max_split_size)));

        // reset result
        partitions.clear();
        partitions.resize(config.hibf_config.tmax);
        idx = 0;
        number_of_split_tbs = 0;
        number_of_merged_tbs = config.hibf_config.tmax;
        max_split_size = 0;
        max_merged_size = 0;

        partition_split_bins();
        // All user bins can be assigned as split bins (idx == sorted_positions.size()).
        // In that case there are no remaining bins to distribute via partition_merged_bins.
        if (idx < sorted_positions.size())
            partition_merged_bins();
    }
}

} // namespace chopper::layout
