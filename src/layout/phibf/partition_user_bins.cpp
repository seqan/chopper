// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/fast_layout_cluster.hpp>
#include <chopper/layout/lsh.hpp>
#include <chopper/layout/phibf/multi_cluster.hpp>
#include <chopper/layout/phibf/partition_user_bins.hpp>
#include <chopper/next_multiple_of_64.hpp>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>
#include <hibf/sketch/toolbox.hpp>
#include <hibf/misc/iota_vector.hpp>

namespace chopper::layout::phibf
{

/*!\brief Clusters similar user bins by iterative MinHash LSH (ADD+OR scheme).
 * \param[in] minHash_sketches           The MinHash tables of all user bins.
 * \param[in] cardinalities              The cardinality of each user bin.
 * \param[in] average_technical_bin_size Stop clustering as soon as the largest cluster's cardinality (sum of its user
 *                                       bins' cardinalities) reaches this.
 * \param[in] number_of_bands            LSH parameter b: the maximum number of rounds, one MinHash sketch per round.
 * \param[in] number_of_hashes_per_band  LSH parameter r: the number of hashes combined per LSH key.
 */
std::vector<Cluster> LSH_clustering(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                    std::vector<size_t> const & cardinalities,
                                    size_t const average_technical_bin_size,
                                    size_t const number_of_bands,
                                    size_t const number_of_hashes_per_band)
{

    size_t const number_of_user_bins{cardinalities.size()};
    size_t const number_of_max_minHash_sketches{number_of_bands};     // LSH ADD+OR parameter b
    size_t const minHash_sketche_size{number_of_hashes_per_band};     // LSH ADD+OR parameter r

    assert(number_of_user_bins == minHash_sketches.size());
    assert(!minHash_sketches.empty());
    assert(minHash_sketches[0].table.size() >= number_of_max_minHash_sketches);
    assert(minHash_sketches[0].table[0].size() >= minHash_sketche_size);

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
    size_t current_max_cluster_size{0};
    size_t current_number_of_sketch_hashes{minHash_sketche_size};
    size_t current_sketch_index{0};

    for (size_t user_bin_idx = 0; user_bin_idx < number_of_user_bins; ++user_bin_idx)
    {
        clusters.emplace_back(user_bin_idx);
        current_cluster_cardinality[user_bin_idx] = cardinalities[user_bin_idx];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[user_bin_idx]);
    }

    // refine clusters
    while (current_max_cluster_size < average_technical_bin_size
           && current_sketch_index < number_of_max_minHash_sketches)
    {

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // LSH_fill_hashtable returns sorted lists without duplicates.
            auto const begin = list.begin();
            auto const end = list.end();

            if (list.size() <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            size_t const representative_cluster_id = LSH_find_representative_cluster(clusters, *begin);
            auto & representative_cluster = clusters[representative_cluster_id];
            assert(representative_cluster.is_valid(representative_cluster_id));
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            for (auto current = begin + 1; current < end; ++current)
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                size_t const next_cluster_id = LSH_find_representative_cluster(clusters, *current);
                auto & next_cluster = clusters[next_cluster_id];

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.size() == 0);
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now
                assert(representative_cluster.is_valid(representative_cluster_id)); // and it should still be valid

                current_cluster_cardinality[representative_cluster.id()] += current_cluster_cardinality[next_cluster.id()];
            }

            current_max_cluster_size = *std::ranges::max_element(current_cluster_cardinality);
        }

        ++current_sketch_index;

    }

    return clusters;
}

// At the end of this function I have clusters of clusters (Multi clusters).
// E.g. [1:[A,B], 2:[C,D,E,F], 3:[G]]
// Let A-G be (initial_)clusters of very similar user bins obtained by previous clustering
// Then 1,2, and 3 Are a MultiCluster each. The clusters within each MultiCluster are rather similar
// Now clusters from 1,2 and 3 respectively can be used to initialise a partition because between the MultiClusters
// there is not much similarity
std::vector<MultiCluster> most_distant_LSH_partitioning(std::vector<Cluster> const & initial_clusters,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].table.empty());

    size_t const number_of_user_bins{initial_clusters.size()}; // empty clusters are still left in the vector
    size_t const number_of_max_minHash_sketches{minHash_sketches[0].table.size()}; // LSH ADD+OR parameter b
    size_t current_number_of_sketch_hashes{5};

    assert(number_of_user_bins == minHash_sketches.size());
    assert(minHash_sketches[0].table[0].size() >= current_number_of_sketch_hashes);

    size_t current_sketch_index{0};

    size_t current_number_of_clusters = [&initial_clusters] ()
    {
        size_t number{0};
        for (auto const & cluster : initial_clusters)
        {
            if (cluster.size())
                ++number;
        }
        return number;
    }();

    std::vector<MultiCluster> clusters;
    clusters.reserve(number_of_user_bins);

    for (size_t user_bin_idx = 0; user_bin_idx < number_of_user_bins; ++user_bin_idx)
    {
        clusters.emplace_back(initial_clusters[user_bin_idx]);
        assert(clusters[user_bin_idx].is_valid(user_bin_idx));
    }

    // refine clusters. I want to cluster as much as possible to initialise partitions with very distant clusters
    while (current_number_of_clusters > (config.hibf_config.tmax * 2) &&
           current_sketch_index < number_of_max_minHash_sketches)
    {

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // LSH_fill_hashtable returns sorted lists without duplicates.
            auto const begin = list.begin();
            auto const end = list.end();

            if (list.size() <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            size_t const representative_cluster_id = LSH_find_representative_cluster(clusters, *begin);
            auto & representative_cluster = clusters[representative_cluster_id];
            assert(representative_cluster.is_valid(representative_cluster_id));
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            for (auto current = begin + 1; current < end; ++current)
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                size_t const next_cluster_id = LSH_find_representative_cluster(clusters, *current);
                auto & next_cluster = clusters[next_cluster_id];

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.size() == 0);
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now
                assert(representative_cluster.is_valid(representative_cluster_id)); // and it should still be valid

                --current_number_of_clusters;
            }
        }

        ++current_sketch_index;

        // The LSH matching looser over time: Every 4 rounds, the number of hashes that must agree drops by one.
        // More and more dissimilar clusters merge, until there are at most 2·tmax clusters or the rounds run out.
        if (current_sketch_index % 4 == 0 && current_number_of_sketch_hashes > 1)
            --current_number_of_sketch_hashes;

    }

    return clusters;
}

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
    auto cluster_size_cmp = [&cardinalities](auto const & v1, auto const & v2)
    {
        if (v2.size() == v1.size() && !v2.empty())
            return  cardinalities[v2.contained_user_bins()[0]] < cardinalities[v1.contained_user_bins()[0]];
        return v2.size() < v1.size();
    };
    std::partial_sort(clusters.begin(), clusters.begin() + config.number_of_partitions, clusters.end(), cluster_size_cmp);

    // after filling up the partitions with the biggest clusters, sort the clusters by cardinality of the biggest ub
    // s.t. that euqally sizes ub are assigned after each other and the small stuff is added at last.
    // the largest ub is already at the start because of former sorting.
    auto compare_cardinality_and_move_empty_clusters_to_the_end = [&cardinalities](auto const & v1, auto const & v2)
    {
        if (v1.empty())
            return false; // v1 can never be larger than v2 then

        if (v2.empty()) // and v1 is not, since the first if would catch
            return true;

        return cardinalities[v2.contained_user_bins()[0]] < cardinalities[v1.contained_user_bins()[0]];
    };
    std::sort(clusters.begin() + config.number_of_partitions, clusters.end(), compare_cardinality_and_move_empty_clusters_to_the_end);

    assert(clusters.size() < 2 || clusters[0].size() >= clusters[1].size()); // sanity check
    // sanity check. Note that a cluster's id() is not its largest user bin; contained_user_bins()[0] is.
    assert(std::is_sorted(clusters.begin() + config.number_of_partitions,
                          clusters.end(),
                          compare_cardinality_and_move_empty_clusters_to_the_end));

// debug
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
// debug
}

void post_process_clusters(std::vector<MultiCluster> & clusters,
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
    auto cluster_size_cmp = [&cardinalities](auto const & mc1, auto const & mc2)
    {
        if (mc1.empty())
            return false; // mc1 can never be larger than mc2 then
        if (mc2.empty()) // and mc1 is not, since the first if would catch
            return true;

        auto const & largest_c1 = mc1.contained_user_bins()[0]; // first cluster is the largest because of sorting above
        auto const & largest_c2 = mc2.contained_user_bins()[0]; // first cluster is the largest because of sorting above

        assert(!largest_c1.empty()); // should never be empty
        assert(!largest_c2.empty()); // should never be empty

        if (largest_c2.size() == largest_c1.size())
            return  cardinalities[largest_c2[0]] < cardinalities[largest_c1[0]];
        return largest_c2.size() < largest_c1.size();
    };
    std::partial_sort(clusters.begin(), clusters.begin() + config.number_of_partitions, clusters.end(), cluster_size_cmp);

    auto move_empty_clusters_to_the_end = [](auto const & v1, auto const & v2)
    {
        return v2.size() < v1.size();
    };
    std::sort(clusters.begin() + config.number_of_partitions, clusters.end(), move_empty_clusters_to_the_end);

// debug sanity check
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
// debug
}

std::vector<MultiCluster> sim_dist_LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                      std::vector<size_t> const & cardinalities,
                                      size_t const average_technical_bin_size,
                                      chopper::configuration const & config)
{
    std::vector<Cluster> clusters = LSH_clustering(minHash_sketches,
                                                   cardinalities,
                                                   average_technical_bin_size,
                                                   3,                                     // b
                                                   minHash_sketches[0].table[0].size()); // r
    std::vector<MultiCluster> multi_clusters = most_distant_LSH_partitioning(clusters, minHash_sketches, config);
    post_process_clusters(multi_clusters, cardinalities, config);
    return multi_clusters;
}

bool find_best_partition(chopper::configuration const & config,
                         size_t const cardinality_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                         std::vector<size_t> & partition_cardinality)
{
    seqan::hibf::sketch::hyperloglog const current_sketch = [&sketches, &cluster, &config]()
    {
        seqan::hibf::sketch::hyperloglog result{config.hibf_config.sketch_bits};

        for (size_t const user_bin_idx : cluster)
            result.merge(sketches[user_bin_idx]);

        return result;
    }();

    size_t const current_cardinality = [&cardinalities, &cluster]()
    {
        size_t result{0};

        for (size_t const user_bin_idx : cluster)
            result += cardinalities[user_bin_idx];

        return result;
    }();

    // edge case: empty sequence
    if (current_cardinality == 0) // no kmers to be added
    {
        // just add empty user bin to the first partition
        for (size_t const user_bin_idx : cluster)
        {
            positions[0].push_back(user_bin_idx);
            partition_cardinality[0] += cardinalities[user_bin_idx];
        }
        partition_sketches[0].merge(current_sketch);
        return true;
    }

    // TODO: afterwads check if I should again merge by how little the effective text ratio grows

    // search best partition fit by similarity
    // similarity here is defined as:
    // "whose (<-partition) effective text size is subsumed most by the current user bin"
    // or in other words:
    // "which partition has the largest intersection with user bin b compared to its own (partition) size."
    double best_subsume_ratio{0.0};
    size_t best_p{0};
    bool best_p_found{false};

    for (size_t p = 0; p < config.number_of_partitions; ++p)
    {
        seqan::hibf::sketch::hyperloglog tmp = current_sketch;
        tmp.merge(partition_sketches[p]);
        size_t const union_size = tmp.estimate();
        size_t const current_partition_size = partition_sketches[p].estimate();
        assert(union_size >= current_partition_size); // partition cannot shrink by adding kmers
        size_t const change = union_size - current_partition_size;
        size_t const intersection = current_sketch.estimate() - change;
        double const subsume_ratio = static_cast<double>(intersection) / current_partition_size;

        if (subsume_ratio > best_subsume_ratio &&
            partition_cardinality[p] < cardinality_per_part &&
            (current_cardinality + partition_cardinality[p]) < (cardinality_per_part * 1.2) )
        {
            best_subsume_ratio = subsume_ratio;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        return false;

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        partition_cardinality[best_p] += cardinalities[user_bin_idx];
    }
    partition_sketches[best_p].merge(current_sketch);

    return true;
}

/*!\brief Assigns a single user bin with find_best_partition.
 * \throws std::runtime_error If the user bin does not fit into any partition.
 */
void assign_user_bin(chopper::configuration const & config,
                     size_t const cardinality_per_part,
                     size_t const user_bin_idx,
                     std::vector<size_t> const & cardinalities,
                     std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                     std::vector<std::vector<size_t>> & partitions,
                     std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                     std::vector<size_t> & partition_cardinality)
{
    if (!find_best_partition(config,
                             cardinality_per_part,
                             {user_bin_idx},
                             cardinalities,
                             sketches,
                             partitions,
                             partition_sketches,
                             partition_cardinality))
    {
        throw std::runtime_error{"User bin " + std::to_string(user_bin_idx) + " (cardinality "
                                 + std::to_string(cardinalities[user_bin_idx])
                                 + ") does not fit into any partition (target cardinality per partition: "
                                 + std::to_string(cardinality_per_part)
                                 + "). Use fewer partitions or another partitioning approach."};
    }
}

void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<std::vector<size_t>> & partitions)
{
    if (config.partitioning_approach < partitioning_scheme::blocked
        || config.partitioning_approach > partitioning_scheme::lsh_sim)
    {
        throw std::invalid_argument{"Unknown partitioning approach " + std::to_string(config.partitioning_approach)
                                    + "."};
    }

    if (cardinalities.size() < config.number_of_partitions)
    {
        throw std::invalid_argument{"The number of partitions (" + std::to_string(config.number_of_partitions)
                                    + ") must not exceed the number of user bins ("
                                    + std::to_string(cardinalities.size()) + ")."};
    }

    // all approaches need sorted positions
    std::vector<size_t> const sorted_positions = [&cardinalities]()
    {
        std::vector<size_t> ps;
        ps.resize(cardinalities.size());
        std::iota(ps.begin(), ps.end(), 0);
        seqan::hibf::sketch::toolbox::sort_by_cardinalities(cardinalities, ps);
        return ps;
    }();

    size_t const sum_of_cardinalities = [&cardinalities]()
    {
        size_t sum{0};
        for (size_t const card : cardinalities)
            sum += card;
        return sum;
    }();

    if (config.partitioning_approach == partitioning_scheme::blocked)
    {
        size_t const u_bins_per_part = cardinalities.size() / config.number_of_partitions;
        // for small edge cases, just use u_bins_per_part,
        // otherwise, use next_multiple_of_64(ceil(sqrt(u_bins_per_part))) to have am ibf friendly block size.
        // That way, we loop over the partitions, filling them block wise, again and again.
        size_t const block_size = (u_bins_per_part < 64) ? u_bins_per_part :
                                   chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part))));

        size_t current_part{0u};
        size_t current_block_count{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            partitions[current_part].push_back(current_user_bin_id);
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
        size_t const cardinality_per_part = sum_of_cardinalities /  config.number_of_partitions;

        size_t current_cardinality{0u};
        size_t current_part{0};

        for (size_t i = 0; i < sorted_positions.size(); ++i)
        {
            partitions[current_part].push_back(sorted_positions[i]);
            current_cardinality += cardinalities[sorted_positions[i]];

            if (current_cardinality >= cardinality_per_part /*partition is full*/||
                (sorted_positions.size() - i <= config.number_of_partitions - current_part) /* small number edge case:not enough UBs left*/)
            {
                ++current_part;
                // If cardinality_per_part was not chosen optimally
                // we might need start at position 0 again.
                // If all paritions are full, this will add one user bin at a time
                // to each partition until no user bins remain.
                if (current_part >= config.number_of_partitions)
                    current_part = 0;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::folded)
    {
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
            partitions[parts[current_part]].push_back(current_user_bin_id);
            current_cardinality += cardinalities[current_user_bin_id];

            if (current_cardinality >= cardinality_per_part_halved)
            {
                current_cardinality = 0;
                // Remaining user bins (with cardinality 0) stay in the last part.
                if (current_part + 1 < parts.size())
                    ++current_part;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::weighted_fold)
    {
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);

        // The unassigned user bins are sorted_positions[big, small_end). Large user bins are taken from the front,
        // small ones from the back.
        size_t big{0};
        size_t small_end{sorted_positions.size()};

        // How far a partition with the given cardinality and number of user bins is off the targets. Smaller is better.
        auto compute_score = [&](size_t const cardinality, size_t const number_of_user_bins)
        {
            double const weight = static_cast<double>(cardinality) / cardinality_per_part;
            double const amount = static_cast<double>(number_of_user_bins) / u_bins_per_part;
            return std::abs(1.0 - weight) + std::abs(1.0 - amount);
        };

        for (size_t current_part = 0; current_part + 1 < config.number_of_partitions; ++current_part)
        {
            std::vector<size_t> large_bins;
            std::vector<size_t> small_bins;
            size_t current_cardinality{0};

            // first add all large bins that fit
            while (current_cardinality < cardinality_per_part && big < small_end)
            {
                large_bins.push_back(sorted_positions[big]);
                current_cardinality += cardinalities[sorted_positions[big]];
                ++big;
            }

            double local_optimum = compute_score(current_cardinality, large_bins.size());

            // then remove large bins and add small bins until a local optimum is reached
            while (!large_bins.empty())
            {
                // remove the last large bin ...
                size_t trial_cardinality = current_cardinality - cardinalities[large_bins.back()];
                size_t trial_number_of_user_bins = large_bins.size() - 1 + small_bins.size();
                size_t trial_small_end = small_end;
                double trial_score = compute_score(trial_cardinality, trial_number_of_user_bins);

                // ... and add small bins as long as this improves the score. The removed large bin is
                // sorted_positions[big - 1] and must not be added as a small bin.
                while (trial_small_end > big)
                {
                    size_t const candidate = sorted_positions[trial_small_end - 1];
                    double const score =
                        compute_score(trial_cardinality + cardinalities[candidate], trial_number_of_user_bins + 1);

                    if (score >= trial_score) // smaller is better
                        break;

                    trial_cardinality += cardinalities[candidate];
                    ++trial_number_of_user_bins;
                    --trial_small_end;
                    trial_score = score;
                }

                if (trial_number_of_user_bins == 0 || local_optimum < trial_score) // score would increase. Stop
                    break;

                // update
                large_bins.pop_back();
                --big;
                for (size_t pos = small_end; pos > trial_small_end; --pos)
                    small_bins.push_back(sorted_positions[pos - 1]);
                small_end = trial_small_end;
                current_cardinality = trial_cardinality;
                local_optimum = trial_score;
            }

            partitions[current_part] = std::move(large_bins);
            partitions[current_part].insert(partitions[current_part].end(), small_bins.begin(), small_bins.end());
        }

        // remaining user bins go to last partition
        for (; big < small_end; ++big)
            partitions[config.number_of_partitions - 1].push_back(sorted_positions[big]);
    }
    else if (config.partitioning_approach == partitioning_scheme::similarity)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const block_size = [&]()
        {
            size_t const u_bins_per_part = cardinalities.size() / config.number_of_partitions;
            return std::min(u_bins_per_part,
                            chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));
        }();
        size_t const number_of_blocks = seqan::hibf::divide_and_ceil(cardinalities.size(), block_size);

        // don't move from largest to smallest but pick the next block to process randomly.
        // this probably leads to more evenly distributed partitions (evenly in terms of number of user bins)
        std::vector<size_t> indices(number_of_blocks);
        std::iota(indices.begin(), indices.end(), 0);
        std::random_device shuffle_random_device;
        std::mt19937 shuffle_engine(shuffle_random_device());
        std::shuffle(indices.begin(), indices.end(), shuffle_engine);

        // initialise partitions with the first random config.number_of_partitions blocks
        assert(number_of_blocks >= config.number_of_partitions);
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            size_t const i = indices[p];
            size_t const end = (i == (indices.size() - 1) && (cardinalities.size() % block_size != 0)) ? (cardinalities.size() % block_size) : block_size;
            for (size_t x = 0; x < end; ++x)
            {
                assert(block_size * i + x < sorted_positions.size());
                partition_sketches[p].merge(sketches[sorted_positions[block_size * i + x]]);
                partition_cardinality[p] += cardinalities[sorted_positions[block_size * i + x]];
                partitions[p].push_back(sorted_positions[block_size * i + x]);
            }
        }

        auto has_been_processed_in_init = [&config,&indices, &block_size](size_t const i)
        {
            bool result{false};
            for (size_t p = 0; p < config.number_of_partitions; ++p)
                result |= ((indices[p] * block_size) == i);
            return result;
        };

        // assign the rest by similarity, user in by user bin (single bins, no blocks)
        for (size_t i = 0; i < sorted_positions.size(); ++i)
        {
            if (has_been_processed_in_init(i))
            {
                i += (block_size - 1); // -1 because there will be an increment after continue
                continue;
            }

            assign_user_bin(config,
                            cardinality_per_part,
                            sorted_positions[i],
                            cardinalities,
                            sketches,
                            partitions,
                            partition_sketches,
                            partition_cardinality);
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        // initial partitioning using locality sensitive hashing (LSH)
        std::vector<Cluster> raw_clusters = LSH_clustering(minHash_sketches,
                                                           cardinalities,
                                                           cardinality_per_part,
                                                           seqan::hibf::sketch::minhashes::num_sketches,  // b
                                                           seqan::hibf::sketch::minhashes::sketch_size); // r
        post_process_clusters(raw_clusters, cardinalities, config);
        std::vector<Cluster> const & clusters = raw_clusters; // clusters do not change from this point on

        // initialise partitions with the first p largest clusters (post_processing sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            // There may be fewer (non-empty) clusters than partitions. The remaining partitions stay empty, which is
            // reported at the end.
            if (clusters[cidx].empty())
                break;

            bool p_has_been_incremented{false};

            for (size_t uidx = 0; uidx < clusters[cidx].size(); ++uidx)
            {
                size_t const user_bin_idx = clusters[cidx].contained_user_bins()[uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if (partition_cardinality[p] > cardinality_per_part)
                {
                    ++p;
                    p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sketches[p].merge(sketches[user_bin_idx]);
                partition_cardinality[p] += cardinalities[user_bin_idx];
                partitions[p].push_back(user_bin_idx);
            }

            if (p_has_been_incremented && partition_cardinality[p] < cardinality_per_part)
            {
                // p will be incremeted again in the next iteration of this for loop
                // since we have incremented p previously because of some user bins that spilled over in this currnet p
                // we want to decrease p here, s.t. we stay at the same p.
                // this ensures that the current p doesn't only have a small rest that just didn't fit into the
                //  previous partition
                assert(p > 0);
                --p;
            }

            ++cidx;
        }

        // assign the rest by similarity
        assert(config.hibf_config.number_of_user_bins == clusters.size());
        for (; cidx < config.hibf_config.number_of_user_bins; ++cidx)
        {
            auto const & cluster = clusters[cidx];

            if (cluster.empty()) // non valid clusters have been removed. Since list was sorted, we can abort
                break;

            bool const cluster_could_be_added = find_best_partition(config, cardinality_per_part, cluster.contained_user_bins(), cardinalities, sketches, partitions, partition_sketches, partition_cardinality);

            // if the cluster is too large to be added. Add single user bins
            if (!cluster_could_be_added)
            {
                for (size_t const user_bin_idx : cluster.contained_user_bins())
                {
                    assign_user_bin(config, cardinality_per_part, user_bin_idx, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);
                }
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh_sim)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        // initial partitioning using locality sensitive hashing (LSH)
        std::vector<MultiCluster> const clusters = sim_dist_LSH_partitioning(minHash_sketches, cardinalities, cardinality_per_part, config);

        // initialise partitions with the first p largest clusters (post processing sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            // There may be fewer (non-empty) clusters than partitions. The remaining partitions stay empty, which is
            // reported at the end.
            if (clusters[cidx].empty())
                break;

            bool p_has_been_incremented{false};
            auto const & cluster = clusters[cidx].contained_user_bins()[0];

            for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
            {
                size_t const user_bin_idx = cluster[uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if (partition_cardinality[p] > cardinality_per_part)
                {
                    ++p;
                    p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sketches[p].merge(sketches[user_bin_idx]);
                partition_cardinality[p] += cardinalities[user_bin_idx];
                partitions[p].push_back(user_bin_idx);
            }

            if (p_has_been_incremented && partition_cardinality[p] < cardinality_per_part)
            {
                // p will be incremeted again in the next iteration of this for loop
                // since we have incremented p previously because of some user bins that spilled over in this currnet p
                // we want to decrease p here, s.t. we stay at the same p.
                // this ensures that the current p doesn't only have a small rest that just didn't fit into the
                //  previous partition
                assert(p > 0);
                --p;
            }

            ++cidx;
        }

        // assign the rest, cluster by cluster (because the clusters are very similar), by similarity
        assert(config.hibf_config.number_of_user_bins == clusters.size());

        std::vector<std::vector<size_t>> remaining_clusters{};

        for (size_t i = 0; i < clusters.size(); ++i)
        {
            if (clusters[i].empty())
                break;
            // if i < cidx the first cluster was already used in initialising above
            auto start = clusters[i].contained_user_bins().begin() + ((i < cidx) ? 1 : 0);
            remaining_clusters.insert(remaining_clusters.end(), start, clusters[i].contained_user_bins().end());
        }
        auto cmp = [&cardinalities] (auto const & c1, auto const & c2) { return cardinalities[c2[0]] < cardinalities[c1[0]]; };
        std::sort(remaining_clusters.begin(), remaining_clusters.end(), cmp);

        for (auto const & cluster : remaining_clusters)
        {
            bool const cluster_could_be_added = find_best_partition(config, cardinality_per_part, cluster, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);

            // if the cluster is too large to be added. Add single user bins
            if (!cluster_could_be_added)
            {
                for (size_t const user_bin_idx : cluster)
                {
                    assign_user_bin(config, cardinality_per_part, user_bin_idx, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);
                }
            }
        }
    }

#ifndef NDEBUG
    // sanity check in debug
    std::vector<size_t> user_bins{};
    for (auto const & partition : partitions) {
        user_bins.insert(user_bins.end(), partition.begin(), partition.end());
    }
    if (!std::ranges::is_permutation(user_bins, seqan::hibf::iota_vector(sketches.size())))
    {
        std::string str{"Not all user bins have been assigned to the "};
        str += std::to_string(partitions.size());
        str += " partitions! (";
        str += std::to_string(user_bins.size());
        str += "/";
        str += std::to_string(sketches.size());
        str += ")\n";
        for (auto const & p : partitions)
        {
            str += "[";
            for (auto const h : p)
            {
                str += std::to_string(h);
                str += ",";
            }
            str.back() = ']';
            str += '\n';
        }

        throw std::logic_error{str};
    }
#endif

    // Each partition is laid out as a separate HIBF, which needs at least one user bin.
    for (size_t p = 0; p < partitions.size(); ++p)
    {
        if (partitions[p].empty())
        {
            throw std::runtime_error{"Partition " + std::to_string(p) + " of " + std::to_string(partitions.size())
                                     + " is empty. Use fewer partitions or another partitioning approach."};
        }
    }
}

} // namespace chopper::layout::phibf
