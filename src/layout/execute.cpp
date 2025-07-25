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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/determine_split_bins.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/lsh.hpp>
#include <chopper/next_multiple_of_64.hpp>
#include <chopper/sketch/output.hpp>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp> // for compute_relaxed_fpr_correction
#include <hibf/layout/layout.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/toolbox.hpp>

namespace chopper::layout
{

uint64_t lsh_hash_the_sketch(std::vector<uint64_t> const & sketch, size_t const number_of_hashes_to_consider)
{
    // lets just compute the sum
    uint64_t sum{0};

    for (size_t i = 0; i < number_of_hashes_to_consider; ++i)
        sum += sketch[i];

    return sum;
}

auto LSH_fill_hashtable(std::vector<Cluster> const & clusters,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        size_t const current_sketch_index,
                        size_t const current_number_of_sketch_hashes)
{
    robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

    size_t processed_user_bins{0}; // only for sanity check
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

    return table;
}

auto LSH_fill_hashtable(std::vector<MultiCluster> const & clusters,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        size_t const current_sketch_index,
                        size_t const current_number_of_sketch_hashes)
{
    robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

    size_t processed_user_bins{0}; // only for sanity check
    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        auto const & current = clusters[pos];
        assert(current.is_valid(pos));

        if (current.has_been_moved()) // cluster has been moved somewhere else, don't process
            continue;

        for (auto const & similarity_cluster : current.contained_user_bins())
        {
            for (size_t const user_bin_idx : similarity_cluster)
            {
                ++processed_user_bins;
                uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx].table[current_sketch_index],
                                                         current_number_of_sketch_hashes);
                table[key].push_back(current.id()); // insert representative for all user bins
            }
        }
    }
    assert(processed_user_bins == clusters.size()); // all user bins should've been processed by one of the clusters

    return table;
}

// minHash_sketches data structure:
// Vector L1 : number of user bins
// Vector L2 : number_of_max_minHash_sketches (LSH ADD+OR parameter b)
// Vector L3 : minHash_sketche_size (LSH ADD+OR parameter r)

std::vector<Cluster> very_similar_LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                                   std::vector<size_t> const & positions,
                                                   std::vector<size_t> const & cardinalities,
                                                   std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                                   size_t const average_technical_bin_size,
                                                   chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].table.empty());
    assert(!minHash_sketches[0].table[0].empty());

    size_t const number_of_user_bins{positions.size()};
    assert(number_of_user_bins <= minHash_sketches.size());
    size_t const number_of_max_minHash_sketches{3}; // LSH ADD+OR parameter b
    // size_t const minHash_sketche_size{minHash_sketches[0].table[0].size()};   // LSH ADD+OR parameter r
    size_t const minHash_sketche_size{5}; // LSH ADD+OR parameter r
    seqan::hibf::sketch::hyperloglog const empty_sketch{config.hibf_config.sketch_bits};
    // std::cout << "sketch size available: " << minHash_sketches[0].table[0].size() << std::endl;
    // std::cout << "sketch size used here: " << minHash_sketche_size << std::endl;
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
    size_t current_number_of_sketch_hashes{minHash_sketche_size}; // start with high r but decrease it iteratively
    size_t current_sketch_index{0};
    [[maybe_unused]] size_t current_number_of_clusters{number_of_user_bins}; // initially, each UB is a separate cluster

    for (size_t pos = 0; pos < number_of_user_bins; ++pos)
    {
        clusters.emplace_back(pos, positions[pos]);
        current_cluster_cardinality[pos] = cardinalities[positions[pos]];
        current_cluster_sketches[pos] = sketches[positions[pos]];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[positions[pos]]);
    }

    // refine clusters
    //std::cout << "Start clustering with threshold average_technical_bin_size: " << average_technical_bin_size << std::endl;
    while (current_max_cluster_size < average_technical_bin_size
           && /*number_of_clusters / static_cast<double>(number_of_user_bins) > 0.5 &&*/
           current_sketch_index < number_of_max_minHash_sketches) // I want to cluster 10%?
    {
        //std::cout << "Current number of clusters: " << current_number_of_clusters;

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // uniquify list. Since I am inserting representative_idx's into the table, the same number can
            // be inserted into multiple splots, and multiple times in the same slot.
            std::sort(list.begin(), list.end());
            auto const end = std::unique(list.begin(), list.end());
            auto const begin = list.begin();

            if (end - begin <= 1) // nothing to do here
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

                current_cluster_sketches[representative_cluster.id()].merge(
                    current_cluster_sketches[next_cluster.id()]);
                current_cluster_cardinality[representative_cluster.id()] =
                    current_cluster_sketches[representative_cluster.id()].estimate();

                --current_number_of_clusters;
            }

            current_max_cluster_size = *std::ranges::max_element(current_cluster_cardinality);
        }

        ++current_sketch_index;

        //std::cout << " and after this clustering step there are: " << current_number_of_clusters << "with max cluster size" << current_max_cluster_size << std::endl;
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
            return cardinalities[v2.contained_user_bins()[0]] < cardinalities[v1.contained_user_bins()[0]];
        return v2.size() < v1.size();
    };
    std::partial_sort(clusters.begin(),
                      clusters.begin() + std::min<size_t>(clusters.size(), config.hibf_config.tmax),
                      clusters.end(),
                      cluster_size_cmp);

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
    std::sort(clusters.begin() + std::min<size_t>(clusters.size(), config.hibf_config.tmax),
              clusters.end(),
              compare_cardinality_and_move_empty_clusters_to_the_end);

    assert(clusters[0].size() >= clusters[1].size()); // sanity check
    // assert(cardinalities[clusters[std::min<size_t>(clusters.size(), config.hibf_config.tmax)].contained_user_bins()[0]] >= cardinalities[clusters[std::min<size_t>(clusters.size(), config.hibf_config.tmax) + 1].contained_user_bins()[0]]); // sanity check

    // debug
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty()
            && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
    // debug
}

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

    // TODO: afterwads check if I should again merge by how little the effective text ratio grows

    // search best partition fit by similarity
    // similarity here is defined as:
    // "whose (<-partition) effective text size is subsumed most by the current user bin"
    // or in other words:
    // "which partition has the largest intersection with user bin b compared to its own (partition) size."
    // double best_subsume_ratio{0.0};
    size_t smallest_change{std::numeric_limits<size_t>::max()};
    size_t best_p{0};
    bool best_p_found{false};

    auto penalty_lower_level = [&](size_t const additional_number_of_user_bins, size_t const p)
    {
        size_t min = std::min(max_card, min_partition_cardinality[p]);
        size_t max = std::max(max_card, max_partition_cardinality[p]);

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
            size_t const penalty = min * config.hibf_config.tmax;
            return penalty;
        }
        else // positions[p].size() + additional_number_of_user_bins < tmax
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

        return (size_t)0u;
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

        // size_t const intersection = current_sketch.estimate() - change;
        // double const subsume_ratio = static_cast<double>(intersection) / current_partition_size;
        //std::cout << "p:" << p << " p-#UBs" << positions[p].size() << " penalty:" <<  penalty(cluster.size(), p) << " change:" << change << " union-current_p:" << (union_estimate - current_partition_size) << " union:" << union_estimate << " current_p:" << current_partition_size << " t:" << corrected_estimate_per_part << std::endl;
        if (change == 0 || /* If there is no penalty at all, this is a best fit even if the partition is "full"*/
            (smallest_change > change /*&& subsume_ratio > best_subsume_ratio &&*/
             /*current_partition_size < corrected_estimate_per_part*/))
        {
            //std::cout << "smaller!" << std::endl;
            // best_subsume_ratio = subsume_ratio;
            smallest_change = change;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        throw "currently there are no safety measures if a partition is not found because it is very unlikely";

    //std::cout << "best_p:" << best_p << std::endl<< std::endl;

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        max_partition_cardinality[best_p] = std::max(max_partition_cardinality[best_p], cardinalities[user_bin_idx]);
        min_partition_cardinality[best_p] = std::min(min_partition_cardinality[best_p], cardinalities[user_bin_idx]);
    }
    partition_sketches[best_p].merge(current_sketch);
    corrected_estimate_per_part =
        std::max(corrected_estimate_per_part, static_cast<size_t>(partition_sketches[best_p].estimate()));

    return true;
}

size_t split_bins(chopper::configuration const & config,
                  std::vector<size_t> const & sorted_positions,
                  std::vector<size_t> const & cardinalities,
                  std::vector<std::vector<size_t>> & partitions,
                  size_t const sum_of_cardinalities,
                  size_t const end_idx)
{
    // assign split bins to the end. makes it easier to process merged bins afterwards
    size_t pos{partitions.size() - 1};
    for (size_t idx = 0; idx < end_idx; ++idx)
    {
        size_t const ub_idx{sorted_positions[idx]};
        size_t const number_of_split_tbs =
            std::max((size_t)1, config.hibf_config.tmax * cardinalities[ub_idx] / sum_of_cardinalities);
        std::cout << "number_of_split_tbs: " << number_of_split_tbs
                  << " cardinalities[ub_idx]:" << cardinalities[ub_idx]
                  << " sum_of_cardinalities:" << sum_of_cardinalities << std::endl;
        // fill partitions from behind to ensure an easier layouting
        for (size_t i = 0; i < number_of_split_tbs; ++i)
        {
            assert(pos > 0);
            assert(partitions[pos].empty());
            partitions[pos].push_back(ub_idx);
            --pos;
        }
    }
    std::cout << "pos: " << pos << std::endl;
    return pos + 1;
}

std::pair<size_t, size_t> find_bins_to_be_split(std::vector<size_t> const & sorted_positions,
                                                std::vector<size_t> const & cardinalities,
                                                size_t threshold,
                                                size_t const max_bins)
{
    size_t idx{0};
    size_t sum{0};
    size_t number_of_split_bins = max_bins + 1; // initialise to something bigger to trigger while loop below

    auto find_idx_and_sum = [&]()
    {
        while (idx < sorted_positions.size() && cardinalities[sorted_positions[idx]] > threshold)
        {
            sum += cardinalities[sorted_positions[idx]];
            ++idx;
        }
    };

    while (number_of_split_bins > max_bins)
    {
        idx = 0;
        sum = 0;
        find_idx_and_sum();
        number_of_split_bins = sum / threshold;
        // max(1.01, ...) because a tie would end in an endless loop.
        threshold = static_cast<double>(threshold) * std::max(1.01, static_cast<double>(sum) / (threshold * max_bins));
    }

    // Maybe there are more user bins (position of idx) than available split bins.
    // Thus, clamp the number of user bins to split to maximally the number_of_split_bins
    idx = std::min(idx, number_of_split_bins);

    return {idx, number_of_split_bins};
}

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
    //std::cout << "LSH partitioning into " << config.hibf_config.tmax << std::endl;
    std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(number_of_remaining_tbs,
                                                                     seqan::hibf::sketch::hyperloglog(sketch_bits));

    std::vector<size_t> max_partition_cardinality(number_of_remaining_tbs, 0u);
    std::vector<size_t> min_partition_cardinality(number_of_remaining_tbs, std::numeric_limits<size_t>::max());

    // initial partitioning using locality sensitive hashing (LSH)
    config.lsh_algorithm_timer.start();
    std::vector<Cluster> clusters = very_similar_LSH_partitioning(minHash_sketches,
                                                                  sorted_positions2,
                                                                  cardinalities,
                                                                  sketches,
                                                                  technical_bin_size_threshold,
                                                                  config);
    post_process_clusters(clusters, cardinalities, config);
    config.lsh_algorithm_timer.stop();

    std::ofstream ofs{"/tmp/final.clusters"};
    for (size_t i = 0; i < clusters.size(); ++i)
    {
        seqan::hibf::sketch::hyperloglog sketch(config.hibf_config.sketch_bits);
        for (size_t j = 0; j < clusters[i].size(); ++j)
        {
            sketch.merge(sketches[clusters[i].contained_user_bins()[j]]);
        }
        ofs << i << ":" << sketch.estimate() << ":" << clusters[i].size() << std::endl;
    }

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
        config.search_partition_algorithm_timer.start();
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

    auto const [sum_of_cardinalities, max_cardinality, joint_estimate] = [&]()
    {
        size_t sum{0};
        size_t max{0};
        seqan::hibf::sketch::hyperloglog sketch{config.hibf_config.sketch_bits};

        for (size_t const pos : positions)
        {
            sum += cardinalities[pos];
            max = std::max(max, cardinalities[pos]);
            sketch.merge(sketches[pos]);
        }

        return std::tuple<size_t, size_t, size_t>{sum, max, sketch.estimate()};
    }();

    // If the effective text size is very low, it can happen that the joint_estimate divided by the number of partitions
    // is lower than the largest single user bin. But of course, we can never reach a smaller max technical bin size
    // then that of the largest user user bin. Thus we can correct the estimate_per_part beforehand.
    // This way we make sure there is at least 1 LSH clustering step.
    [[maybe_unused]] size_t const estimate_per_part =
        std::max(seqan::hibf::divide_and_ceil(joint_estimate, config.hibf_config.tmax), max_cardinality + 1);

    double const relaxed_fpr_correction = seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config.hibf_config.maximum_fpr, //
         .relaxed_fpr = config.hibf_config.relaxed_fpr,
         .hash_count = config.hibf_config.number_of_hash_functions});
    //std::cout << "sum_of_cardinalities:" << sum_of_cardinalities << " joint_estimate:" << joint_estimate << std::endl;

    size_t split_threshold{};
    size_t idx{0}; // start in sorted positions
    size_t number_of_split_tbs{0};
    size_t number_of_merged_tbs{config.hibf_config.tmax};
    size_t max_split_size{0};
    size_t max_merged_size{0};

    // split_threshold = seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.hibf_config.tmax);
    split_threshold = seqan::hibf::divide_and_ceil(joint_estimate, config.hibf_config.tmax);

    auto parition_split_bins = [&]()
    {
        size_t number_of_potential_split_bins{0};
        std::tie(idx, number_of_potential_split_bins) =
            find_bins_to_be_split(sorted_positions, cardinalities, split_threshold, config.hibf_config.tmax - 1);
        std::tie(number_of_split_tbs, max_split_size) =
            chopper::layout::determine_split_bins(config,
                                                  sorted_positions,
                                                  cardinalities,
                                                  number_of_potential_split_bins,
                                                  idx,
                                                  partitions);
        number_of_merged_tbs = config.hibf_config.tmax - number_of_split_tbs;
    };

    auto partitions_merged_bins = [&]()
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

    parition_split_bins();
    partitions_merged_bins();

    int64_t const difference =
        static_cast<int64_t>(max_merged_size * relaxed_fpr_correction) - static_cast<int64_t>(max_split_size);

    std::cout << "number_of_split_tbs:" << number_of_split_tbs << " difference:" << difference << std::endl;
    std::cout << "Reconfiguring threshold.  from:" << split_threshold;

    if (number_of_split_tbs == 0)
        split_threshold = (split_threshold + max_merged_size) / 2; // increase threshold
    else if (difference > 0)                                       // need more merged bins -> increase threshold
        split_threshold =
            static_cast<double>(split_threshold)
            * ((static_cast<double>(max_merged_size) * relaxed_fpr_correction) / static_cast<double>(max_split_size));
    else // need more split bins -> decrease threshold
        split_threshold =
            static_cast<double>(split_threshold)
            * ((static_cast<double>(max_merged_size) * relaxed_fpr_correction) / static_cast<double>(max_split_size));

    std::cout << " to:" << split_threshold << std::endl;

    // reset result
    partitions.clear();
    partitions.resize(config.hibf_config.tmax);
    idx = 0;
    number_of_split_tbs = 0;
    number_of_merged_tbs = config.hibf_config.tmax;
    max_split_size = 0;
    max_merged_size = 0;

    parition_split_bins();
    partitions_merged_bins();

    // sanity check:
    size_t sum{0};
    size_t last_index{positions.size()}; // non existing user bin idx
    for (auto const & p : partitions)
    {
        if (p.size() == 1)
        {
            // for split bins, avoid additional counts
            if (p[0] != last_index)
                ++sum;

            last_index = p[0];
        }
        else
        {
            sum += p.size();
        }
    }

    if (sum != positions.size())
    {
        std::string str{"Not all user bins have been assigned to the "};
        str += std::to_string(partitions.size());
        str += " partitions! (";
        str += std::to_string(sum);
        str += "/";
        str += std::to_string(positions.size());
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

    // assert([&](){ bool x{false}; for (auto const & p : partitions) { x &= !p.empty(); }; return x; });
}

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
    auto config_copy = config;

    std::vector<std::vector<size_t>> tmax_partitions(config.hibf_config.tmax);

    // here we assume that we want to start with a fast layout
    config.intital_partition_timer.start();
    partition_user_bins(config_copy, positions, cardinalities, sketches, minHash_sketches, tmax_partitions);
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
    hibf_layout.user_bins.resize(config_copy.hibf_config.number_of_user_bins);

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
#pragma omp parallel
#pragma omp single
    {
#pragma omp taskloop
        for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
        {
            auto const & partition = tmax_partitions[partition_idx];

            if (partition.empty() || partition.size() == 1) // nothing to merge
                continue;

            if (do_I_need_a_fast_layout(config_copy, partition, cardinalities))
            {
                fast_layout_recursion(config_copy,
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
}

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches)
{
    config.hibf_config.validate_and_set_defaults();

    std::vector<size_t> cardinalities;
    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);

    if (true) // 0 == unset == single HIBF, 1 == single HIBF
    {
        seqan::hibf::layout::layout hibf_layout;

        if (config.determine_best_tmax)
        {
            hibf_layout = determine_best_number_of_technical_bins(config, cardinalities, sketches);
        }
        else
        {
            config.dp_algorithm_timer.start();
            fast_layout(config,
                        seqan::hibf::iota_vector(sketches.size()),
                        cardinalities,
                        sketches,
                        minHash_sketches,
                        hibf_layout);
            config.dp_algorithm_timer.stop();

            // hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
            //                                                   cardinalities,
            //                                                   sketches,
            //                                                   seqan::hibf::iota_vector(sketches.size()),
            //                                                   config.union_estimation_timer,
            //                                                   config.rearrangement_timer);

            // sort records ascending by the number of bin indices (corresponds to the IBF levels)
            // GCOVR_EXCL_START
            std::ranges::sort(hibf_layout.max_bins,
                              [](auto const & r, auto const & l)
                              {
                                  if (r.previous_TB_indices.size() == l.previous_TB_indices.size())
                                      return std::ranges::lexicographical_compare(r.previous_TB_indices,
                                                                                  l.previous_TB_indices);
                                  else
                                      return r.previous_TB_indices.size() < l.previous_TB_indices.size();
                              });
            // GCOVR_EXCL_STOP

            if (config.output_verbose_statistics)
            {
                size_t dummy{};
                chopper::layout::hibf_statistics global_stats{config, sketches, cardinalities};
                global_stats.hibf_layout = hibf_layout;
                global_stats.print_header_to(std::cout);
                global_stats.print_summary_to(dummy, std::cout);
            }
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);
        hibf_layout.write_to(fout);
    }
    else
    {
        std::vector<std::vector<size_t>> positions(config.hibf_config.tmax); // asign positions for each partition

        std::vector<size_t> ps;
        ps.resize(cardinalities.size());
        std::iota(ps.begin(), ps.end(), 0);
        partition_user_bins(config, ps, cardinalities, sketches, minHash_sketches, positions);

        std::vector<seqan::hibf::layout::layout> hibf_layouts(config.hibf_config.tmax); // multiple layouts

#pragma omp parallel for schedule(dynamic) num_threads(config.hibf_config.threads)
        for (size_t i = 0; i < config.hibf_config.tmax; ++i)
        {
            // reset tmax to fit number of user bins in layout
            auto local_hibf_config = config.hibf_config; // every thread needs to set individual tmax
            local_hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(positions[i].size()))));

            config.dp_algorithm_timer.start();
            hibf_layouts[i] = seqan::hibf::layout::compute_layout(local_hibf_config,
                                                                  cardinalities,
                                                                  sketches,
                                                                  std::move(positions[i]),
                                                                  config.union_estimation_timer,
                                                                  config.rearrangement_timer);
            config.dp_algorithm_timer.stop();
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);

        for (size_t i = 0; i < config.hibf_config.tmax; ++i)
            hibf_layouts[i].write_to(fout);
    }

    return 0;
}

} // namespace chopper::layout
