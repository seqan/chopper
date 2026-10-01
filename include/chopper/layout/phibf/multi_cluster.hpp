// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::phibf::MultiCluster.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <optional>
#include <vector>

#include <chopper/layout/fast_layout_cluster.hpp>

namespace chopper::layout::phibf
{

/*!\brief A cluster of clusters of user bins, used by the `lsh_sim` partitioning of the partitioned HIBF.
 *
 * A MultiCluster is created from a Cluster (the result of a first, very similar, LSH clustering). The second LSH
 * clustering then joins MultiClusters, such that a MultiCluster holds several clusters of very similar user bins.
 *
 * The validity invariants are the same as for chopper::layout::Cluster, with `contained_user_bins()` being a vector of
 * clusters instead of a vector of user bins. Hence, chopper::layout::LSH_find_representative_cluster and
 * chopper::layout::LSH_fill_hashtable work on MultiClusters, too.
 */
struct MultiCluster
{
private:
    //!\brief The id of the cluster.
    size_t representative_id{};

    //!\brief The clusters of user bins contained in this cluster.
    std::vector<std::vector<size_t>> user_bins{};

    //!\brief The id of the cluster that this cluster's user bins were moved to, if they were moved.
    std::optional<size_t> moved_id{std::nullopt};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    MultiCluster() = default;                                 //!< Defaulted.
    MultiCluster(MultiCluster const &) = default;             //!< Defaulted.
    MultiCluster(MultiCluster &&) = default;                  //!< Defaulted.
    MultiCluster & operator=(MultiCluster const &) = default; //!< Defaulted.
    MultiCluster & operator=(MultiCluster &&) = default;      //!< Defaulted.
    ~MultiCluster() = default;                                //!< Defaulted.

    /*!\brief Creates a MultiCluster from a Cluster. Not explicit, such that `{Cluster{0}}` converts.
     * \param[in] clust The cluster. If it has been moved, the MultiCluster is moved to the same id and is empty.
     *                  Otherwise, the MultiCluster contains the user bins of `clust` as its single cluster.
     */
    MultiCluster(Cluster const & clust) : representative_id{clust.id()}
    {
        if (clust.has_been_moved())
            moved_id = clust.moved_to_cluster_id();
        else
            user_bins.push_back(clust.contained_user_bins());
    }
    //!\}

    //!\brief Returns the id of the cluster.
    size_t id() const
    {
        return representative_id;
    }

    //!\brief Returns the clusters of user bins contained in the cluster.
    std::vector<std::vector<size_t>> const & contained_user_bins() const
    {
        return user_bins;
    }

    //!\brief Returns whether the user bins of this cluster were moved to another cluster (see move_to).
    bool has_been_moved() const
    {
        return moved_id.has_value();
    }

    //!\brief Returns whether the cluster contains no clusters.
    bool empty() const
    {
        return user_bins.empty();
    }

    //!\brief Returns the number of clusters in the cluster.
    size_t size() const
    {
        return user_bins.size();
    }

    /*!\brief Checks that the cluster at position `id` is either valid or moved, see chopper::layout::Cluster.
     * \param[in] id The position of the cluster, i.e., the id it must have.
     */
    bool is_valid(size_t const id) const
    {
        bool const ids_equal = representative_id == id;
        bool const properly_moved = has_been_moved() && empty();
        bool const not_moved = !has_been_moved() && !empty();

        return ids_equal && (properly_moved || not_moved);
    }

    //!\brief Returns the id of the cluster that this cluster's user bins were moved to. The cluster must be moved.
    size_t moved_to_cluster_id() const
    {
        assert(moved_id.has_value());
        assert(is_valid(representative_id));
        return moved_id.value();
    }

    /*!\brief Moves all clusters of this cluster to `target_cluster` and marks this cluster as moved.
     * \param[in,out] target_cluster The cluster to move the clusters to. Must be a different cluster.
     */
    void move_to(MultiCluster & target_cluster)
    {
        target_cluster.user_bins.insert(target_cluster.user_bins.end(), this->user_bins.begin(), this->user_bins.end());
        this->user_bins.clear();
        moved_id = target_cluster.id();
    }

    /*!\brief Sorts the user bins within each cluster by descending cardinality, and the clusters by descending size.
     * \param[in] cardinalities The cardinality of each user bin, indexed by user bin.
     */
    void sort_by_cardinality(std::vector<size_t> const & cardinalities)
    {
        auto cmp = [&cardinalities](auto const & v1, auto const & v2)
        {
            return cardinalities[v2] < cardinalities[v1];
        };
        for (auto & user_bin_cluster : user_bins)
            std::sort(user_bin_cluster.begin(), user_bin_cluster.end(), cmp);

        auto cmp_clusters = [](auto const & c1, auto const & c2)
        {
            return c2.size() < c1.size();
        };
        std::sort(user_bins.begin(), user_bins.end(), cmp_clusters);
    }
};

} // namespace chopper::layout::phibf
