// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::adjust_seed.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>

namespace chopper
{

/*\brief foo
 */
struct Cluster
{
protected:
    size_t representative_id{}; // representative id of the cluster; identifier;

    std::vector<size_t> user_bins{}; // the user bins contained in thus cluster

    std::optional<size_t> moved_id{std::nullopt}; // where this Clusters user bins where moved to

public:
    Cluster() = default;
    Cluster(const Cluster&) = default;
    Cluster(Cluster&&) = default;
    Cluster& operator=(const Cluster&) = default;
    Cluster& operator=(Cluster&&) = default;
    ~Cluster() = default;

    Cluster(size_t const id)
    {
        representative_id = id;
        user_bins = {id};
    }

    Cluster(size_t const id, size_t const user_bins_id)
    {
        representative_id = id;
        user_bins = {user_bins_id};
    }

    size_t id() const
    {
        return representative_id;
    }

    std::vector<size_t> const & contained_user_bins() const
    {
        return user_bins;
    }

    bool has_been_moved() const
    {
        return moved_id.has_value();
    }

    bool empty() const
    {
        return user_bins.empty();
    }

    size_t size() const
    {
        return user_bins.size();
    }

    bool is_valid(size_t id) const
    {
        bool const ids_equal = representative_id == id;
        bool const properly_moved = has_been_moved() && empty() && moved_id.has_value();
        bool const not_moved = !has_been_moved() && !empty();

        return ids_equal && (properly_moved || not_moved);
    }

    size_t moved_to_cluster_id() const
    {
        assert(moved_id.has_value());
        assert(is_valid(representative_id));
        return moved_id.value();
    }

    void move_to(Cluster & target_cluster)
    {
        target_cluster.user_bins.insert(target_cluster.user_bins.end(), this->user_bins.begin(), this->user_bins.end());
        this->user_bins.clear();
        moved_id = target_cluster.id();
    }

    void sort_by_cardinality(std::vector<size_t> const & cardinalities)
    {
        std::sort(user_bins.begin(), user_bins.end(), [&cardinalities](auto const & v1, auto const & v2){ return cardinalities[v2] < cardinalities[v1]; });
    }
};

struct MultiCluster : Cluster
{
protected:
    std::vector<std::vector<size_t>> user_bins{}; // the user bins contained in this cluster

public:
    MultiCluster() = default;
    MultiCluster(const MultiCluster&) = default;
    MultiCluster(MultiCluster&&) = default;
    MultiCluster& operator=(const MultiCluster&) = default;
    MultiCluster& operator=(MultiCluster&&) = default;
    ~MultiCluster() = default;

    MultiCluster(Cluster const & clust)
    {
        representative_id = clust.id();
        if (clust.has_been_moved())
        {
            moved_id = clust.moved_to_cluster_id();
        }
        else
        {
            user_bins.push_back(clust.contained_user_bins());
        }
    }

    std::vector<std::vector<size_t>> const & contained_user_bins() const
    {
        return user_bins;
    }

    // needs to be defined again because of redefinition of `user_bins`
    bool empty() const
    {
        return user_bins.empty();
    }

    // needs to be defined again because of redefinition of `user_bins`
    size_t size() const
    {
        return user_bins.size();
    }

    bool is_valid(size_t id) const
    {
        bool const ids_equal = representative_id == id;
        bool const properly_moved = has_been_moved() && empty() && moved_id.has_value();
        bool const not_moved = !has_been_moved() && !empty();

        return ids_equal && (properly_moved || not_moved);
    }

    void move_to(MultiCluster & target_cluster)
    {
        target_cluster.user_bins.insert(target_cluster.user_bins.end(), this->user_bins.begin(), this->user_bins.end());
        this->user_bins.clear();
        moved_id = target_cluster.id();
    }

    // sort user bins within a Cluster by cardinality and the clusters themselves by size
    void sort_by_cardinality(std::vector<size_t> const & cardinalities)
    {
        auto cmp = [&cardinalities](auto const & v1, auto const & v2){ return cardinalities[v2] < cardinalities[v1]; };
        for (auto & user_bin_cluster : user_bins)
            std::sort(user_bin_cluster.begin(), user_bin_cluster.end(), cmp);

        auto cmp_clusters = [](auto const & c1, auto const & c2){ return c2.size() < c1.size(); };
        std::sort(user_bins.begin(), user_bins.end(), cmp_clusters);
    }
};


// A valid cluster is one that hasn't been moved but actually contains user bins
// A valid cluster at position i is identified by the following equality: cluster[i].size() >= 1 && cluster[i][0] == i
// A moved cluster is one that has been joined and thereby moved to another cluster
// A moved cluster i is identified by the following: cluster[i].size() == 1 && cluster[i][0] != i
// returns position of the representative cluster
template <typename cluster_type>
size_t LSH_find_representative_cluster(std::vector<cluster_type> const & clusters, size_t current_id)
{
    std::reference_wrapper<cluster_type const> representative = clusters[current_id];

    assert(representative.get().is_valid(current_id));

    while (representative.get().has_been_moved())
    {
        current_id = representative.get().moved_to_cluster_id();
        representative = clusters[current_id]; // replace by next cluster
        assert(representative.get().is_valid(current_id));
    }

    return current_id;
}

} // namespace chopper
