// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cassert>
#include <optional>
#include <vector>

namespace chopper::layout
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
    Cluster(Cluster const &) = default;
    Cluster(Cluster &&) = default;
    Cluster & operator=(Cluster const &) = default;
    Cluster & operator=(Cluster &&) = default;
    ~Cluster() = default;

    Cluster(size_t const id, size_t const user_bins_id) : representative_id{id}, user_bins({user_bins_id})
    {}

    Cluster(size_t const id) : Cluster{id, id}
    {}

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

    size_t pop_back()
    {
        size_t last = user_bins.back();
        user_bins.pop_back();
        return last;
    }

    void add_user_bin(size_t user_bin)
    {
        user_bins.push_back(user_bin);
    }

    bool is_valid(size_t id) const
    {
        bool const ids_equal = representative_id == id;
        bool const properly_moved = has_been_moved() && empty();
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
        std::ranges::sort(user_bins,
                          [&cardinalities](auto const & v1, auto const & v2)
                          {
                              return cardinalities[v1] > cardinalities[v2];
                          });
    }
};

// A valid cluster is one that hasn't been moved but actually contains user bins
// A valid cluster at position i is identified by the following equality: cluster[i].size() >= 1 && cluster[i][0] == i
// A moved cluster is one that has been joined and thereby moved to another cluster
// A moved cluster i is identified by the following: cluster[i].size() == 1 && cluster[i][0] != i
// returns position of the representative cluster
inline size_t LSH_find_representative_cluster(std::vector<Cluster> const & clusters, size_t current_id)
{
    std::reference_wrapper<Cluster const> representative = clusters[current_id];

    assert(representative.get().is_valid(current_id));

    while (representative.get().has_been_moved())
    {
        current_id = representative.get().moved_to_cluster_id();
        representative = clusters[current_id]; // replace by next cluster
        assert(representative.get().is_valid(current_id));
    }

    return current_id;
}

} // namespace chopper::layout