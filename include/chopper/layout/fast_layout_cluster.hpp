// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <optional>
#include <vector>

namespace chopper::layout
{

/*!\brief A cluster of user bins, used by the LSH clustering of the fast layout.
 *
 * In a vector of clusters in which the cluster at position `i` was created with id `i`, the cluster at position `i`
 * keeps `id() == i` and is either
 * - **valid**: it has not been moved and contains at least one user bin, or
 * - **moved**: its user bins were moved into another cluster (move_to), it is empty, and moved_to_cluster_id() is the
 *   id of that cluster. That cluster may itself have been moved, so moves can form a chain.
 *   LSH_find_representative_cluster follows the chain to the valid cluster.
 *
 * Note that `is_valid(i)` is true in both states: it checks that `id() == i` and that the cluster is in one of them.
 */
struct Cluster
{
private:
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

    explicit Cluster(size_t const id) : Cluster{id, id}
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

    void add_user_bin(size_t const user_bin)
    {
        user_bins.push_back(user_bin);
    }

    bool is_valid(size_t const id) const
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
        auto & target = target_cluster.user_bins;
        auto & source = this->user_bins;
#ifdef __cpp_lib_containers_ranges
        target.append_range(source);
#else
        target.insert(target.end(), source.cbegin(), source.cend());
#endif
        source = std::vector<size_t>{}; // .clear() AND release memory

        moved_id = target_cluster.id();
    }

    void sort_by_cardinality(std::vector<size_t> const & cardinalities)
    {
        std::ranges::sort(user_bins,
                          std::ranges::greater{},
                          [&cardinalities](size_t const i)
                          {
                              return cardinalities[i];
                          });
    }
};

// Follows the chain of moves, starting at clusters[current_id], and returns the position of the representative
// cluster, i.e., the valid cluster that holds the user bins now. See Cluster for valid and moved clusters.
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
