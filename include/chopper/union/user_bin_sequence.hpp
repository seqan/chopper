#pragma once

#include <seqan3/std/filesystem>
#include <seqan3/utility/views/to.hpp>
#include <seqan3/std/ranges>

#include <fstream>
#include <queue>
#include <random>

#include <omp.h>

#include <robin_hood.h>

#include <chopper/union/hyperloglog.hpp>

class user_bin_sequence
{
private:
    //!\brief type for a node in the clustering tree when for the rearrangement
    struct clustering_node
    {
        // children in the tree
        size_t left;
        size_t right;
        // hll sketch of the union if the node is still a root
        hyperloglog hll;
    };

    //!\brief element of the second priority queue layer of the distance matrix
    struct neighbor
    {
        size_t id;
        double dist;

        bool operator>(neighbor const & other) const
        {
            return dist > other.dist;
        }
    };

    //!\brief type of a min heap based priority queue
    using prio_queue = std::priority_queue<neighbor, std::vector<neighbor>, std::greater<neighbor>>;

    //!\brief entry of the distance matrix that has the id of a cluster with its neighbors in a prio queue
    struct entry
    {
        size_t id;
        prio_queue pq;
    };

    //!\brief type of the distance matrix for the clustering for the rearrangement
    using distance_matrix = std::vector<entry>;

    //!\brief A reference to the filenames of the user input sequences.
    std::vector<std::string> & filenames;

    //!\brief A referece to kmer counts associated with the above files used to pack user bin into technical bins.
    std::vector<size_t> & user_bin_kmer_counts;

    //!\brief HyperLogLog sketches on the k-mer sets of the sequences from the files of filenames.
    std::vector<hyperloglog> sketches;

public:
    user_bin_sequence() = delete; //!< Deleted.
    user_bin_sequence(user_bin_sequence const &) = default; //!< Defaulted.
    user_bin_sequence & operator=(user_bin_sequence const &) = default; //!< Defaulted.
    user_bin_sequence(user_bin_sequence &&) = default; //!< Defaulted.
    user_bin_sequence & operator=(user_bin_sequence &&) = default; //!< Defaulted.
    ~user_bin_sequence() = default; //!< Defaulted.

    /*!\brief A sequence of user bins for which filenames and counts are given.
     * \param[in] filenames_ filenames of the sequence files for the user bins
     * \param[in] user_bin_kmer_counts_ counts of the k-mer sets of the bins corresponding to filenames
     */
    user_bin_sequence(std::vector<std::string> & filenames_,
                      std::vector<size_t> & user_bin_kmer_counts_) :
        filenames{filenames_},
        user_bin_kmer_counts{user_bin_kmer_counts_}
    {}

    //!\brief Sorts filenames and cardinalities by looking only at the cardinalities.
    void sort_by_cardinalities()
    {
        // generate permutation of indices sorted in descinding order by cardinalities
        auto permutation = std::views::iota(0ul, user_bin_kmer_counts.size()) | seqan3::views::to<std::vector>;
        assert(permutation.size() == user_bin_kmer_counts.size());
        auto cardinality_compare = [this] (auto const i1, auto const i2)
                                        { return user_bin_kmer_counts[i2] < user_bin_kmer_counts[i1]; };
        std::sort(permutation.begin(), permutation.end(), cardinality_compare);

        apply_permutation(permutation);
    }

    /*!\brief Restore the HLL sketches from the files in hll_dir
    * \param[in] hll_dir path to the directory where hll caches will be found 
    */
    void read_hll_files(std::filesystem::path const & hll_dir)
    {
        if (hll_dir.empty())
        {
            throw std::runtime_error("A directory where the HyperLogLog sketches are stored must be given "
                                     "when union estimates are enabled");
        }

        sketches.reserve(filenames.size());

        try
        {
            for (auto & filename : filenames)
            {
                std::filesystem::path path = hll_dir / std::filesystem::path(filename).stem();
                path += ".hll";
                std::ifstream hll_fin(path, std::ios::binary);

                // the sketch bits will be automatically read from the files
                sketches.emplace_back().restore(hll_fin);
            }
        }
        catch (std::runtime_error const & err)
        {
            std::cerr << "[CHOPPER PACK ERROR] Something went wrong trying to read the HyperLogLog sketches from files:\n"
                      << err.what() << '\n';
            exit(1);
        }
    }

    /*!\brief For all intervals of filenames: estimate the cardinality of the union
     * of k-mer sets of all sequences in the files of the interval.
     * estimates[i][j] will be the union cardinality estimate of the interval j, ..., i.
     * This unintuitive convention is chosen for cache efficiency in the hierarchical binning.
     * \param[in] num_threads_ the number of threads to use
     * \param[out] estimates output table
     */
    void estimate_interval_unions(std::vector<std::vector<uint64_t>> & estimates, size_t const num_threads_) const
    {
        estimates.clear();
        size_t const n = filenames.size();
        estimates.resize(n);

        size_t const chunk_size_ = std::floor(std::sqrt(n));

        #pragma omp parallel num_threads(num_threads_)
        {
            // initialize estimates
            #pragma omp for
            for (size_t i = 0; i < n; ++i)
            {
                estimates[i].resize(i + 1);
            }

            // fill estimates
            #pragma omp for schedule(nonmonotonic: dynamic, chunk_size_)
            for (size_t i = 0; i < n; ++i)
            {
                hyperloglog temp_hll = sketches[i];
                estimates[i][i] = user_bin_kmer_counts[i];

                for (size_t j = i + 1; j < n; ++j)
                {
                    estimates[j][i] = static_cast<uint64_t>(temp_hll.merge_and_estimate_SIMD(sketches[j]));
                }
            }
        }
    }

    /*!\brief Rearrange filenames, sketches and counts such that similar bins are close to each other
     * \param[in] max_ratio the maximal cardinality ratio in the clustering intervals (must be <= 1 and >= 0)
     * \param[in] num_threads the number of threads to use
     */
    void rearrange_bins(double const max_ratio, size_t const num_threads)
    {
        std::vector<size_t> permutation;

        size_t first = 0;
        size_t last = 1;

        while (first < filenames.size())
        {
            // size difference is too large or sequence is over -> do the clustering
            if (last == filenames.size() || user_bin_kmer_counts[first] * max_ratio > user_bin_kmer_counts[last])
            {
                // if this is not the first group, we want one bin overlap
                cluster_bins(permutation, first, last, num_threads);
                first = last;
            }
            ++last;
        }

        apply_permutation(permutation);
    }

private:
    /*!\brief Perform an agglomerative clustering variant on the index range [first:last)
     * \param[in] first id of the first cluster of the interval
     * \param[in] last id of the last cluster of the interval plus one
     * \param[in] num_threads_ the number of threads to use
     * \param[out] permutation append the new order to this
     */
    void cluster_bins(std::vector<size_t> & permutation,
                      size_t const first,
                      size_t const last,
                      size_t const num_threads_)
    {
        assert(num_threads_ >= 1);

        size_t const n = filenames.size();
        size_t const chunk_size_ = std::floor(std::sqrt(n));

        size_t const prune_steps = chunk_size_;
        size_t steps_without_prune = 0;

        size_t const none = std::numeric_limits<size_t>::max();
        /* internal map that stores the distances
        *
        * The first layer is a hash map with the ids of active clusters as keys.
        * The values (second layer) are priority queues with neighbors of the cluster
        * with the respective key in the first layer.
        * These neighbors are themselves clusters with an id and store a distance to the
        * cluster of the first layer.
        */
        distance_matrix dist;
        dist.reserve(n + 1);

        // map that indicates which ids of clusters are still in the distance matrix
        // the values are the indices where the priority queue for the given id as key can be found in dist
        robin_hood::unordered_flat_map<size_t, size_t> remaining_ids;

        // clustering tree stored implicitly in a vector
        std::vector<clustering_node> clustering;
        clustering.reserve(2 * n);

        // cache for hll cardinality estimates
        std::vector<double> estimates;
        estimates.reserve(2 * n);

        // every thread will write its observed id with minimal distance to some other here
        // id == none means that the thread observed only empty or no priority queues
        std::vector<size_t> min_ids(num_threads_, none);

        // these will be the new ids for new clusters
        // the first one is invalid, but it will be incremented before it is used for the first time
        size_t new_id = last - 1;

        // initialize clustering and estimates
        for (size_t id = first; id < last; ++id)
        {
            // id i is at the index i - first
            clustering.push_back({none, none, sketches[id]});
            estimates.emplace_back(sketches[id].estimate());
        }

        // if this is not the first group, we want to have one overlapping bin
        size_t previous_rightmost = none;
        if (first != 0)
        {
            // For all other clusters, their id is also their index in filesnames, sketches etc. .
            // This is important, because their id is then inserted into the clustering.
            // This does not work for previous rightmost, because its index does not necessarily lie on
            // the continuous spectrum from first to last. We run into a problem, because the entries are
            // stored in vectors. Therefore we give previous_rightmost a different id (==last). This is
            // fine, because we only need the HLL sketch of the actual index. previous_rightmost will be ignored
            // in the traceback anyway and won't be added to the permutation in this step.
            size_t actual_previous_rightmost = permutation.back();
            ++new_id;
            previous_rightmost = new_id;

            clustering.push_back({none, none, sketches[actual_previous_rightmost]});
            estimates.emplace_back(sketches[actual_previous_rightmost].estimate());
        }

        // initialize priority queues in the distance matrix (sequentially)
        for (size_t id = first; id < first + clustering.size(); ++id)
        {
            // empty priority queue for every item in clustering
            dist.push_back({id, prio_queue{}});
            remaining_ids[id] = id - first;
        }

        #pragma omp parallel num_threads(num_threads_)
        {
            double min_dist = std::numeric_limits<double>::max();
            // minimum distance exclusively for this thread

            // initialize all the priority queues of the distance matrix
            // while doing that, compute the first min_id
            #pragma omp for schedule(nonmonotonic: dynamic, chunk_size_)
            for (size_t i = 0; i < clustering.size(); ++i)
            {
                for (size_t j = 0; j < clustering.size(); ++j)
                {
                    // we only want one diagonal of the distance matrix
                    if (i < j)
                    {
                        // this must be a copy, because merging changes the hll sketch
                        hyperloglog temp_hll = clustering[i].hll;
                        double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering[j].hll);
                        // Jaccard distance estimate
                        double const distance = 2 - (estimates[i] + estimates[j]) / estimate_ij;
                        dist[i].pq.push({j + first, distance});
                    }
                }
                if (dist[i].pq.empty()) continue;

                // check if the just initialized priority queue contains the minimum value for this thread
                neighbor const & curr = dist[i].pq.top();
                if (curr.dist < min_dist)
                {
                    min_dist = curr.dist;
                    min_ids[omp_get_thread_num()] = dist[i].id;
                }
            } // implicit barrier

            // a single thread shuffles dist to approximately balance loads in static scheduling
            #pragma omp single
            random_shuffle(dist, remaining_ids);

            // main loop of the clustering
            // keep merging nodes until we have a complete tree
            while (remaining_ids.size() > 1)
            {
                // Wait for all threads to have evaluated remaining_ids.size() as remaining_ids
                // may be modified by the following pragma omp single.
                #pragma omp barrier

                #pragma omp single
                {
                    // perform critical update
                    // increment id for the new cluster (must be done at the beginning)
                    ++new_id;

                    // compute the final min_id from the min_ids of the worker threads
                    size_t min_id = min_ids[0];
                    double min_dist = std::numeric_limits<double>::max();
                    for (auto candidate_id : min_ids)
                    {
                        // check if the thread saw any id
                        if (candidate_id == none) continue;

                        size_t const dist_index = remaining_ids.at(candidate_id);
                        neighbor const & curr = dist[dist_index].pq.top();
                        if (curr.dist < min_dist)
                        {
                            min_dist = curr.dist;
                            min_id = candidate_id;
                        }
                    }

                    size_t const min_index = remaining_ids.at(min_id); // how can min_id be none?
                    size_t const neighbor_id = dist[min_index].pq.top().id;

                    // merge the two nodes with minimal distance together insert the new node into the clustering
                    clustering.push_back({min_id, neighbor_id, std::move(clustering[min_id - first].hll)});
                    estimates.emplace_back(clustering.back().hll
                                                .merge_and_estimate_SIMD(clustering[neighbor_id - first].hll));

                    // remove old ids
                    remaining_ids.erase(min_id);
                    remaining_ids.erase(neighbor_id);

                    // overwrite one of the removed entries with the new one
                    remaining_ids[new_id] = min_index;
                    dist[min_index] = {new_id, prio_queue{}};

                    // prune the distance matrix to reduce overhead due to inactive entries
                    ++steps_without_prune;
                    if (steps_without_prune > prune_steps)
                    {
                        prune(dist, remaining_ids);
                        steps_without_prune = 0;
                    }
                } // implicit barrier

                // reset values for the computation of the new minimum
                min_ids[omp_get_thread_num()] = none;
                min_dist = std::numeric_limits<double>::max();

                hyperloglog const new_hll = clustering.back().hll;

                // update distances in dist
                // while doing that, compute the new min_id
                #pragma omp for schedule(static)
                for (size_t i = 0; i < dist.size(); ++i)
                {
                    size_t other_id = dist[i].id;
                    if (other_id == new_id || !remaining_ids.contains(other_id)) continue;

                    // this must be a copy, because merge_and_estimate_SIMD() changes the hll
                    hyperloglog temp_hll = new_hll;
                    double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering[other_id - first].hll);
                    // Jaccard distance estimate
                    double const distance = 2 - (estimates[other_id - first] + estimates.back()) / estimate_ij;
                    dist[i].pq.push({new_id, distance});

                    // make sure the closest neighbor is not yet deleted (this is a lazy update)
                    while (!remaining_ids.contains(dist[i].pq.top().id))
                    {
                        dist[i].pq.pop();
                    }

                    // check if the just updated priority queue contains the minimum value for this thread
                    neighbor const & curr = dist[i].pq.top();
                    if (curr.dist < min_dist)
                    {
                        min_dist = curr.dist;
                        min_ids[omp_get_thread_num()] = other_id;
                    }
                } // implicit barrier
            }
        } // end of the parallel region

        size_t final_root_index = remaining_ids.begin()->second;
        size_t final_root_id = dist[final_root_index].id;

        // rotate the previous rightmost to the left so that it has the correct place in the permutation
        if (first != 0)
        {
            rotate(clustering, previous_rightmost, first, final_root_id);
        }

        // traceback into permutation and ignore the previous rightmost
        trace(clustering, permutation, previous_rightmost, first, final_root_id);
    }

    /*!\brief Randomly entries in dist while keeping track of the changes of indices
     * \param[in] dist the distance matrix (vector of priority queues) to shuffle
     * \param[in] remaining_ids the map with information about which ids remain at which index
     */
    void random_shuffle(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids)
    {
        size_t const n = dist.size();

        std::random_device rd;
        // random generator seeded with device's random bit source
        std::mt19937 gen(rd());

        for (size_t i = 0; i < n - 1; ++i)
        {
            std::uniform_int_distribution<size_t> distrib(i, n - 1);
            size_t const swap_i = distrib(gen);

            size_t const id = dist[i].id;
            size_t const swap_id = dist[swap_i].id;

            // swap entries and update the reming ids, because the indices in dist changed
            std::swap(dist[i], dist[swap_i]);
            std::swap(remaining_ids[id], remaining_ids[swap_id]);
        }
    }

    /*!\brief Delete inactive entries out of dist and shrink to fit its size while keeping track of the changes of indices
     * \param[in] dist the distance matrix (vector of priority queues) to prune
     * \param[in] remaining_ids the map with information about which ids remain at which index
     */
    void prune(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids)
    {
        if (dist.empty()) return;

        // index of the first entry after the valid range
        size_t valid_range_end = 0;
        // index of the first entry before the invalid range
        size_t invalid_range_start = dist.size() - 1;

        while (valid_range_end != invalid_range_start)
        {
            if (remaining_ids.contains(dist[valid_range_end].id))
            {
                ++valid_range_end;
            }
            else if (!remaining_ids.contains(dist[invalid_range_start].id))
            {
                --invalid_range_start;
            }
            else
            {
                // If we arrive here, then valid_range_end has an invalid id
                // and invalid_range_start has a valid id. The correspoding entries should be swapped
                std::swap(dist[valid_range_end], dist[invalid_range_start]);

                // update the index of the valid entry
                remaining_ids.at(dist[valid_range_end].id) = valid_range_end;
            }
        }

        // check the last element between the valid and invalid range
        if (remaining_ids.contains(dist[valid_range_end].id))
        {
            ++valid_range_end;
        }

        // cut off invalid values
        dist.resize(valid_range_end);
    }

    /*!\brief Rotate the previous rightmost bin to the left of the clustering tree
     * \param[in, out] clustering the tree to do the rotation on
     * \param[in] previous_rightmost the id of the node to be rotated to the left
     * \param[in] first the id of the first node in the interval to shift the index
     * \param[in] id the id of the current node
     *
     * \return whether previous rightmost was in the subtree rooted at id
     */
    bool rotate(std::vector<clustering_node> & clustering,
                size_t const previous_rightmost,
                size_t const first,
                size_t const id)
    {
        if (id == previous_rightmost) return true;

        clustering_node & curr = clustering[id - first];
        if (curr.left == std::numeric_limits<size_t>::max())
        {
            return false;
        }

        // nothing to do if previous_rightmost is in the left subtree
        if(rotate(clustering, previous_rightmost, first, curr.left)) return true;

        // rotate if previous_rightmost is in the right subtree
        if(rotate(clustering, previous_rightmost, first, curr.right))
        {
            std::swap(curr.left, curr.right);
            return true;
        }

        // previous_rightmost is not in this subtree
        return false;
    }

    /*!\brief Do a recursive traceback to find the order of leaves in the clustering tree
     * \param[in] clustering the tree to do the traceback on
     * \param[out] permutation append the new order to this
     * \param[in] previous_rightmost the id of the node on the left which should be ignored
     * \param[in] first the id of the first node in the interval to shift the index
     * \param[in] id the id of the current node
     */
    void trace(std::vector<clustering_node> const & clustering,
               std::vector<size_t> & permutation,
               size_t const previous_rightmost,
               size_t const first,
               size_t const id)
    {
        clustering_node const & curr = clustering[id - first];

        if (curr.left == std::numeric_limits<size_t>::max())
        {
            if (id != previous_rightmost) permutation.push_back(id);
            return;
        }

        trace(clustering, permutation, previous_rightmost, first, curr.left);
        trace(clustering, permutation, previous_rightmost, first, curr.right);
    }

    /*!\brief Apply a given permutation to filenames, user_bin_kmer_counts and sketches
     * \param[in] permutation the permutation to apply
     */
    void apply_permutation(std::vector<size_t> & permutation)
    {
        for (size_t i = 0; i < permutation.size(); ++i)
        {
            size_t current = i;
            while (i != permutation[current])
            {
                size_t next = permutation[current];
                std::swap(filenames[current], filenames[next]);
                std::swap(user_bin_kmer_counts[current], user_bin_kmer_counts[next]);
                if (!sketches.empty()) std::swap(sketches[current], sketches[next]);
                permutation[current] = current;
                current = next;
            }
            permutation[current] = current;
        }
    }
};
