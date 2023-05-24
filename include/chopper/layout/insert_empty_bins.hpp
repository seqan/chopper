#pragma once

#include <cassert>
#include <cmath>

#include <chopper/configuration.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::layout
{
/*!\brief insert empty bins in various datastructures based on k-mer counts.
    * \param[in] empty_bin_fraction Currently a maximum of 1 is supported.
    * \param[in] hll Whether HLL sketches are used in the layout algorithm
    * \author Myrthe Willemsen
    */
inline void insert_empty_bins(std::vector<size_t> const & insertion_indices,
                       bool hll,
                       uint8_t sketch_bits,
                       std::vector<std::string> & filenames,
                       std::vector<sketch::hyperloglog> & sketches,
                       std::vector<size_t> & user_bin_kmer_counts,
                       std::vector<bool> & empty_bins)
{
    for (size_t idx = 0; idx < insertion_indices.size(); ++idx)
    {
        // add `idx` because the indices will be shifted because the array grows longer upon inserting.
        size_t insertion_idx = insertion_indices[idx] + idx;
        filenames.insert(filenames.begin() + insertion_idx,
                          std::to_string(user_bin_kmer_counts.at(insertion_idx)) + ".empty_bin"); // +size of UB?

        // insert in the back of the list. or kmer_counts[idx] - kmer_counts[idx+1] to interpolate.
        user_bin_kmer_counts.insert(user_bin_kmer_counts.begin() + insertion_idx, user_bin_kmer_counts.at(insertion_idx));
        empty_bins.insert(empty_bins.begin() + insertion_idx, user_bin_kmer_counts.at(insertion_idx));

        if (hll)
        {
            chopper::sketch::hyperloglog empty_sketch(sketch_bits);
            sketches.insert(sketches.begin() + insertion_idx, empty_sketch); //Insert an empty sketch.
        }
    }
}
    /*!\brief Return indices of where to insert empty bins in various datastructures based on k-mer counts.
     * \details sample evenly among sorted kmer counts, according to a certain percentage.
     * K-mer counts should already be sorted before calling this function.
     * \param[in] empty_bin_fraction Currently a maximum of 1 is supported.
     * \param[in] original_size Size of the datastructures before inserting.
     * \param[in] order an array of the descending ordering of the k-mer counts
     * \author Myrthe Willemsen
     */
    inline std::vector<size_t> empty_bin_indices(double empty_bin_fraction,
                                                 size_t original_size,
                                                 std::vector<size_t> & order) {
        int stepsize = 1 / empty_bin_fraction; //this way, it works uptill 100%.
        assert(stepsize > 0);
        std::vector<size_t> insertion_indices;
        for (size_t idx = 0; idx < original_size; idx += stepsize) {
            insertion_indices.push_back(order[std::round(idx)]); // map the indices to their indices in the ordered version of the k-mer counts.
        }
        std::sort(insertion_indices.begin(), insertion_indices.end()); // sort the indices such that later the smallest indices are inserted first.
        return insertion_indices;
    }

    inline void insert_empty_bins(std::vector<bool> &empty_bins,
                                   std::vector<size_t> &empty_bin_cum_sizes,
                                   std::vector<size_t> &kmer_counts,
                                   std::vector<chopper::sketch::hyperloglog> &sketches,
                                   std::vector<std::string> &filenames,
                                   chopper::configuration &config) {
        empty_bins.resize(sketches.size());
        empty_bin_cum_sizes.resize(sketches.size());
        if (config.update_ubs != 0) // insert empty bins
        {
            std::vector<size_t> order(kmer_counts.size()); //order of the k-mer counts if they were sorted.
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return kmer_counts[i] > kmer_counts[j];
            }); // make sure that the empty bins are inserted based on the sorted k-mer counts.
            std::vector<size_t> insertion_indices = empty_bin_indices(config.update_ubs, sketches.size(), order);

            chopper::layout::insert_empty_bins(insertion_indices, !config.disable_estimate_union, config.sketch_bits,
                                               filenames, sketches, kmer_counts, empty_bins);

            empty_bins.resize(sketches.size());
            empty_bin_cum_sizes.resize(sketches.size());

            // create empty_bin_cum_sizes with cumulative sizes // This is currently not used
            if (empty_bins[0])
                empty_bin_cum_sizes[0] = kmer_counts.at(0);

            for (size_t j = 1; j < sketches.size(); ++j) {
                if (empty_bins[j])
                    empty_bin_cum_sizes[j] = kmer_counts.at(j);     //but this should also be done according to the positions array right? The empty bin cum sizes only works if it has the cumulative sizes are calculated according to the order of the order of the (re)arranged bins. Hence, it must be created by looping over the positions array.
                empty_bin_cum_sizes[j] += empty_bin_cum_sizes[j - 1];
            }

        }
    }

} // namespace chopper::layout
