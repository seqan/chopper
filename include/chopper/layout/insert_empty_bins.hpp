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
        // insert in the back of the list. or kmer_counts[idx] - kmer_counts[idx+1] to interpolate.
        empty_bins.insert(empty_bins.begin() + insertion_idx, user_bin_kmer_counts.at(insertion_idx));

        if (hll)
        {
            chopper::sketch::hyperloglog empty_sketch(sketch_bits);
            sketches.insert(sketches.begin() + insertion_idx, empty_sketch); //Insert an empty sketch.
        }
    }
}

} // namespace chopper::layout
