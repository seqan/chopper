#pragma once

#include <chopper/configuration.hpp>
#include <chopper/layout/data_store.hpp>

namespace chopper::layout
{

/*!\brief Return indices of where to insert empty bins in various datastructures based on k-mer counts.
 * \details sample evenly among sorted kmer counts, according to a certain percentage.
 * K-mer counts should already be sorted before calling this function.
 * \param[in] empty_bin_fraction Currently a maximum of 1 is supported.
 * \param[in] original_size Size of the datastructures before inserting.
 * \author Myrthe Willemsen
 */

std::vector<size_t> empty_bin_indices (double empty_bin_fraction, size_t original_size) {
    int stepsize = 1 / empty_bin_fraction; //this way, it works uptill 100%.
    assert(stepsize > 0);
    std::vector<size_t> insertion_indices;
    for (size_t idx = 0; idx < original_size; idx += stepsize) {
        insertion_indices.push_back(std::round(idx));
    }
    return insertion_indices;
}

/*!\brief Depending on cli flags given, use HyperLogLog estimates and/or rearrangement algorithms
 * \details If empty bins need to be inserted, then first sort the k-mer counts by cardinality and then insert empty bins.
 * \author Adapted by Myrthe Willemsen
 */
// TODO myrthe bin 07 disapears. in the halved example it is for instance _01.
inline void arrange_user_bins(data_store & data, configuration const & config)
{
    if (!data.user_bins_arranged)
    {
        data.sketch_toolbox = sketch::user_bin_sequence{data.filenames, data.kmer_counts};
        data.sketch_toolbox.sort_by_cardinalities(); // sort filenames by k-mer count
        if (config.estimate_union){data.sketch_toolbox.read_hll_files(config.sketch_directory);} // read hll files

        if (config.update_ubs != 0) { // insert empty bins
            std::vector<size_t> insertion_indices = empty_bin_indices(config.update_ubs, data.filenames.size());
            data.sketch_toolbox.insert_empty_bins(insertion_indices, config.estimate_union, config.sketch_bits);
            data.insert_empty_bins(insertion_indices); // filenames and kmer counts are already updated, by updating it in the sketch toolbox.
        }


        if (config.estimate_union)
        {
            if (config.rearrange_user_bins)data.sketch_toolbox.rearrange_bins(config.max_rearrangement_ratio, config.threads);

        }

        data.user_bins_arranged = true;
    }
}




} // namespace chopper::layout
