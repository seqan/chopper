#pragma once

#include <chopper/layout/configuration.hpp>
#include <chopper/layout/data_store.hpp>
#include <chopper/sketch/user_bin_sequence.hpp>

namespace chopper::layout
{

//!\brief Depending on cli flags given, use HyperLogLog estimates and/or rearrangement algorithms
inline void arrange_user_bins(data_store & data, configuration const & config)
{
    if (!data.user_bins_arranged)
    {
        chopper::sketch::user_bin_sequence bin_sequence{data.filenames, data.kmer_counts};
        bin_sequence.sort_by_cardinalities(); // myrthe 1) sort filenames sort_by_cardinalities

        if (config.estimate_union){bin_sequence.read_hll_files(config.hll_dir);} // myrthe 2) read hll filenames here

        if (config.update_ubs != 0) { // myrthe 3) insert empty bins
            bin_sequence.insert_empty_bins(config.update_ubs, config.estimate_union);
        }

        if (config.estimate_union)
        {
            if (config.rearrange_bins)
                bin_sequence.rearrange_bins(config.max_ratio, config.num_threads);

            bin_sequence.estimate_interval_unions(data.union_estimates, config.num_threads);
        }

        data.user_bins_arranged = true;
    }
}

} // namespace chopper::layout
