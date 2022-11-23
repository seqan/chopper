#pragma once

#include <chopper/configuration.hpp>
#include <chopper/layout/data_store.hpp>

namespace chopper::layout
{

//!\brief Depending on cli flags given, use HyperLogLog estimates and/or rearrangement algorithms
inline void arrange_user_bins(data_store & data, configuration const & config)
{
    if (!data.user_bins_arranged)
    {
        data.sketch_toolbox = sketch::user_bin_sequence{data.filenames, data.kmer_counts};
        data.sketch_toolbox.sort_by_cardinalities();

        if (config.estimate_union)
        {
            data.sketch_toolbox.read_hll_files(config.sketch_directory);
            if (config.rearrange_user_bins)
                data.sketch_toolbox.rearrange_bins(config.max_rearrangement_ratio, config.threads);
        }

        data.user_bins_arranged = true;
    }
}

} // namespace chopper::layout
