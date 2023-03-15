#pragma once

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>

namespace chopper::layout
{

//!\brief Depending on cli flags given, use HyperLogLog estimates and/or rearrangement algorithms
inline void arrange_user_bins(data_store & data, configuration const & config)
{
    if (!data.user_bins_arranged)
    {
        data.sketch_toolbox = sketch::user_bin_sequence{data.filenames, data.kmer_counts, data.sketches};
        data.sketch_toolbox.sort_by_cardinalities();

        if (!config.disable_estimate_union)
        {
            if (!config.disable_rearrangement)
                data.sketch_toolbox.rearrange_bins(config.max_rearrangement_ratio, config.threads);
        }

        data.user_bins_arranged = true;
    }
}

} // namespace chopper::layout
