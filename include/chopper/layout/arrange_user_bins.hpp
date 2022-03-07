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
        bin_sequence.sort_by_cardinalities();

        if (config.estimate_union)
        {
            bin_sequence.read_hll_files(config.sketch_directory);
            if (config.rearrange_user_bins)
                bin_sequence.rearrange_bins(config.max_rearrangement_ratio, config.threads);

            bin_sequence.precompute_interval_union_estimations(data.union_estimates, config.threads);
        }

        data.user_bins_arranged = true;
    }
}

} // namespace chopper::layout
