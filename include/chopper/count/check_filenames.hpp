#pragma once

#include <robin_hood.h>

#include <chopper/count/configuration.hpp>

namespace chopper::count
{

//!\brief Checks the `filename_clusters` for consistent files, wither precomputed or sequence files.
inline void check_filenames(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                            configuration & config)
{
    assert(!filename_clusters.empty());
    assert(!filename_clusters.begin()->second.empty());

    config.precomputed_files = std::filesystem::path{filename_clusters.begin()->second[0]}.extension() == ".minimizer";

    for (auto const & [key, filenames] : filename_clusters)
    {
        for (auto const & filename : filenames)
        {
            std::filesystem::path path{filename};

            if (config.precomputed_files && path.extension() != ".minimizer")
            {
                throw std::invalid_argument{"You are providing precomputed files but the file " + path.string() +
                                            " does not have the correct file extension (.minimizer)."
                                            " Mixing non-/precomputed files is not allowed."};
            }
            else if (!config.precomputed_files && path.extension() == ".minimizer")
            {
                throw std::invalid_argument{"You are providing sequence files but the file " + path.string() +
                                            " was identified as a precomputed file (.minimizer)."
                                            " Mixing non-/precomputed files is not allowed."};
            }
        }
    }
}

} // namespace chopper::count
