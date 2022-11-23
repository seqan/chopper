#pragma once

#include <cassert>

#include <robin_hood.h>

#include <chopper/configuration.hpp>

namespace chopper::count
{

//!\brief Checks the `filename_clusters` for consistent files, either precomputed or sequence files.
inline void check_filenames(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                            configuration & config)
{
    assert(!filename_clusters.empty());
    assert(!filename_clusters.begin()->second.empty());

    auto case_insensitive_string_ends_with = [] (std::string_view str, std::string_view suffix)
    {
        size_t const suffix_length{suffix.size()};
        size_t const str_length{str.size()};

        if (suffix_length > str_length)
            return false;

        for (size_t j = 0, s_start = str_length - suffix_length; j < suffix_length; ++j)
            if (std::tolower(str[s_start + j]) != std::tolower(suffix[j]))
                return false;

        return true;
    };

    // If the first filename ends in .minimizer we expect all files to end in .minimizer
    config.precomputed_files = case_insensitive_string_ends_with(filename_clusters.begin()->second[0], ".minimizer");

    for (auto const & [key, filenames] : filename_clusters)
    {
        for (auto const & filename : filenames)
        {
            if (config.precomputed_files && !case_insensitive_string_ends_with(filename, ".minimizer"))
            {
                throw std::invalid_argument{"You are providing precomputed files but the file " + filename +
                                            " does not have the correct file extension (.minimizer)."
                                            " Mixing non-/precomputed files is not allowed."};
            }
            else if (!config.precomputed_files && case_insensitive_string_ends_with(filename, ".minimizer"))
            {
                throw std::invalid_argument{"You are providing sequence files but the file " + filename +
                                            " was identified as a precomputed file (.minimizer)."
                                            " Mixing non-/precomputed files is not allowed."};
            }
        }
    }
}

} // namespace chopper::count
