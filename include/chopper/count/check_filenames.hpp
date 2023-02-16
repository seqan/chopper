#pragma once

#include <cassert>
#include <robin_hood.h>

#include <chopper/configuration.hpp>

namespace chopper::count
{

//!\brief Checks the `filenames` for consistent files, either precomputed or sequence files.
inline void check_filenames(std::vector<std::string> & filenames, configuration & config)
{
    assert(!filenames.empty());

    auto case_insensitive_string_ends_with = [](std::string_view str, std::string_view suffix)
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

    // If the first filename ends in .minimiser we expect all files to end in .minimiser
    config.precomputed_files = case_insensitive_string_ends_with(filenames[0], ".minimiser");

    for (auto const & filename : filenames)
    {
        if (config.precomputed_files && !case_insensitive_string_ends_with(filename, ".minimiser"))
        {
            throw std::invalid_argument{"You are providing precomputed files but the file " + filename
                                        + " does not have the correct file extension (.minimiser)."
                                            " Mixing non-/precomputed files is not allowed."};
        }
        else if (!config.precomputed_files && case_insensitive_string_ends_with(filename, ".minimiser"))
        {
            throw std::invalid_argument{"You are providing sequence files but the file " + filename
                                        + " was identified as a precomputed file (.minimiser)."
                                            " Mixing non-/precomputed files is not allowed."};
        }
    }
}

} // namespace chopper::count
