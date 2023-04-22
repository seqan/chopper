#pragma once

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>

#include <chopper/configuration.hpp>
#include <chopper/count/check_filenames.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/prefixes.hpp>
#include <chopper/print_peak_memory_usage.hpp>

namespace chopper::count
{

inline int execute(configuration & config)
{
    auto filename_clusters = read_data_file(config);

    if (filename_clusters.empty())
        throw sharg::parser_error{
            sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

    chopper::count::check_filenames(filename_clusters, config);

    count_kmers(filename_clusters, config);

    print_peak_memory_usage();

    return 0;
}

} // namespace chopper::count