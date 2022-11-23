#pragma once

#include <chopper/configuration.hpp>
#include <chopper/prefixes.hpp>
#include <chopper/count/check_filenames.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/print_peak_memory_usage.hpp>

namespace chopper::count
{

int execute(chopper::configuration & config)
{
    auto filename_clusters = chopper::count::read_data_file(config);

    chopper::count::check_filenames(filename_clusters, config);

    chopper::count::count_kmers(filename_clusters, config);

    chopper::print_peak_memory_usage();

    return 0;
}

} // namespace chopper::count
