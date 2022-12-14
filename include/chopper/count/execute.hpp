#pragma once

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>

#include <chopper/configuration.hpp>
#include <chopper/count/check_filenames.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/data_store.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/prefixes.hpp>
#include <chopper/print_peak_memory_usage.hpp>

namespace chopper::count
{

inline int execute(configuration & config, data_store & data)
{
    read_data_file(config, data);

    if (data.filenames.empty())
        throw sharg::parser_error{
            sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

    chopper::count::check_filenames(data.filenames, config);

    count_kmers(config, data);

    print_peak_memory_usage();

    return 0;
}

} // namespace chopper::count
