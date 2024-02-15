// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <iostream>
#include <set>

#include <robin_hood.h>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/sketch/hyperloglog.hpp>

struct config
{
    std::filesystem::path input_index{};
    std::filesystem::path general_stats_output{"general.stats"};
};

struct stats
{
    std::vector<size_t> ibf_sizes;
    std::vector<size_t> ibf_levels;
    std::vector<double> ibf_load_factor;
};

inline void set_up_parser(sharg::parser & parser, config & cfg)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Compute a top-level HIBF layout figure file";

    parser.info.description.emplace_back("Computes an table to display the top-level layout.");

    parser.add_subsection("Main options:");
    parser.add_option(cfg.input_index,
                      sharg::config{.short_id = '\0',
                                    .long_id = "index",
                                    .description = "The input must be an index computed via raptor layout/build. ",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(cfg.general_stats_output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"layout_stats", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    config cfg{};
    set_up_parser(parser, cfg);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    auto index = raptor::raptor_index<index_structure::hibf>{};
    raptor::load_index(index, arguments);
}
