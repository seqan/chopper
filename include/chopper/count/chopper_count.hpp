#pragma once

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/count/count_config.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, count_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Count all kmers of each file in a directory.";
    parser.info.version = "1.0.0";

    parser.add_option(config.data_file, 'f', "data_file", "Give me a filename to a seqinfo file.", seqan3::option_spec::REQUIRED);
    parser.add_option(config.column_index_to_cluster, 'c', "column-index", "The column index by which to cluster.");
    parser.add_option(config.num_threads, 't', "threads", "Number of threads.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to count kmers.");
    parser.add_option(config.w, 'w', "window-size", "The window size for minimisers.");
}

int chopper_count(seqan3::argument_parser & parser)
{
    count_config config{};
    initialize_argument_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[CHOPPER COUNT ERROR] " << ext.what() << "\n";
        return -1;
    }

    auto filename_clusters = read_data_file(config);

    count_kmers(filename_clusters, config);

    return 0;
}

