#pragma once

#include <chopper/split/split_data.hpp>
#include <chopper/split/split_config.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/minimizer_msa.hpp>
#include <chopper/split/sequence_input.hpp>
#include <chopper/split/traverse_graph.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser,
                                     split_config & config)
{
    parser.info.version = "1.0.0";
    parser.add_option(config.seqfiles, 's', "seq", "Name of multi-fasta input file.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(config.out_path, 'o', "outfile", "Name of the traversal output file.");
    parser.add_option(config.kmer_size, 'k', "kmer-size", "The kmer size to compute minimizer.");
    parser.add_option(config.window_size, 'w', "window-size", "The window size to compute minimizer.");
    parser.add_option(config.bins, 'b', "technical-bins", "How many technical bins do you want you sequences to be split?.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER SPLIT ERROR] " << ext.what() << '\n'; // customize your error message
        return -2;
    }
}

// runs splitting of sequences into technical bins
int chopper_split(seqan3::argument_parser & parser)
{
    split_config config;

    if (auto r = set_up_and_parse_subparser_split(parser, config); r != 0)
        return r;

    // Load data
    // -------------------------------------------------------------------------
    split_data data;

    auto start = std::chrono::steady_clock::now();
    for (auto const & file_name : config.seqfiles)
        if (!load_minimizer_sequences(data, config, file_name.c_str()))
            throw std::runtime_error{"Something went wrong when reading file " + file_name};
    auto end = std::chrono::steady_clock::now();
    seqan3::debug_stream << ">>> Loading " << seqan::length(data.sequences)
                         << " sequences and computing minimizers complete "
                         << distance_matrix_initialiser::secs(start, end) << std::endl;

    // Compute minimizer MSA
    // -------------------------------------------------------------------------

    config.output_graph_file = "/tmp/graph.dot";
    minimizer_msa(data, config);
    // currently writes graph to an output file.

    // Traverse graph
    // -------------------------------------------------------------------------

    config.write_out_graph = false;
    config.write_out_weights = false;

    traverse_graph(config);

    return 0;
}
