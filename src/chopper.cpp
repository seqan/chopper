#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "chopper_data.hpp"
#include "minimizer.hpp"
#include "minimizer_msa.hpp"
#include "chopper_config.hpp"
#include "sequence_input.hpp"
#include "traverse_graph.hpp"

void set_up_argument_parser(seqan3::argument_parser & parser,
                            chopper_config & config,
                            cmd_arguments & traverse_config)
{
    parser.info.version = "1.0.0";
    parser.add_option(config.seqfiles, 's', "seq", "Name of multi-fasta input file.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(traverse_config.out_path, 'o', "outfile", "Name of the traversal output file.");
    parser.add_option(config.kmer_size, 'k', "kmer-size", "The kmer size to compute minimizer.");
    parser.add_option(config.window_size, 'w', "window-size", "The window size to compute minimizer.");
    parser.add_option(traverse_config.bins, 'b', "technical-bins", "How many technical bins do you want you sequences to be split?.");
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv [])
{
    chopper_config config;
    cmd_arguments traverse_config;

    // Command line parsing
    seqan3::argument_parser parser{"chopper", argc, argv};
    set_up_argument_parser(parser, config, traverse_config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    // Load data
    // -------------------------------------------------------------------------
    chopper_data data;

    auto start = std::chrono::steady_clock::now();
    for (auto const & file_name : config.seqfiles)
        if (!load_minimizer_sequences(data, config, file_name.c_str()))
            throw std::runtime_error{"Something went wrong when reading file " + file_name};
    auto end = std::chrono::steady_clock::now();
    seqan3::debug_stream << ">>> Loading " << seqan::length(data.sequences) << " sequences and computing minimizers complete "
                         << distance_matrix_initialiser::secs(start, end) << std::endl;

    // Compute minimizer MSA
    // -------------------------------------------------------------------------

    config.output_graph_file = "/tmp/graph.dot";
    minimizer_msa(data, config);
    // currently writes graph to an output file.

    // Traverse graph
    // -------------------------------------------------------------------------

    traverse_config.in_path = config.output_graph_file;
    traverse_config.write_out_graph = false;
    traverse_config.write_out_weights = false;

    traverse_graph(traverse_config);

    return 0;
}
