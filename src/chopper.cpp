#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "chopper_data.hpp"
#include "minimizer.hpp"
#include "minimizer_msa.hpp"
#include "chopper_config.hpp"
#include "sequence_input.hpp"
#include "traverse_graph.hpp"

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser,
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

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv [])
{
    chopper_config config;
    cmd_arguments traverse_config;

    seqan3::argument_parser top_level_parser{"chopper", argc, argv, false, {"split"}};
    top_level_parser.info.version = "1.0.0";

    try
    {
        top_level_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

    if (sub_parser.info.app_name == std::string_view{"chopper-split"})
    {
        if (auto r = set_up_and_parse_subparser_split(sub_parser, config, traverse_config); r != 0)
            return r;
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
