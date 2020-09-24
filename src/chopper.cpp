#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "minimizer.hpp"
#include "minimizer_msa.hpp"
#include "segment_generation_config.hpp"
#include "sequence_input.hpp"
#include "traverse_graph.hpp"

template <typename TSize>
void set_up_argument_parser(seqan3::argument_parser & parser,
                            segment_generation_config<TSize> & seg_gen_config,
                            cmd_arguments & traverse_config)
{
    parser.info.version = "1.0.0";
    parser.add_option(seg_gen_config.seqfiles, 's', "seq", "Name of multi-fasta input file.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(traverse_config.out_path, 'o', "outfile", "Name of the traversal output file.");
    parser.add_option(seg_gen_config.kmer_size, 'k', "kmer-size", "The kmer size to compute minimizer.");
    parser.add_option(seg_gen_config.window_size, 'w', "window-size", "The window size to compute minimizer.");
    parser.add_option(traverse_config.bins, 'b', "technical-bins", "How many technical bins do you want you sequences to be split?.");
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv [])
{
    using TSize = typename seqan::Size<seqan::StringSet<seqan::String<minimizer>, seqan::Dependent<> >>::Type;
    segment_generation_config<TSize> seg_gen_config;
    cmd_arguments traverse_config;

    // Command line parsing
    seqan3::argument_parser parser{"chopper", argc, argv};
    set_up_argument_parser(parser, seg_gen_config, traverse_config);

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

    typedef seqan::String<minimizer> TSequence;
    seqan::StringSet<TSequence, seqan::Owner<> > sequenceSet;
    seqan::StringSet<seqan::String<char> > sequenceNames;
    seqan::String<size_t> sequenceLengths;

    auto start = std::chrono::steady_clock::now();
    for (auto const & file_name : seg_gen_config.seqfiles)
        if (!load_minimizer_sequences(sequenceSet, sequenceNames, sequenceLengths, seg_gen_config, file_name.c_str()))
            throw std::runtime_error{"Something went wrong when reading file " + file_name};
    auto end = std::chrono::steady_clock::now();
    seqan3::debug_stream << ">>> Loading " << seqan::length(sequenceSet) << " sequences and computing minimizers complete "
                         << distance_matrix_initialiser::secs(start, end) << std::endl;

    // Compute minimizer MSA
    // -------------------------------------------------------------------------

    seg_gen_config.output_graph_file = "/tmp/graph.dot";
    minimizer_msa(sequenceSet, sequenceNames, sequenceLengths, seg_gen_config);
    // currently writes graph to an output file.

    // Traverse graph
    // -------------------------------------------------------------------------

    traverse_config.in_path = seg_gen_config.output_graph_file;
    traverse_config.write_out_graph = false;
    traverse_config.write_out_weights = false;

    traverse_graph(traverse_config);

    return 0;
}
