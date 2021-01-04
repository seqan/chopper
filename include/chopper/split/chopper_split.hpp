#pragma once

#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <chopper/split/distance_matrix_initialiser.hpp>
#include <chopper/split/filename_batches_range.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/seqan2_msa_alignment.hpp>
#include <chopper/split/sequence_input.hpp>
#include <chopper/split/split_config.hpp>
#include <chopper/split/split_data.hpp>
#include <chopper/split/transform_graphs.hpp>
#include <chopper/split/traverse_graph.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, split_config & config)
{
    parser.info.version = "1.0.0";
    parser.add_flag(config.verbose, 'v', "verbose", "Ooutput verbose logging information.");
    parser.add_option(config.data_filename, 'f', "binning-file", "A high_level_ibf.binning or low_level_ibfs.binning file.");
    parser.add_option(config.seqfiles, 's', "seq", "Name of multi-fasta input file.");
    parser.add_option(config.out_path, 'o', "outfile", "Name of the output file (file extension should be .split).");
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

    return 0;
}

// runs splitting of sequences into technical bins
int chopper_split(seqan3::argument_parser & parser)
{
    split_config config;

    if (auto r = set_up_and_parse_subparser_split(parser, config); r != 0)
        return r;

    // Check config
    // -------------------------------------------------------------------------
    if (!config.data_filename.empty() && !config.seqfiles.empty())
        throw std::runtime_error{"[CHOPPER SPLIT ERROR] You may EITHER specify files with -s OR give a data file "
                                 "with -f!"};
    if (config.data_filename.empty() && config.seqfiles.empty())
        throw std::runtime_error{"[CHOPPER SPLIT ERROR] You must specify EITHER files with -s OR give a data file "
                                 "with -f!"};

    std::ofstream fout{config.out_path};

    // transfer header from binning file to chopper_split
    {
        std::ifstream data_file{config.data_filename};
        std::string line;
        while (std::getline(data_file, line) && line.substr(0, 7) != "#BIN_ID")
            fout << line << '\n';
    }

    fout << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n";

    for (auto const & current_batch_config : filename_batches_range{config})
    {
        std::cout << "Processing file " << current_batch_config.bin_name << " with " << current_batch_config.bins << " number of bins." << std::endl;

        if (current_batch_config.bins == 1) // nothing to split here
        {
            for (auto const & filename : current_batch_config.seqfiles)
            {
                seqan3::sequence_file_input fin{filename, seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

                for (auto const & [id, seq] : fin)
                {
                    fout << filename << '\t'
                         << id << '\t'
                         << 0 << '\t'
                         << seq.size() << '\t'
                         << current_batch_config.hibf_bin_idx_offset << '\t';
                    if (current_batch_config.merged_bin)
                        fout << current_batch_config.libf_bin_idx_offset;
                    else
                        fout << '-';
                    fout << '\n';
                }
            }

            continue;
        }

        // Load data
        // -------------------------------------------------------------------------
        split_data data;
        data.outstream = &fout;

        auto start = std::chrono::steady_clock::now();
        for (auto const & file_name : current_batch_config.seqfiles)
            if (!load_minimizer_sequences(data, current_batch_config, file_name.c_str()))
                throw std::runtime_error{"Something went wrong when reading file " + file_name};
        auto end = std::chrono::steady_clock::now();

        if (config.verbose)
        {
            seqan3::debug_stream << ">>> Loading " << seqan::length(data.sequences)
                                 << " sequences and computing minimizers complete "
                                 << distance_matrix_initialiser<>::secs(start, end) << std::endl;
        }

        // Compute minimizer MSA
        // -------------------------------------------------------------------------

        // Alignment of the sequences
        typedef seqan::Graph<seqan::Alignment<seqan::StringSet<seqan::String<minimizer>, seqan::Dependent<> >, void, seqan::WithoutEdgeId> > TGraph;
        TGraph seqan2_graph;

        distance_matrix_initialiser<seqan::String<double>> initialiser{seqan::length(data.sequences)};
        auto distance_matrix = initialiser.mash_distance(data, config);

        // MSA
        try
        {
            seqan2_msa_alignment(seqan2_graph, data.sequences, distance_matrix, config);
        }
        catch (const std::bad_alloc & exception)
        {
            std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                      << exception.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // Transform seqan2 graph to lemon graph
        // -------------------------------------------------------------------------
        lemon::ListDigraph g;
        std::vector<lemon::ListDigraph::Node> nodes;
        lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> node_map{g};

        transform_graphs(g, nodes, node_map, seqan2_graph, data, current_batch_config);

        // Traverse graph
        // -------------------------------------------------------------------------
        traverse_graph(g, nodes, node_map, data, current_batch_config);
    }

    return 0;
}
