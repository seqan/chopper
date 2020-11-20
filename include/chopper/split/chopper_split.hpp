#pragma once

#include <seqan3/core/debug_stream.hpp>

#include <chopper/split/split_data.hpp>
#include <chopper/split/split_config.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/minimizer_msa.hpp>
#include <chopper/split/sequence_input.hpp>
#include <chopper/split/traverse_graph.hpp>
#include <chopper/split/filename_batches_range.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, split_config & config)
{
    parser.info.version = "1.0.0";
    parser.add_option(config.data_filename, 'f', "binning-file", "A high_level_ibf.binning or low_level_ibfs.binning file.");
    parser.add_option(config.seqfiles, 's', "seq", "Name of multi-fasta input file.");
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

    return 0;
}

// runs splitting of sequences into technical bins
int chopper_split(seqan3::argument_parser & parser)
{
    split_config config;

    if (auto r = set_up_and_parse_subparser_split(parser, config); r != 0)
    return r;

    // defaults set for every config
    config.output_graph_file = "/tmp/graph.dot";
    config.write_out_graph = false;
    config.write_out_weights = false;

    // Check config
    // -------------------------------------------------------------------------
    if (!config.data_filename.empty() && !config.seqfiles.empty())
        throw std::runtime_error{"[CHOPPER SPLIT ERROR] You may EITHER specify files with -s OR give a data file "
                                 "with -f!"};
    if (config.data_filename.empty() && config.seqfiles.empty())
        throw std::runtime_error{"[CHOPPER SPLIT ERROR] You must specify EITHER files with -s OR give a data file "
                                 "with -f!"};

    for (auto const & batch_config : filename_batches_range{config})
    {
        std::cout << "Processing file " << batch_config.out_path << " with " << batch_config.bins << " number of bins." << std::endl;

        if (batch_config.bins == 1) // nothing to split here
        {
            std::string const low_level_prefix{config.out_path.string() + "LOW_LEVEL_IBF_"};

            if (starts_with(batch_config.out_path, low_level_prefix))
            {
                bool const output_file_exists = std::filesystem::exists(batch_config.out_path);

                std::ofstream fout{batch_config.out_path, std::ios::binary | std::ios::app}; // append to file

                if (!output_file_exists)
                    fout << "FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER\n"; // header

                for (auto const & filename : batch_config.seqfiles)
                {
                    seqan3::sequence_file_input fin{filename, seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

                    for (auto const & [id, seq] : fin)
                    {
                        fout << filename << '\t'
                             << id << '\t'
                             << 0 << '\t'
                             << seq.size() << '\t'
                             << batch_config.bin_index_offset << '\n';
                    }
                }
            }

            continue;
        }

        // Load data
        // -------------------------------------------------------------------------
        split_data data;

        auto start = std::chrono::steady_clock::now();
        for (auto const & file_name : batch_config.seqfiles)
            if (!load_minimizer_sequences(data, batch_config, file_name.c_str()))
                throw std::runtime_error{"Something went wrong when reading file " + file_name};
        auto end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Loading " << seqan::length(data.sequences)
                             << " sequences and computing minimizers complete "
                             << distance_matrix_initialiser::secs(start, end) << std::endl;

        // Compute minimizer MSA
        // -------------------------------------------------------------------------

        minimizer_msa(data, batch_config);
        // currently writes graph to an output file.

        // Traverse graph
        // -------------------------------------------------------------------------
        traverse_graph(data, batch_config);
    }

    return 0;
}
