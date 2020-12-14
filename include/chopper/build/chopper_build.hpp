#pragma once

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/build/build_config.hpp>
#include <chopper/build/create_ibfs_from_data_file.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, build_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Build IBF on results from chopper-split.";
    parser.info.version = "1.0.0";

    parser.add_option(config.binning_filename, 'f', "binning_filename", "Give me a filename to a seqinfo file.",
                      seqan3::option_spec::REQUIRED);
    parser.add_option(config.traversal_path_prefix, 'p', "traversal-prefix", "Give the prefix were the traversal files are stored.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to build kmers.");
    parser.add_option(config.overlap, 'l', "overlap", "The overlap between split regions of the same sequence.");
    parser.add_option(config.FPR, 'r', "false-positive-rate", "The minimum false positive rate of every IBF.");
    parser.add_option(config.output_prefix, 'o', "out-prefix", "Prefix of the output files.");
    parser.add_flag(config.verbose, 'v', "verbose", "Output logging/progress information.");
}

int chopper_build(seqan3::argument_parser & parser)
{
    build_config config{};
    initialize_argument_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[CHOPPER BIULD ERROR] " << ext.what() << "\n";
        return -1;
    }

    auto && [high_level_ibf, low_level_ibf_ids, low_level_ibfs] = create_ibfs_from_data_file(config);

    {
        std::string const out_filename{config.output_prefix + "high_level.ibf"};
        std::ofstream fout(out_filename, std::ios::binary);

        if (!fout.good() || !fout.is_open())
            throw std::runtime_error{"Could not open " + out_filename + " for writing."};

        cereal::BinaryOutputArchive archive(fout);
        archive(high_level_ibf);
    }

    assert(low_level_ibfs.size() == low_level_ibf_ids.size());
    for (size_t i = 0; i < low_level_ibf_ids.size(); ++i)
    {
        std::string const out_filename{config.output_prefix + "low_level_" + low_level_ibf_ids[i] + ".ibf"};
        std::ofstream fout(out_filename, std::ios::binary);

        if (!fout.good() || !fout.is_open())
            throw std::runtime_error{"Could not open " + out_filename + " for writing."};

        cereal::BinaryOutputArchive archive(fout);
        archive(low_level_ibfs[i]);
    }

    return 0;
}

