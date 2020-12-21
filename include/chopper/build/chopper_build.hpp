#pragma once

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/build/build_config.hpp>
#include <chopper/build/create_ibfs.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, build_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Build IBF on results from chopper-split.";
    parser.info.version = "1.0.0";

    parser.add_option(config.traversal_filename, 'p', "traversal", "Provide the traversal files from chopper split.");
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

    auto && [high_level_ibf, low_level_ibfs] = create_ibfs(config);

    {
        std::string const out_filename{config.output_prefix + "high_level.ibf"};
        std::ofstream fout(out_filename, std::ios::binary);

        if (!fout.good() || !fout.is_open())
            throw std::runtime_error{"Could not open " + out_filename + " for writing."};

        cereal::BinaryOutputArchive archive(fout);
        archive(high_level_ibf);
    }

    for (size_t i = 0; i < low_level_ibfs.size(); ++i)
    {
        if (low_level_ibfs[i].bin_size() != 1) // no dummy
        {
            std::string const out_filename{config.output_prefix + "low_level_" + std::to_string(i) + ".ibf"};
            std::ofstream fout(out_filename, std::ios::binary);

            if (!fout.good() || !fout.is_open())
                throw std::runtime_error{"Could not open " + out_filename + " for writing."};

            cereal::BinaryOutputArchive archive(fout);
            archive(low_level_ibfs[i]);
        }
    }

    return 0;
}

