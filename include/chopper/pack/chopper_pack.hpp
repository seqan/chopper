#pragma once

#include <chopper/pack/aggregate_by.hpp>
#include <chopper/pack/hierarchical_binning.hpp>
#include <chopper/pack/filenames_data_input.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, pack_config & config)
{
    parser.info.version = "1.0.0";
    parser.add_option(config.data_file, 'f', "filenames", "A file with filenames.", seqan3::option_spec::REQUIRED);
    parser.add_option(config.bins, 'b', "technical-bins",
                      "In how many technical bins do you (hierarchical) IBF?.");
    parser.add_option(config.aggregate_by_column, 'y', "aggregate-by",
                      "Which column do you want to aggregate your files by? Start counting your columns from 1!",
                      seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{3, std::numeric_limits<int>::max()});

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER PACK ERROR] " << ext.what() << '\n'; // customize your error message
        return -2;
    }

    return 0;
}

int chopper_pack(seqan3::argument_parser & parser)
{
    pack_config config;

    if (auto r = set_up_and_parse_subparser_split(parser, config); r != 0)
        return r;

    pack_data data;
    read_filename_data_file(data, config);

    if (data.filenames.empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] File Error: you passed an empty file."};

    if (config.aggregate_by_column != -1 &&  data.extra_information[0].empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] Aggregate Error: You want to aggregate by something but your "
                                 "file does not contain any extra information columns."};

    if (config.aggregate_by_column > data.extra_information[0].size())
        throw std::runtime_error{"[CHOPPER PACK ERROR]Aggregate Error: You want to aggregate by a column index that is "
                                 "larger than the number of extra information columns."};

    aggregate_by(data, config.aggregate_by_column);

    hierarchical_binning algo{data.filenames, data.kmer_counts, config.bins};

    return 0;
}
