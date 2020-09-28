#pragma once

#include <chopper/pack/filenames_data_input.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, pack_config & config)
{
    parser.info.version = "1.0.0";
    parser.add_option(config.data_file, 'f', "filenames", "A file with filenames.", seqan3::option_spec::REQUIRED);
    parser.add_option(config.bins, 'b', "technical-bins",
                      "In how many technical bins do you (hierarchical) IBF?.");

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

    return 0;
}
