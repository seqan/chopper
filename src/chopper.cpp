// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <sharg/parser.hpp>

#include <chopper/configuration.hpp>
#include <chopper/input_functor.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/check_filenames.hpp>
#include <chopper/sketch/read_data_file.hpp>

int main(int argc, char const * argv[])
{
    sharg::parser parser{"chopper", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    chopper::configuration config;
    set_up_parser(parser, config);
    parser.info.synopsis.front().insert(0, "chopper");

    try
    {
        parser.parse();

        if (!parser.is_option_set("window"))
            config.window_size = config.k;
        else if (config.k > config.window_size)
            throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

        config.disable_sketch_output = !parser.is_option_set("output-sketches-to");
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    int exit_code{};

    std::vector<std::string> filenames{};

    chopper::sketch::read_data_file(config, filenames);

    try
    {
        if (filenames.empty())
            throw sharg::parser_error{
                sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

        chopper::sketch::check_filenames(filenames, config);

        config.hibf_config.input_fn =
            chopper::input_functor{filenames, config.precomputed_files, config.k, config.window_size};
        config.hibf_config.number_of_user_bins = filenames.size();

        exit_code |= chopper::layout::execute(config, filenames);
    }
    catch (std::exception const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
        return -1;
    }

    return exit_code;
}
