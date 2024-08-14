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

#include <chopper/chopper_layout.hpp>
#include <chopper/configuration.hpp>
#include <chopper/input_functor.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/sketch/check_filenames.hpp>
#include <chopper/sketch/output.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include <hibf/sketch/compute_sketches.hpp>

namespace chopper
{

int chopper_layout(chopper::configuration & config, sharg::parser & parser)
{
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

    std::vector<std::vector<std::string>> filenames{};

    chopper::sketch::read_data_file(config, filenames);

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;

    try
    {
        if (filenames.empty())
            throw sharg::parser_error{
                sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

        chopper::sketch::check_filenames(filenames, config);

        config.hibf_config.input_fn =
            chopper::input_functor{filenames, config.precomputed_files, config.k, config.window_size};
        config.hibf_config.number_of_user_bins = filenames.size();

        config.compute_sketches_timer.start();
        seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches);
        config.compute_sketches_timer.stop();

        exit_code |= chopper::layout::execute(config, filenames, sketches);
    }
    catch (std::exception const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
        return -1;
    }

    if (!config.disable_sketch_output)
    {
        if (!std::filesystem::exists(config.sketch_directory))
            std::filesystem::create_directory(config.sketch_directory);

        assert(filenames.size() == sketches.size());
        for (size_t i = 0; i < filenames.size(); ++i)
            chopper::sketch::write_sketch_file(filenames[i][0], sketches[i], config);
    }

    if (!config.output_timings.empty())
    {
        std::ofstream output_stream{config.output_timings};
        output_stream << std::fixed << std::setprecision(2);
        output_stream << "sketching_in_seconds\t"
                      << "layouting_in_seconds\t"
                      << "union_estimation_in_seconds\t"
                      << "rearrangement_in_seconds\n";
        output_stream << config.compute_sketches_timer.in_seconds() << '\t';
        output_stream << config.dp_algorithm_timer.in_seconds() << '\t';
        output_stream << config.union_estimation_timer.in_seconds() << '\t';
        output_stream << config.rearrangement_timer.in_seconds() << '\t';
    }

    return exit_code;
}

} // namespace chopper
