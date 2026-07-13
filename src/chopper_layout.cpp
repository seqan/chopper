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
#include <chopper/sketch/sketch_file.hpp>

#include <hibf/sketch/compute_sketches.hpp>

namespace chopper
{

/*!\brief Validates the chopper configuration of the current chopper call with that from the sketch file
 *
 * If a sketch file is given as input, certain configurations have been already used for generating the sketches but
 * the user can still give the same arguments to the current command line.
 * This function checks for consistency als overwrites default values within the chopper configuration with those from the
 * sketch file. Overwriting is needed because the current chopper config is written to the layout file and should
 * be consistent with the sketch file.
 */
void validate_configuration(sharg::parser & parser,
                            chopper::configuration & config,
                            chopper::configuration const & sketch_config)
{
    if (parser.is_option_set("sketch-bits"))
        throw sharg::parser_error{"You cannot set --sketch-bits when using a sketch file as input."};

    if (parser.is_option_set("kmer") && config.k != sketch_config.k)
    {
        std::cerr << sharg::detail::to_string(
            "[WARNING] Given k-mer size (",
            config.k,
            ") differs from k-mer size in the sketch file (",
            sketch_config.k,
            "). The results may be suboptimal. If this was a conscious decision, you can ignore this warning.\n");
    }
    else
    {
        config.k = sketch_config.k; // make sure default is overwritten
    }

    if (parser.is_option_set("window") && config.window_size != sketch_config.window_size)
    {
        throw sharg::parser_error{"You are using a sketch file as input but want to use a --window size that differs "
                                  "from the one used for sketching. This will result in a layout (and subsequently "
                                  "an index that does not represent your data correctly."};
    }
    else
    {
        config.window_size = sketch_config.window_size; // make sure default is overwritten
    }
}

int chopper_layout(chopper::configuration & config, sharg::parser & parser)
{
    parser.parse();

    if (!parser.is_option_set("window"))
        config.window_size = config.k;
    else if (config.k > config.window_size)
        throw sharg::parser_error{"The k-mer size cannot be bigger than the window size."};

    auto has_sketch_file_extension = [](std::filesystem::path const & path)
    {
        return path.string().ends_with(".sketch") || path.string().ends_with(".sketches");
    };

    config.disable_sketch_output = !parser.is_option_set("output-sketches-to");
    if (!config.disable_sketch_output && !has_sketch_file_extension(config.sketch_directory))
        throw sharg::parser_error{"The sketch output file must have the extension \".sketch\" or \".sketches\"."};

    bool const input_is_a_sketch_file = has_sketch_file_extension(config.data_file);

    int exit_code{};

    std::vector<std::vector<std::string>> filenames{};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};

    if (input_is_a_sketch_file)
    {
        chopper::sketch::sketch_file sin{};

        { // Deserialization is guaranteed to be complete when going out of scope.
            std::ifstream is{config.data_file};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(sin);
        }

        filenames = std::move(sin.filenames); // No need to call check_filenames because the files are not read.
        sketches = std::move(sin.hll_sketches);
        minHash_sketches = std::move(sin.minHash_sketches);
        validate_configuration(parser, config, sin.chopper_config);
    }
    else
    {
        chopper::sketch::read_data_file(config, filenames);

        if (filenames.empty())
            throw sharg::parser_error{
                sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

        // Files need to exist because they will be read for sketching.
        chopper::sketch::check_filenames(filenames, config);
    }

    config.hibf_config.input_fn =
        chopper::input_functor{filenames, config.precomputed_files, config.k, config.window_size};
    config.hibf_config.number_of_user_bins = filenames.size();
    config.hibf_config.validate_and_set_defaults();

    if (!input_is_a_sketch_file)
    {
        config.compute_sketches_timer.start();
        seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);
        config.compute_sketches_timer.stop();
    }

    exit_code |= chopper::layout::execute(config, filenames, sketches, minHash_sketches);

    if (!config.disable_sketch_output)
    {
        chopper::sketch::sketch_file sout{.chopper_config = config,
                                          .filenames = std::move(filenames),
                                          .hll_sketches = std::move(sketches),
                                          .minHash_sketches = std::move(minHash_sketches)};
        std::ofstream os{config.sketch_directory, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(sout);
    }

    if (!config.output_timings.empty())
    {
        std::ofstream output_stream{config.output_timings};
        output_stream << std::fixed << std::setprecision(2);
        output_stream << "sketching_in_seconds\t"
                      << "layouting_in_seconds\t"
                      << "union_estimation_in_seconds\t"
                      << "rearrangement_in_seconds\t"
                      << "lsh_in_seconds\t"
                      << "intital_partition_timer_in_seconds\t"
                      << "small_layouts_timer_in_seconds\t"
                      << "search_best_p_in_seconds\n";
        output_stream << config.compute_sketches_timer.in_seconds() << '\t';
        output_stream << config.dp_algorithm_timer.in_seconds() << '\t';
        output_stream << config.union_estimation_timer.in_seconds() << '\t';
        output_stream << config.rearrangement_timer.in_seconds() << '\t';
        output_stream << config.lsh_algorithm_timer.in_seconds() << '\t';
        output_stream << config.intital_partition_timer.in_seconds() << '\t';
        output_stream << config.small_layouts_timer.in_seconds() << '\t';
        output_stream << config.search_partition_algorithm_timer.in_seconds() << '\n';
    }

    return exit_code;
}

} // namespace chopper
