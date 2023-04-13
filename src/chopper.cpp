#include <sharg/parser.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>

int main(int argc, char const * argv[])
{
    sharg::parser parser{"chopper", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    chopper::configuration config;
    set_up_parser(parser, config);

    try
    {
        parser.parse();

        config.disable_sketch_output = !parser.is_option_set("output-sketches-to");
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    int exit_code{};

    chopper::layout::layout hibf_layout{};

    chopper::data_store store{.false_positive_rate = config.false_positive_rate, .hibf_layout = &hibf_layout};

    try
    {
        exit_code |= chopper::sketch::execute(config, store.filenames, store.sketches);
        chopper::sketch::estimate_kmer_counts(store.sketches, store.kmer_counts);
        exit_code |= chopper::layout::execute(config, store);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
        return -1;
    }

    return exit_code;
}
