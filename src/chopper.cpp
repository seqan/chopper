#include <sharg/parser.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <chopper/configuration.hpp>
#include <chopper/count/execute.hpp>
#include <chopper/data_store.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>

int main(int argc, char const * argv[])
{
    sharg::parser parser{"chopper", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    chopper::configuration config;
    set_up_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    config.input_prefix = config.output_prefix;

    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    int exit_code{};


    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .output_buffer = &output_buffer,
                             .header_buffer = &header_buffer};

    try
    {
        exit_code |= chopper::count::execute(config);
        exit_code |= chopper::layout::execute(config, data);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
        return -1;
    }

    return exit_code;
}
