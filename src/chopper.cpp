#include <sharg/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/count/execute.hpp>
#include <chopper/layout/execute.hpp>

int main(int argc, const char *argv [])
{
    sharg::parser top_level_parser{"chopper", argc, argv, sharg::update_notifications::off, {"count", "layout"}};
    top_level_parser.info.version = "1.0.0";

    try
    {
        top_level_parser.parse();
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    sharg::parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

    int error_code{};

    if (sub_parser.info.app_name == std::string_view{"chopper-layout"})
        error_code = chopper::layout::execute(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"chopper-count"})
        error_code = chopper::count::execute(sub_parser);

    return error_code;
}
