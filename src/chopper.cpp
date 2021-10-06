#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/count/chopper_count.hpp>
#include <chopper/pack/chopper_pack.hpp>
#include <chopper/split/chopper_split.hpp>

int main(int argc, const char *argv [])
{
    seqan3::argument_parser top_level_parser{"chopper", argc, argv, seqan3::update_notifications::off,
                                             {"count", "pack", "split"}};
    top_level_parser.info.version = "1.0.0";

    try
    {
        top_level_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

    int error_code{};

    if (sub_parser.info.app_name == std::string_view{"chopper-split"})
        error_code = chopper_split(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"chopper-pack"})
        error_code = chopper_pack(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"chopper-count"})
        error_code = chopper_count(sub_parser);

    return error_code;
}
