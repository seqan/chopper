#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "minimizer.hpp"
#include "minimizer_msa.hpp"
#include "segment_generation_config.hpp"

template <typename TSize>
void set_up_argument_parser(seqan3::argument_parser & parser, segment_generation_config<TSize> & seg_gen_config)
{
    parser.info.version = "1.0.0";
    parser.add_option(seg_gen_config.seqfiles, 's', "seq", "Name of multi-fasta input file.",
                      seqan3::option_spec::REQUIRED);
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv [])
{
    using TSize = typename seqan::Size<seqan::StringSet<seqan::String<minimizer>, seqan::Dependent<> >>::Type;
    segment_generation_config<TSize> seg_gen_config;

    // Command line parsing
    seqan3::argument_parser parser{"chopper", argc, argv};
    set_up_argument_parser(parser, seg_gen_config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    minimizer_msa(seg_gen_config);

    return 0;
}
