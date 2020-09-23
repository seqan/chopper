#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/std/ranges>

#include "include/minimizer.hpp"
#include "include/ClusterAlgorithm.hpp"
#include "include/seqan2_msa_alignment.hpp"
#include "include/graph_output.hpp"
#include "include/segment_generation_config.hpp"

using namespace seqan;

template <typename TSize>
inline void
customizedMsaAlignment(segment_generation_config<TSize> & seg_gen_config)
{
    typedef String<minimizer> TSequence;
    StringSet<TSequence, Owner<> > sequenceSet;
    StringSet<String<char> > sequenceNames;
    String<size_t> sequenceLengths;

    // Alignment of the sequences
    typedef Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > TGraph;
    TGraph gAlign;

    ClusterAlgorithm alg{};
    alg.fill_mash_distance_matrix(sequenceSet, sequenceNames, sequenceLengths, seg_gen_config);

    // MSA
    try
    {
        seqan2_msa_alignment(gAlign, sequenceSet, seg_gen_config);
    }
    catch (const std::bad_alloc & exception)
    {
        std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                  << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    String<String<char> > nodeMap;
    String<String<char> > edgeMap;
    _createEdgeAttributes(gAlign, edgeMap);

    // create node attributes
    typedef typename Id<TGraph>::Type TIdType;
    resizeVertexMap(nodeMap, gAlign);

    typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
    TConstIter it(gAlign);
    for(;!atEnd(it);++it) {
        TIdType id = sequenceId(gAlign, *it);
        std::ostringstream outs;
        outs << "label = \"";
        outs << "[";
        auto regStart = fragmentBegin(gAlign, *it);
        if (regStart == 0)
            outs << "0"; // if it is the very first minimizer, include beginning of the sequence
        else
            outs << sequenceSet[id][regStart].position;
        outs << ",";
        auto regEnd = fragmentBegin(gAlign, *it) + fragmentLength(gAlign, *it);
        if (regEnd >= length(sequenceSet[id]))
            outs << sequenceLengths[id];
        else
            outs << sequenceSet[id][regEnd].position;
        outs << ")";
        outs << "\", group = ";
        outs << id;
        append(property(nodeMap, *it), outs.str().c_str());
        //std::cout << property(nodeMap, *it) << std::endl;
    }

    std::ofstream dotFile("graph.dot");
    write_graph(dotFile, gAlign, nodeMap, edgeMap, sequenceNames, sequenceLengths, DotDrawing());
    dotFile.close();
}

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
    using TSize = typename Size<StringSet<String<minimizer>, Dependent<> >>::Type;
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

    customizedMsaAlignment(seg_gen_config);

    return 0;
}
