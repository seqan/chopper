#pragma once

#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/std/ranges>

#include <chopper/split/split_config.hpp>
#include <chopper/split/split_data.hpp>
#include <chopper/split/map_distance_matrix.hpp>
#include <chopper/split/distance_matrix_initialiser.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/seqan2_msa_alignment.hpp>
#include <chopper/split/graph_output.hpp>

inline void minimizer_msa(split_data & data, batch_config const & config)
{
    // Alignment of the sequences
    typedef seqan::Graph<seqan::Alignment<seqan::StringSet<seqan::String<minimizer>, seqan::Dependent<> >, void, seqan::WithoutEdgeId> > TGraph;
    TGraph gAlign;

    distance_matrix_initialiser initialiser{};
    auto distance_matrix = initialiser.mash_distance(data, config);

    // MSA
    try
    {
        seqan2_msa_alignment(gAlign, data.sequences, distance_matrix, config);
    }
    catch (const std::bad_alloc & exception)
    {
        std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                  << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    seqan::String<seqan::String<char> > nodeMap;
    seqan::String<seqan::String<char> > edgeMap;
    seqan::_createEdgeAttributes(gAlign, edgeMap);

    // create node attributes
    typedef typename seqan::Id<TGraph>::Type TIdType;
    seqan::resizeVertexMap(nodeMap, gAlign);

    typedef typename seqan::Iterator<TGraph, seqan::VertexIterator>::Type TConstIter;
    TConstIter it(gAlign);
    for(;!seqan::atEnd(it);++it) {
        TIdType id = seqan::sequenceId(gAlign, *it);
        std::ostringstream outs;
        outs << "label = \"";
        outs << "[";
        auto regStart = seqan::fragmentBegin(gAlign, *it);
        if (regStart == 0)
            outs << "0"; // if it is the very first minimizer, include beginning of the sequence
        else
            outs << data.sequences[id][regStart].position;
        outs << ",";
        auto regEnd = seqan::fragmentBegin(gAlign, *it) + seqan::fragmentLength(gAlign, *it);
        if (regEnd >= seqan::length(data.sequences[id]))
            outs << data.lengths[id];
        else
            outs << data.sequences[id][regEnd].position;
        outs << ")";
        outs << "\", group = ";
        outs << id;
        seqan::append(seqan::property(nodeMap, *it), outs.str().c_str());
        //std::cout << property(nodeMap, *it) << std::endl;
    }

    std::ofstream dotFile(config.output_graph_file);
    write_graph(dotFile, gAlign, nodeMap, edgeMap, data.ids, data.lengths, seqan::DotDrawing());
    dotFile.close();
}
