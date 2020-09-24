#pragma once

#include <cstdlib>

#include "segment_generation_config.hpp" // include before #include <seqan/graph_msa.h> !!!

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/std/ranges>

#include "minimizer.hpp"
#include "ClusterAlgorithm.hpp"
#include "seqan2_msa_alignment.hpp"
#include "graph_output.hpp"

template <typename TNameSet, typename TSequence, typename TSize>
inline void minimizer_msa(seqan::StringSet<TSequence, seqan::Owner<>> & sequenceSet,
                          TNameSet & sequenceNames,
                          seqan::String<size_t> & sequenceLengths,
                          segment_generation_config<TSize> & seg_gen_config)
{
    // Alignment of the sequences
    typedef seqan::Graph<seqan::Alignment<seqan::StringSet<TSequence, seqan::Dependent<> >, void, seqan::WithoutEdgeId> > TGraph;
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
            outs << sequenceSet[id][regStart].position;
        outs << ",";
        auto regEnd = seqan::fragmentBegin(gAlign, *it) + seqan::fragmentLength(gAlign, *it);
        if (regEnd >= seqan::length(sequenceSet[id]))
            outs << sequenceLengths[id];
        else
            outs << sequenceSet[id][regEnd].position;
        outs << ")";
        outs << "\", group = ";
        outs << id;
        seqan::append(seqan::property(nodeMap, *it), outs.str().c_str());
        //std::cout << property(nodeMap, *it) << std::endl;
    }

    std::ofstream dotFile(seg_gen_config.output_graph_file);
    write_graph(dotFile, gAlign, nodeMap, edgeMap, sequenceNames, sequenceLengths, seqan::DotDrawing());
    dotFile.close();
}
