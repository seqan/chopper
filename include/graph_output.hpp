#pragma once

template <typename TTarget, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes, typename TNames, typename TSeq>
void write_graph(TTarget & target,
                 seqan::Graph<TSpec> const& g,
                 TNodeAttributes const& nodeMap,
                 TEdgeAttributes const& edgeMap,
                 TNames const & names,
                 TSeq const & seq_lenths,
                 seqan::DotDrawing)
{
    typedef seqan::Graph<TSpec> TGraph;
    typedef typename seqan::VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typename seqan::DirectionIterator<TTarget, seqan::Output>::Type iter = directionIterator(target, seqan::Output());

    seqan::_writeGraphType(iter,g,seqan::DotDrawing());
    seqan::write(iter, " G {\n");
    seqan::writeValue(iter, '\n');
    seqan::write(iter, "/* Graph Attributes */\n");
    seqan::write(iter, "graph [rankdir = LR];\n");
    seqan::writeValue(iter, '\n');
    seqan::write(iter, "/* Node Attributes */\n");
    seqan::write(iter, "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
    seqan::writeValue(iter, '\n');
    seqan::write(iter, "/* Edge Attributes */\n");
    seqan::write(iter, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
    seqan::writeValue(iter, '\n');
    seqan::write(iter, "/* Sequence Lengths */\n");
    for (size_t i = 0; i < seqan::length(names); ++i)
    {
        seqan::write(iter, names[i]);
        seqan::writeValue(iter, '\t');
        seqan::write(iter, seq_lenths[i]);
        seqan::writeValue(iter, '\n');
    }
    seqan::writeValue(iter, '\n');

    seqan::write(iter, "/* Nodes */\n");
    typedef typename seqan::Iterator<TGraph, seqan::VertexIterator>::Type TConstIter;
    TConstIter it(g);
    for(;!seqan::atEnd(it);++it) {
        seqan::appendNumber(iter, (int)*it);
        seqan::write(iter, " [");
        seqan::write(iter, seqan::getProperty(nodeMap, *it));
        seqan::write(iter, "];\n");
    }
    seqan::writeValue(iter, '\n');

    seqan::write(iter, "/* Edges */\n");
    typedef typename seqan::Iterator<TGraph, seqan::EdgeIterator>::Type TConstEdIter;
    TConstEdIter itEd(g);
    for(;!seqan::atEnd(itEd);++itEd) {
        TVertexDescriptor sc = seqan::sourceVertex(itEd);
        TVertexDescriptor tr = seqan::targetVertex(itEd);
        seqan::appendNumber(iter, sc);
        seqan::_writeEdgeType(iter, g, seqan::DotDrawing());
        seqan::appendNumber(iter, tr);
        seqan::write(iter, " [");
        seqan::write(iter, seqan::getProperty(edgeMap, *itEd));
        seqan::write(iter, "];\n");
    }
    seqan::writeValue(iter, '\n');

    seqan::_writeGraphFooter(iter,g,seqan::DotDrawing());

    seqan::write(iter, "}\n");
}
