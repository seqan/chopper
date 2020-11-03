#pragma once

#include <map>

#include <seqan/graph_msa.h>

namespace seqan
{

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPos>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type
find_vertex(Graph<Alignment<TStringSet, TCargo, TSpec>> const & g, TSeqId const id, TPos const pos)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename TGraph::TKey_ TKey;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TIdType;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    return (pos >= (TPos) length(getValueById(stringSet(g), id))) ?
                getNil<TVertexDescriptor>() :
                g.data_pvMap.upper_bound(TKey((TIdType)id, (TSize)pos))->second;
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TPosition, typename TSequence>
inline void build_leaf_string(Graph<Alignment<TStringSet, TCargo, TSpec>> const & g,
                              TPosition const pos,
                              TSequence & alignSeq)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TId;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TSequence>::Type TVertexString;

    TStringSet const & str = stringSet(g);
    TId const seqId = positionToId(str, pos);
    TSize const lenRoot = length(str[pos]);
    TSize i = 0;

    while (i < lenRoot)
    {
        TVertexDescriptor nextVertex = find_vertex(g, seqId, i);
        TVertexString vs;
        appendValue(vs, nextVertex);
        appendValue(alignSeq, vs, Generous());
        i += fragmentLength(g, nextVertex);
    }
}

template<typename TString, typename TWeightMap, typename TPositions>
inline void heaviest_increasing_subsequence(TString const & str,
                                            TWeightMap const & weights,
                                            TPositions& pos)
{
    typedef typename Size<TString>::Type TSize;
    typedef typename Value<TString>::Type TValue;
    typedef typename Value<TPositions>::Type TPos;
    typedef typename Value<TWeightMap>::Type TWeight;

    // The list of decreasing covers, only the smallest element of each member must be remembered
    typedef std::pair<TValue, std::pair<TWeight, TPos> > TKey;
    typedef std::set<TKey, std::less<TKey> > TSortedSequence;
    typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
    TSortedSequence list;

    // The trace-back graph
    typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TGraph g;

    // Walk through the sequence and build the decreasing covers
    typedef typename Iterator<TString const, Standard>::Type TStringIter;
    TStringIter it = begin(str, Standard());
    TStringIter endIt = end(str, Standard());
    TSize pos_of_iterator = 0;

    for (; it != endIt; ++it, ++pos_of_iterator)
    {
        TWeight w = weights[pos_of_iterator];
        // Letters that do not contribute a weight (e.g., w = 0) are excluded!
        // Weights must increase!
        if (w == 0)
        {
            addVertex(g);  // Note: The vertex id corresponds to the position
            continue;
        }

        // Get previous element
        TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, std::make_pair(0, 0)));

        // Get next element
        TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);

        // Determine new weight
        if (a_k_it != list.end())
            w += a_k_it->second.first;

        // Delete from list
        while ((b_l_it != list.end()) && (w >= b_l_it->second.first))
        {
            TSortedSequenceIter tmp = b_l_it;
            b_l_it = _nextInSortedSequence(list, b_l_it);
            list.erase(*tmp);
        }

        // Insert new list element
        if ((b_l_it == list.end()) || (*it < b_l_it->first))
        {
            list.insert(std::make_pair(*it, std::make_pair(w, pos_of_iterator)));
        }

        // Create the corresponding node, pos_of_iterator == Vertex Descriptor
        addVertex(g);

        // Connect to predecessor
        if (a_k_it != list.end())
            addEdge(g, (TVertexDescriptor) pos_of_iterator, (TVertexDescriptor) a_k_it->second.second);
    }

    // Trace-back
    if (list.rbegin() != list.rend())
    {
        // Last vertex is end of heaviest increasing subsequence
        TVertexDescriptor v = list.rbegin()->second.second;
        while (true)
        {
            appendValue(pos, v, Generous());
            if (g.data_vertex[v])
                v = (*g.data_vertex[v]).data_target;
            else
                break;
        }
    }
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TSize2, typename TSpec2, typename TPositions, typename TString, typename TOutString>
inline void
_heaviest_common_subsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                             String<TSize2, TSpec2> const& slotToPos,
                             TPositions const& positions,
                             TString const& str1,
                             TString const& str2,
                             TOutString& align)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TString>::Type TVertexSet;
    typedef typename Iterator<TString const, Standard>::Type TStringIter;
    typedef typename Iterator<TString, Standard>::Type TSIter;
    typedef typename Iterator<TVertexSet const, Standard>::Type TVertexSetIter;

    TSize const m = length(str1);  // How many sets of vertex descriptors in seq1
    TSize const n = length(str2);  // How many sets of vertex descriptors in seq2
    TSize const numMatches = length(positions);
    TSize const alignLength = numMatches + (n - numMatches) + (m - numMatches);

    // Create the alignment sequence
    clear(align);
    resize(align, alignLength, TVertexSet(), Exact() );

    TSIter pointerAlign = begin(align, Standard());
    TSIter pointerAlignEnd = end(align, Standard());
    TSize posStr1 = 0;
    TSize posStr2 = 0;
    TStringIter pointerStr1 = begin(str1, Standard());
    TStringIter pointerStr2 = begin(str2, Standard());
    int p = length(positions) - 1;

    while (pointerAlign != pointerAlignEnd)
    {
        TSize i = m;
        TSize j = n;

        if (p >= 0)
        {
            i = (TSize) (slotToPos[positions[p]] / (TSize) n);   // Get the index in str1
            j = n - 1 - (TSize) (slotToPos[positions[p]] % (TSize) n); // Get the index in str2
        }

        // In what order do we insert gaps? -> Only important at the beginning and at the end, not between matches
        bool firstI = true;
        if ((i != posStr1) && (j != posStr2))
        {
            TSize tmpPosStr1 = posStr1;
            TSize tmpPosStr2 = posStr2;
            TSize len1 = 0;
            TSize len2 = 0;

            for (TStringIter tmpPointerStr1 = pointerStr1; i != tmpPosStr1; ++tmpPosStr1, ++tmpPointerStr1)
                len1 += fragmentLength(g, value(*tmpPointerStr1, 0));

            for (TStringIter tmpPointerStr2 = pointerStr2; j != tmpPosStr2; ++tmpPosStr2, ++tmpPointerStr2)
                len2 += fragmentLength(g, value(*tmpPointerStr2, 0));

            if (((posStr1 == 0) && (posStr2 == 0) && (len1 > len2)) || ((i == m) && (i == n) && (len1 < len2)))
            {
                firstI = false;
            }
        }

        auto append_to_align = [&pointerAlign] (auto & pointer_str, auto & pos)
        {
            for (TVertexSetIter itV = begin(*pointer_str, Standard()); itV != end(*pointer_str, Standard()); ++itV)
                appendValue(*pointerAlign, *itV, Generous());

            ++pointerAlign;
            ++pointer_str;
            ++pos;
        };

        if (firstI)
        {
            while (i != posStr1) // Gaps in seq 2
                append_to_align(pointerStr1, posStr1);

            while (j != posStr2) // Gaps in seq 1
                append_to_align(pointerStr2, posStr2);
        }
        else
        {
            while (j != posStr2) // Gaps in seq 1
                append_to_align(pointerStr2, posStr2);

            while (i != posStr1) // Gaps in seq 2
                append_to_align(pointerStr1, posStr1);
        }

        // Matches
        if (p >= 0)
        {
            for (TVertexSetIter itV = begin(*pointerStr1, Standard()); itV != end(*pointerStr1, Standard()); ++itV)
                appendValue(*pointerAlign, *itV, Generous());

            for (TVertexSetIter itV2 = begin(*pointerStr2, Standard()); itV2 != end(*pointerStr2, Standard()); ++itV2)
                appendValue(*pointerAlign, *itV2, Generous());

            ++pointerAlign;
            ++pointerStr1; ++posStr1;
            ++pointerStr2; ++posStr2;
            --p;
        }
    }
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
auto get_vertex_descriptor_pos_map(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                                   TString const& str1)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef String<TSize> TMapVertexPos;

    TMapVertexPos map;
    resize(map, getIdUpperBound(_getVertexIdManager(g)), std::numeric_limits<TSize>::max());

    TSize pos = 0;
    for (auto itStr1 = begin(str1, Standard()); itStr1 != end(str1, Standard()); ++itStr1, ++pos)
        for (auto itV = begin(*itStr1, Standard()); itV != end(*itStr1, Standard()); ++itV)
            map[*itV] = pos;

    return map;
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
auto get_slots_and_weights_and_seq(Graph<Alignment<TStringSet, TCargo, TSpec>> const & g,
                                   TString const & str1,
                                   TString const & str2)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Iterator<TString const, Standard>::Type TStringIterConst;
    typedef typename Value<TString>::Type TVertexSet;
    typedef typename Iterator<TVertexSet const, Standard>::Type TVertexSetIterConst;
    typedef String<TSize> TSlotToPos;
    typedef String<TCargo> TWeights;
    typedef String<TSize> TSequenceString;

    auto const map = get_vertex_descriptor_pos_map(g, str1);

    TSize const n = length(str2);  // How many sets of vertex descriptors in seq2

    // We could create the full graph -> too expensive
    // Remember which edges are actually present
    std::map<TSize, TCargo> slots_and_weights;

    TSize posItStr2 = 0;
    for (auto itStr2 = begin(str2, Standard()); itStr2 != end(str2, Standard()); ++itStr2, ++posItStr2)
    {
        for (auto itV = begin(*itStr2, Standard()); itV != end(*itStr2, Standard()); ++itV)
        {
            for (auto itOut = TOutEdgeIterator(g, *itV); !atEnd(itOut); ++itOut)
            {
                // Target vertex must be in the map
                TSize pPos = map[targetVertex(itOut)];
                if (pPos != std::numeric_limits<TSize>::max()) // if edge points to vertex of str1
                {
                    TSize const key = pPos * n + (TSize) (n - posItStr2 - 1);
                    slots_and_weights[key] += (TCargo) cargo(*itOut);
                }
            }
        }
    }

    // Get all occupied positions
    TSlotToPos slotToPos;
    TWeights weights;
    TSequenceString seq;
    reserve(slotToPos, slots_and_weights.size());
    reserve(weights, slots_and_weights.size());
    reserve(seq, slots_and_weights.size());

    for (auto [slot, weight] : slots_and_weights)
    {
        appendValue(slotToPos, slot);
        appendValue(weights, weight);
        appendValue(seq, n - 1 - (slot % n));
    }

    return std::make_tuple(std::move(slotToPos), std::move(weights), std::move(seq));
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TOutString>
inline void heaviest_common_subsequence(Graph<Alignment<TStringSet, TCargo, TSpec>> const & g,
                                        TString const & str1,
                                        TString const & str2,
                                        TOutString & align)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;

    // Note for profile alignments every member of the sequence is a String!!! of vertex descriptors

    // Fill the vertex to position map for str1
    // Remember for each vertex descriptor the position in the sequence
    auto const [slotToPos, weights, seq] = get_slots_and_weights_and_seq(g, str1, str2);

    // Now the tough part: Find the right number for a given position
    // Calculate the heaviest increasing subsequence
    String<TSize> positions;
    heaviest_increasing_subsequence(seq, weights, positions);

    // Retrieve the alignment sequence
    _heaviest_common_subsequence(g, slotToPos, positions, str1, str2, align);
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph>
inline void
progressive_alignment(Graph<Alignment<TStringSet, TCargo, TSpec>> const & g,
                      TGuideTree const & tree,
                      TOutGraph& gOut)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
    typedef typename Iterator<TGuideTree, BfsIterator>::Type TBfsIterator;
    typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
    typedef String<TVertexDescriptor> TVertexString;
    typedef String<TVertexString> TSegmentString;

    // Initialization
    TVertexDescriptor rootVertex = getRoot(tree);
    TSize nVertices = numVertices(tree);

    // Vertices in reversed bfs order
    TVertexString vertices;
    resize(vertices, nVertices);

    // All Strings of Strings of vertices for each node of the guide tree
    String<TSegmentString> segString;
    resize(segString, nVertices);

    // Walk through the tree in bfs order
    typedef typename Iterator<TVertexString, Standard>::Type TVertexIter;
    TVertexIter current_node = begin(vertices, Standard());
    TVertexIter nodes_end = end(vertices, Standard());
    --nodes_end;
    TBfsIterator bfsIt(tree, rootVertex);
    for (; !atEnd(bfsIt); goNext(bfsIt), --nodes_end)
        *nodes_end = *bfsIt;

    // Progressive alignment
    current_node = begin(vertices, Standard());
    nodes_end = end(vertices, Standard());
    for (; current_node != nodes_end; ++current_node)
    {
        if (isLeaf(tree, *current_node))
        {
            build_leaf_string(g, *current_node, segString[*current_node]);
        }
        else
        {
            // Align the two children (Binary tree)
            TAdjacencyIterator child2(tree, *current_node);
            TVertexDescriptor child1 = *child2;
            goNext(child2);
            heaviest_common_subsequence(g, segString[child1], segString[*child2], segString[*current_node]);
            clear(segString[child1]);
            clear(segString[*child2]);
        }
    }

    // Create the alignment graph
    _createAlignmentGraph(g, segString[rootVertex], gOut);
}

}
