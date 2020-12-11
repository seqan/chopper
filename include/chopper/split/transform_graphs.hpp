#pragma once

#include <seqan3/std/charconv>

#include <lemon/color.h>
#include <lemon/concepts/digraph.h>
#include <lemon/core.h>
#include <lemon/dimacs.h>
#include <lemon/list_graph.h>
#include <lemon/static_graph.h>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/take_line.hpp>

template <typename graph_type>
void seqan2_write_graph(graph_type const & gAlign, split_data const & data, batch_config const & config)
{
    seqan::String<seqan::String<char> > edgeMap;
    seqan::_createEdgeAttributes(gAlign, edgeMap);

    // create node attributes
    typedef typename seqan::Id<graph_type>::Type TIdType;

    std::ofstream dotFile(config.output_graph_file);

    typedef typename seqan::VertexDescriptor<graph_type>::Type TVertexDescriptor;
    typename seqan::DirectionIterator<std::ofstream, seqan::Output>::Type iter = directionIterator(dotFile, seqan::Output());

    seqan::write(iter, "/* Nodes */\n");
    typedef typename seqan::Iterator<graph_type, seqan::VertexIterator>::Type TConstIter;

    for(TConstIter it(gAlign);!seqan::atEnd(it);++it) {
        seqan::appendNumber(iter, (int)*it);
        seqan::write(iter, " [");

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

        seqan::write(iter, outs.str());
        seqan::write(iter, "];\n");
    }
    seqan::writeValue(iter, '\n');

    seqan::write(iter, "/*directed edges*/\n");
    seqan::_writeGraphFooter(iter, gAlign, seqan::DotDrawing()); // somehow this writes directed edges :shrug:

    seqan::write(iter, "/* Edges */\n");
    typedef typename seqan::Iterator<graph_type, seqan::EdgeIterator>::Type TConstEdIter;
    TConstEdIter itEd(gAlign);
    for(;!seqan::atEnd(itEd);++itEd) {
        TVertexDescriptor sc = seqan::sourceVertex(itEd);
        TVertexDescriptor tr = seqan::targetVertex(itEd);
        seqan::appendNumber(iter, sc);
        seqan::_writeEdgeType(iter, gAlign, seqan::DotDrawing());
        seqan::appendNumber(iter, tr);
        seqan::write(iter, " [");
        seqan::write(iter, seqan::getProperty(edgeMap, *itEd));
        seqan::write(iter, "];\n");
    }
    seqan::writeValue(iter, '\n');

    seqan::write(iter, "}\n");

    dotFile.close();
}

/* [source] --> [target]            [target]
 *  (2,10)       (10,15)     =>      (2,15)
 *
 * [source] --> [target]            [target]
 *   (0,0)       (10,15)     =>      (10,15)
 *
 * [source] --> [target]            [target]
 *  (2,10)       (0,0)       =>      (2,10)
 */
template <typename node_map_type>
void merge_properties_into(node_map_type            & node_map,
                           typename node_map_type::Key const & source,
                           typename node_map_type::Key const & target)
{
    typename node_map_type::Value & target_value = node_map[target];

    uint16_t idx{};
    for (auto const & [range_start, range_end] : node_map[source])
    {
        assert(range_end >= range_start);
        if (range_start != 0 || range_end != 0) // source non empty
        {
            if (target_value[idx].first == 0 && target_value[idx].second == 0) // target empty
            {
                target_value[idx].second = range_end;
            }
            else
            {
                if (range_end != target_value[idx].first)
                    throw std::runtime_error{"Possibly corrupted input graph (Range idx " + std::to_string(idx) +
                                             "). Source range: [" + std::to_string(range_start) +
                                             ", " + std::to_string(range_end) +
                                             "]Error: source range end was not equal the target range start: " +
                                             std::to_string(target_value[idx].first)};
            }

            target_value[idx].first = range_start;
        }
        ++idx;
    }
}

void read_graph(lemon::ListDigraph & g,
                std::vector<lemon::ListDigraph::Node> & nodes,
                lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> & node_map,
                std::filesystem::path const & graph_file_name,
                split_data const & data,
                batch_config const & config)
{
    std::ifstream stream{graph_file_name};

    auto file_view = std::ranges::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                           decltype(std::istreambuf_iterator<char>{})>
                        {std::istreambuf_iterator<char>{stream},
                         std::istreambuf_iterator<char>{}};

    if (seqan::length(data.lengths) > 65536u)
        throw std::logic_error{"Currently node ids are uint16_t, not all ids can be represented!"};

    // seqan3::debug_stream << data.lengths << std::endl;

    if (!std::ranges::equal(file_view | seqan3::views::take_line, std::string{"/* Nodes */"}))
            throw std::runtime_error{"not in dot format: No nodes :("};

    // read in nodes
    std::set<uint16_t> unique_groups{};
    std::vector<uint16_t> group{};
    std::vector<std::pair<uint32_t, uint32_t>> range{};

    // read in nodes
    while (!seqan3::is_char<'/'>(*std::ranges::begin(file_view)))
    {
        char buffer[100];
        std::ranges::copy(file_view | seqan3::views::take_line, &buffer[0]);
        char * ptr = &buffer[0];

        uint32_t range_start{};
        uint32_t range_end{};
        uint16_t gr{};

        while (!seqan3::is_char<'"'>(*ptr)) ++ptr;
        auto res = std::from_chars(ptr + 2, &buffer[100], range_start);

        if (std::string(res.ptr + 1, res.ptr + 4) == "end")
            throw std::logic_error{"You have an outdated graph. Please run seqan_tcoffee again."};

        res = std::from_chars(res.ptr + 1, &buffer[100], range_end);
        while (!seqan3::is_char<'='>(*ptr)) ++ptr;
        std::from_chars(ptr + 2, &buffer[100], gr);

        if (range_end < range_start)
        {
            seqan3::debug_stream << "wrong range: [" << range_start << "," << range_end << "]" << std::endl;
        }

        range.emplace_back(std::move(range_start), std::move(range_end));

        unique_groups.insert(gr);
        group.emplace_back(std::move(gr));
    }
    // Note: we assume that the order of the node read before is sorted by id!

    nodes.reserve(group.size() + 1);
    g.reserveNode(group.size() + 1);
    g.reserveArc(group.size() + 1);

    lemon::ListDigraph::Node source_node = g.addNode();
    lemon::ListDigraph::Node sink_node = g.addNode();
    nodes.push_back(source_node);
    nodes.push_back(sink_node);
    node_map.set(source_node, std::vector<std::pair<uint32_t, uint32_t>>(unique_groups.size()));
    node_map.set(sink_node, std::vector<std::pair<uint32_t, uint32_t>>(unique_groups.size()));

    for (size_t i = 0; i < group.size(); ++i)
    {
        auto node = g.addNode();
        nodes.push_back(node);

        // add arc from source node to node if it is a start node (start of sequence range)
        if (range[i].first == 0) // [unlikely]
            g.addArc(source_node, node);
        if (range[i].second == data.lengths[group[i]]) // [unlikely]
        {
            // seqan3::debug_stream << "i:" << i << " group[i]:" << group[i] << " data.lengths[group[i]]:" << data.lengths[group[i]] << std::endl;
            g.addArc(node, sink_node);
        }

        // fill node map of 'node'
        std::vector<std::pair<uint32_t, uint32_t>> node_property(unique_groups.size());
        assert(group[i] < unique_groups.size());
        node_property[group[i]] = range[i];
        node_map.set(node, node_property);
    }

    if (config.verbose)
        seqan3::debug_stream << "[LOG] inserted " <<  lemon::countNodes(g) << " nodes into the graph." << std::endl;

    // read directed edges:
    if (!std::ranges::equal(file_view | seqan3::views::take_line, std::string{"/*directed edges*/"}))
        throw std::runtime_error{"please reorder you dot file in case undirected edges were first."};

    while (!seqan3::is_char<'/'>(*std::ranges::begin(file_view)))
    {
        char buffer[100];
        std::ranges::copy(file_view | seqan3::views::take_line, &buffer[0]);

        uint32_t source{};
        uint32_t target{};

        auto res = std::from_chars(&buffer[0], &buffer[100], source);
        while (!seqan3::is_char<'-'>(*res.ptr)) ++res.ptr;
        std::from_chars(res.ptr + 3, &buffer[100], target);

        g.addArc(nodes[source + 2], nodes[target + 2]); // plus 2 because ids were shifted when we inserted a source/sink node
    }

    if (config.verbose)
        seqan3::debug_stream << "[LOG] inserted " <<  lemon::countArcs(g) << " arcs into the graph." << std::endl;

    // read undirected edges:
    if (!std::ranges::equal(file_view | seqan3::views::take_line, std::string{"/* Edges */"}))
        throw std::runtime_error{"No edges :(."};

    while (!seqan3::is_char<'}'>(*std::ranges::begin(file_view)))
    {
        char buffer[100];
        std::ranges::copy(file_view | seqan3::views::take_line, &buffer[0]);

        uint32_t source{};
        uint32_t target{};

        auto res = std::from_chars(&buffer[0], &buffer[100], source);
        while (!seqan3::is_char<'-'>(*res.ptr)) ++res.ptr;
        std::from_chars(res.ptr + 3, &buffer[100], target);

        // plus 2 because ids were shifted when we inserted a source/sink node
        source += 2;
        target += 2;

        // the undirected edges connect minimizer
        if (nodes[source] != nodes[target])                   // if not combined already
        {
            merge_properties_into(node_map, nodes[target], nodes[source]);

            g.contract(nodes[source], nodes[target], true);   // merge target into source (true to remove any loops)
            nodes[target] = nodes[source];                    // merge target into source
        }
    }
}

// currently through an output file
template <typename graph_type>
void transform_graphs(lemon::ListDigraph & g,
                      std::vector<lemon::ListDigraph::Node> & nodes,
                      lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> & node_map,
                      graph_type const & seqan2_graph,
                      split_data const & data,
                      batch_config const & config)
{
    seqan2_write_graph(seqan2_graph, data, config);

    read_graph(g, nodes, node_map, config.output_graph_file, data, config); // also merges nodes along undirected edges
}
