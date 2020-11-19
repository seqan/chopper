#include <iostream>
#include <fstream>
#include <set>
#include <queue>

#include <lemon/color.h>
#include <lemon/concepts/digraph.h>
#include <lemon/core.h>
#include <lemon/dimacs.h>
#include <lemon/list_graph.h>
#include <lemon/static_graph.h>

#define SEQAN3_HAS_ZLIB 1

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/charconv>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include <chopper/split/split_config.hpp>

namespace lemon {

const Color WHITE(1,1,1);

const Color BLACK(0,0,0);
const Color RED(1,0,0);
const Color GREEN(0,1,0);
const Color BLUE(0,0,1);
const Color YELLOW(1,1,0);
const Color MAGENTA(1,0,1);
const Color CYAN(0,1,1);

const Color GREY(0,0,0);
const Color DARK_RED(.5,0,0);
const Color DARK_GREEN(0,.5,0);
const Color DARK_BLUE(0,0,.5);
const Color DARK_YELLOW(.5,.5,0);
const Color DARK_MAGENTA(.5,0,.5);
const Color DARK_CYAN(0,.5,.5);

} //namespace lemon


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
                std::filesystem::path const & graph_file_name)
{
    std::ifstream stream{graph_file_name};

    auto file_view = std::ranges::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                           decltype(std::istreambuf_iterator<char>{})>
                        {std::istreambuf_iterator<char>{stream},
                         std::istreambuf_iterator<char>{}};

    if (!std::ranges::equal(file_view | seqan3::views::take_line, std::string{"graph G {"}))
            throw std::runtime_error{"not in dot format: graph G { not found in the beginning."};

    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>));
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>)); /* Graph Attributes */
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>));
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>)); /* Node Attributes */
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>));
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>)); /* Edge Attributes */
    seqan3::detail::consume(file_view | seqan3::views::take_until_or_throw_and_consume(seqan3::is_char<'/'>));

    if (!std::ranges::equal(file_view | seqan3::views::take_line, std::string{"* Sequence Lengths */"}))
            throw std::runtime_error{"no sequence lengths :( you have on old tcofee-seqan version."};

    std::vector<uint32_t> seq_lengths{};

    while (!seqan3::is_char<'/'>(*std::ranges::begin(file_view)))
    {
        char buffer[1000];
        auto end_ptr = std::ranges::copy(file_view | seqan3::views::take_line, &buffer[0]).out;
        *end_ptr = '\0'; // ensure that string ends here to avoid parsing errors if the previous line was longer
        char * ptr = &buffer[0];

        uint32_t len{};

        while (!seqan3::is_char<'\t'>(*ptr)) ++ptr; // skip name for now
        std::from_chars(ptr + 1, &buffer[1000], len);

        seq_lengths.push_back(len);
    }

    if (seq_lengths.size() > 65536u)
        throw std::logic_error{"Currently node ids are uint16_t, not all ids can be represented!"};

    // seqan3::debug_stream << seq_lengths << std::endl;

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
        if (range[i].second == seq_lengths[group[i]]) // [unlikely]
        {
            // seqan3::debug_stream << "i:" << i << " group[i]:" << group[i] << " seq_lengths[group[i]]:" << seq_lengths[group[i]] << std::endl;
            g.addArc(node, sink_node);
        }

        // fill node map of 'node'
        std::vector<std::pair<uint32_t, uint32_t>> node_property(unique_groups.size());
        assert(group[i] < unique_groups.size());
        node_property[group[i]] = range[i];
        node_map.set(node, node_property);
    }

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

void write_graph(lemon::ListDigraph const & g,
                 lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> const & node_map,
                 std::filesystem::path const & graph_file_name)
{
    std::ofstream sout{graph_file_name};
    seqan3::debug_stream_type fout{sout};
    std::cerr << "[LOG] Writing graph to " << graph_file_name << std::endl;

    std::string header = R"(graph G {
/* Graph Attributes */
graph [rankdir = LR];

/* Node Attributes */
node [shape = rectangle, fillcolor = white, style = filled, fontname = "Times-Italic"];

/* Edge Attributes */
edge [fontname = "Times-Italic", arrowsize = 0.75, fontsize = 16];

/* Nodes */
)";

    fout << header;

    for (lemon::ListDigraph::NodeIt u(g); u != lemon::INVALID; ++u)
        fout << g.id(u) << " [label = \"" << g.id(u) << ":"<< node_map[u] << "\", group = " << 0 << "];" << std::endl;

    fout << std::endl << "/*directed edges*/" << std::endl;

    for (lemon::ListDigraph::ArcIt a(g); a != lemon::INVALID; ++a)
        fout << g.id(g.source(a)) << " -- " << g.id(g.target(a)) << " [];" << std::endl;

    fout << std::endl << "}" << std::endl;
}

void traverse_graph(split_data const & data, split_config const & config)
{
    lemon::ListDigraph g;
    std::vector<lemon::ListDigraph::Node> nodes;
    lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> node_map{g};

    read_graph(g, nodes, node_map, config.output_graph_file); // also merges nodes along undirected edges

    seqan3::debug_stream << "[LOG] " <<  lemon::countNodes(g) << " nodes remain after merging." << std::endl;
    seqan3::debug_stream << "[LOG] " <<  lemon::countArcs(g) << " arcs remain after merging." << std::endl;

    if (config.write_out_graph)
    {
        write_graph(g, node_map, "/tmp/graph_out.dot");
        seqan3::debug_stream << "[LOG] Written graph to /tmp/graph_out.dot" << std::endl;
    }

    std::ofstream weights_file;
    if (config.write_out_weights)
        weights_file.open("/tmp/weights.tsv");

    // sanity check for graph:
    if (lemon::countInArcs(g, nodes[0]) != 0)
        throw std::logic_error{"The source node cannot have any incoming arcs."};
    if (lemon::countOutArcs(g, nodes[0]) != node_map[nodes[0]].size())
        throw std::logic_error{"The source node must have as many outgoing arcs as there are sequences."};
    if (lemon::countOutArcs(g, nodes[1]) != 0)
        throw std::logic_error{"The sink node cannot have any outgoing arcs."};
    if (lemon::countInArcs(g, nodes[1]) != node_map[nodes[1]].size())
        throw std::logic_error{"The sink node must have as many incoming arcs as there are sequences."};

    // When assigning a node to the bin, we want to approximate the kmer_content it adds.
    // This is done by first assigning the node weights, according to the approximate sequence content they span.
    // All weights are summed up divided by the number of bins is an greedy estimate to split bin
    lemon::ListDigraph::NodeMap<uint64_t> node_weigths{g};
    uint64_t weight_per_bin{};
    for (lemon::ListDigraph::NodeIt node_it(g); node_it != lemon::INVALID; ++node_it)
    {
        uint64_t max{};
        for (auto & p : node_map[node_it])
        {
            assert(p.second >= p.first);
            if (p.second - p.first > max)
                max = p.second - p.first;
        }
        node_weigths.set(node_it, max);
        weight_per_bin += max;
    }
    seqan3::debug_stream << "[LOG] total weight:" << weight_per_bin << std::endl;
    weight_per_bin = weight_per_bin / config.bins;
    seqan3::debug_stream << "[LOG] weight per bin :" << weight_per_bin << std::endl;

    // int32_t nodes_per_bin{lemon::countNodes(g) / config.bins}; // #nodes == #minimizers
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> bin_contents{};

    // Since the graph is an alignment graph it is linear.
    // We traverse the graph from left to right adding nodes greedily to a bin.
    // We start at the node with the most incoming edges (shared by most sequences). [todo peek ahead to the best 10]
    // From this, we recursively always chose the next node with the most connections to the current node [todo peak ahead?]
    // and add him to the bin, as well as all nodes to the left connected to the current bin:
    //
    //    -- [node 3] --        If current node is node 1 and we want to add node 2, we also add node 3.
    //    |            |
    // [node 1] == [node 2]
    //
    // We track the number of nodes we added in each iteration and stop with the current bin, if sum of weights > #weight_per_bin
    uint64_t count{};
    std::map<lemon::ListDigraph::Node, uint16_t> arc_count_per_node{};

    lemon::ListDigraph::Node current_node = nodes[0]; // start at source node
    std::vector<std::pair<uint32_t, uint32_t>> current_ranges = node_map[current_node];

    while (current_node != nodes[1]) // did not arrive in sink node yet
    {
        arc_count_per_node.clear();
        for (lemon::ListDigraph::OutArcIt arc_it(g, current_node); arc_it != lemon::INVALID; ++arc_it)
            ++arc_count_per_node[g.target(arc_it)];

        auto next_node = (*std::max_element(arc_count_per_node.begin(),
                                            arc_count_per_node.end(),
                                            [] (auto p1, auto p2) { return p1.second < p2.second; })).first;

        // seqan3::debug_stream << " --> Current Node " << g.id(current_node)
        //                      << " weight: " << node_weigths[current_node]
        //                      << " map:" << node_map[current_node]
        //                      << "\t current range:" << current_ranges
        //                      << " count:" << count
        //                      << "\t Next Node " << g.id(next_node)
        //                      << " weight: " << node_weigths[next_node]
        //                      << " map:" << node_map[next_node] << std::endl;

        g.erase(current_node);

        // remove any loops... although I don't now why there are loops in the first place..
        for (lemon::ListDigraph::InArcIt arc_it(g, next_node); arc_it != lemon::INVALID; ++arc_it)
        {
            if (g.source(arc_it) == next_node)
            {
                g.erase(arc_it);
                break;
            }
        }

        auto incoming_arc = lemon::ListDigraph::InArcIt(g, next_node);
        if (incoming_arc == lemon::INVALID) // no more incoming arc in next_node
        {
            // no op
            // seqan3::debug_stream << " -- --> Node " << g.id(next_node) << " has no incoming arcs -> can be consumed.\n";
        }
        else
        {
            // chose a random out arc, now the first. maybe it is better to chose the node with the most connecting nodes
            // seqan3::debug_stream << " -- --> recurse... ";
            while (incoming_arc != lemon::INVALID)
            {
                next_node = g.source(incoming_arc);
                // seqan3::debug_stream << " -> " << g.id(next_node);
                incoming_arc = lemon::ListDigraph::InArcIt(g, next_node);
            }
            // seqan3::debug_stream << " -- --> recursively chosen node: " << g.id(next_node);

            // seqan3::debug_stream << "NEW NEXT NODE " << g.id(next_node)
            //                      << " weight: " << node_weigths[next_node]
            //                      << " map:" << node_map[next_node] << std::endl;
        }


        if (count + node_weigths[next_node] > weight_per_bin) // then flush current results into bin contents
        {
            // seqan3::debug_stream << std::endl << " =========================================================== count:\t" << count << std::endl;

            if (count + node_weigths[next_node] > weight_per_bin * 1.05) // big fat node coming -> split node
            {
                auto surplus = count + node_weigths[next_node] - weight_per_bin - 1;

                if (config.write_out_weights)
                    weights_file << (count + node_weigths[next_node] - surplus) << "\t";

                for (auto && [curr, next] : seqan3::views::zip(current_ranges, node_map[next_node]))
                {
                    if (next.first != 0 || next.second != 0) // TODO is next.second != 0 sufficient?
                    {
                        uint32_t new_end = std::max<uint32_t>(next.first, ((next.second >= surplus) ? next.second - surplus : 0u));
                        curr.second = new_end;
                        next.first = new_end;
                    }
                }

                // seqan3::debug_stream << " surplus: " << surplus
                //                      << " new_ranges: " << split_ranges
                //                      << " current_ranges:" << current_ranges
                //                      << std::endl;

                bin_contents.push_back(current_ranges);
                for (auto & [begin, end] : current_ranges)
                    begin = end;

                count = surplus;

                for (auto && [curr, next] : seqan3::views::zip(current_ranges, node_map[next_node]))
                {
                    if (curr.first == 0 && curr.second == 0)
                    {
                        curr.first = next.first;
                        curr.second = next.second;
                    }
                    else
                    {
                        if (next.first != 0 || next.second != 0) // TODO is next.second != 0 sufficient?
                        {
                            curr.second = next.second;
                        }
                    }
                }
            }
            else
            {
                if (config.write_out_weights)
                    weights_file << (count  + node_weigths[next_node]) << "\t";

                for (auto && [curr, next] : seqan3::views::zip(current_ranges, node_map[next_node]))
                {
                    if (curr.first == 0 && curr.second == 0)
                    {
                        curr.first = next.first;
                        curr.second = next.second;
                    }
                    else
                    {
                        if (next.first != 0 || next.second != 0) // TODO is next.second != 0 sufficient?
                        {
                            curr.second = next.second;
                        }
                    }
                }

                bin_contents.push_back(current_ranges);
                // seqan3::debug_stream << "count: " << count << " current map: " << node_map[next_node] << std::endl;
                for (auto & [begin, end] : current_ranges)
                    begin = end;

                count = 0;
                // seqan3::debug_stream << " current_ranges:" << current_ranges
                //                      << std::endl;
            }
        }
        else
        {
            count += node_weigths[next_node];

            for (auto && [curr, next] : seqan3::views::zip(current_ranges, node_map[next_node]))
            {
                if (curr.first == 0 && curr.second == 0)
                {
                    curr.first = next.first;
                    curr.second = next.second;
                }
                else
                {
                    if (next.first != 0 || next.second != 0) // TODO is next.second != 0 sufficient?
                    {
                        curr.second = next.second;
                    }
                }
            }
        }

        current_node = next_node;
    }

    bin_contents.push_back(current_ranges); // push_back last bin content

    seqan3::debug_stream << "[LOG] " <<  lemon::countNodes(g) << " nodes remain. (should be 1 only the sink)" << std::endl;

    assert(bin_contents.size() != 0);
    assert(bin_contents[0].size() == seqan::length(data.sequences));

    // remember if this is the first time this file is opened, s.t. we write the header line
    bool const output_file_exists = std::filesystem::exists(config.out_path);

    std::ofstream fout{config.out_path, std::ios::binary | std::ios::app}; // append to file

    if (!output_file_exists)
        fout << "FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER\n"; // header

    std::vector<bool> range_end_sanity_check(seqan::length(data.sequences), 0);
    size_t bin_index{};
    for (auto & bin : bin_contents)
    {
        for (size_t i = 0; i < seqan::length(data.sequences); ++i)
        {
            if (bin[i].first != bin[i].second)
            {
                fout << data.files_of_origin[i] << '\t'
                     << data.ids[i] << '\t'
                     << bin[i].first << '\t'
                     << bin[i].second << '\t'
                     << bin_index + config.bin_index_offset << '\n';
            }

            if (bin[i].second == data.lengths[i])
                range_end_sanity_check[i] = true;
        }
        ++bin_index;
    }

    bool sanity_check_passed{true};
    for (bool end_matched : range_end_sanity_check)
        sanity_check_passed = sanity_check_passed && end_matched;
    if (!sanity_check_passed)
    {
        seqan3::debug_stream << "[WARNING] Sanity check at the end did not pass. "
                             << "Not all sequence ends were contained in the traversed ranges. "
                             << std::endl;
        for (size_t i = 0; i < range_end_sanity_check.size(); ++i)
            if (!range_end_sanity_check[i])
                std::cerr << "[WARNING] End did not match for " << data.ids[i] << " of length " << data.lengths[i] << std::endl;
        seqan3::debug_stream << "[WARNING] This should never happen. Please contact the developer." << std::endl;
    }

    seqan3::debug_stream << "Bin content "
                         << (output_file_exists ? "appended" : "written")
                         << " to : "
                         << config.out_path << std::endl;
}
