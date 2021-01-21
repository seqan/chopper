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
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include <chopper/split/split_config.hpp>
#include <chopper/split/split_data.hpp>

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

void traverse_graph(lemon::ListDigraph & g,
                    std::vector<lemon::ListDigraph::Node> const & nodes,
                    lemon::ListDigraph::NodeMap<std::vector<std::pair<uint32_t, uint32_t>>> & node_map,
                    split_data const & data,
                    batch_config const & config)
{
    if (config.verbose)
    {
        seqan3::debug_stream << "[LOG] " <<  lemon::countNodes(g) << " nodes remain after merging." << std::endl;
        seqan3::debug_stream << "[LOG] " <<  lemon::countArcs(g) << " arcs remain after merging." << std::endl;
    }

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

    if (config.verbose)
        seqan3::debug_stream << "[LOG] total weight:" << weight_per_bin << std::endl;

    weight_per_bin = weight_per_bin / config.bins;

    if (config.verbose)
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
                count = count + node_weigths[next_node];

                while (count > weight_per_bin * 1.05)
                {
                    auto surplus = count - weight_per_bin - 1;

                    if (config.write_out_weights)
                        weights_file << (count - surplus) << "\t";

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
                    //                      << " new_ranges: " << node_map[next_node]
                    //                      << " current_ranges:" << current_ranges
                    //                      << std::endl;

                    bin_contents.push_back(current_ranges);
                    for (auto & [begin, end] : current_ranges)
                        begin = end;

                    count = surplus;
                }

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

    if (config.verbose)
        seqan3::debug_stream << "[LOG] " <<  lemon::countNodes(g) << " nodes remain. (should be 1 only the sink)" << std::endl;

    assert(bin_contents.size() != 0);
    assert(bin_contents[0].size() == seqan::length(data.sequences));

    std::vector<bool> range_end_sanity_check(seqan::length(data.sequences), 0);
    size_t bin_index{};
    for (auto & bin : bin_contents)
    {
        for (size_t i = 0; i < seqan::length(data.sequences); ++i)
        {
            if (bin[i].first != bin[i].second)
            {
                *data.outstream << data.files_of_origin[i] << '\t'
                                << data.ids[i] << '\t'
                                << bin[i].first << '\t'
                                << bin[i].second << '\t'
                                << ((config.merged_bin) ? 0 : bin_index) + config.hibf_bin_idx_offset << '\t'
                                << ((config.merged_bin) ? std::to_string(bin_index + config.libf_bin_idx_offset) : "-")
                                << '\n';
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
}
