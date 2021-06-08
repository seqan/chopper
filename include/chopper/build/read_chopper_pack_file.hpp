#pragma once

#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <seqan3/std/ranges>

#include <chopper/build/build_data.hpp>
#include <chopper/detail_parse_chopper_pack_header_line.hpp>
#include <chopper/detail_parse_chopper_pack_line.hpp>

// data needs to be passed from outside sind the graph in data cannot be moved
inline void read_chopper_pack_file(build_data<chopper_pack_record> & data, std::string const & chopper_pack_filename)
{
    std::ifstream chopper_pack_file{chopper_pack_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_pack_filename + " for reading"};

    // parse header
    // -------------------------------------------------------------------------
    data.number_of_ibfs = parse_chopper_pack_header(data.ibf_graph, data.node_map, chopper_pack_file) + 1;

    // parse lines
    // -------------------------------------------------------------------------
    std::string current_line;
    size_t user_bins{};
    while (std::getline(chopper_pack_file, current_line))
    {
        ++user_bins;
        chopper_pack_record const && record = parse_chopper_pack_line(current_line);

        // go down the tree until you find the matching parent
        lemon::ListDigraph::Node current_node = data.ibf_graph.nodeFromId(0); // start at root

        for (size_t i = 0; i < record.bin_indices.size() - 1; ++i)
        {
            size_t const bin = record.bin_indices[i];
            size_t const num_tbs = record.number_of_bins[i];
            auto & current_data = data.node_map[current_node];

            // update number of technical bins in current_node-IBF
            current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

            bool found_next_node{false}; // sanity check
            for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            {
                auto target = data.ibf_graph.target(arc_it);
                if (data.node_map[target].parent_bin_index == bin)
                {
                    current_node = target;
                    found_next_node = true;
                    break;
                }
            }
            assert(found_next_node);
        }

        size_t const bin = record.bin_indices.back();
        size_t const num_tbs = record.number_of_bins.back();
        auto & current_data = data.node_map[current_node];

        // update number of technical bins in current_node-IBF
        current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + num_tbs);

        if (record.bin_indices.back() == current_data.max_bin_index)
            current_data.remaining_records.insert(current_data.remaining_records.begin(), record);
        else
            current_data.remaining_records.push_back(record);
    }

    data.number_of_user_bins = user_bins;
    data.resize();
};
