#pragma once

#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/std/charconv>

#include <chopper/build/batch.hpp>
#include <chopper/build/build_data.hpp>
#include <chopper/build/region.hpp>
#include <chopper/build/chopper_split_record.hpp>
#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_parse_chopper_pack_header_line.hpp>
#include <chopper/detail_starts_with.hpp>

auto parse_chopper_split_line(std::string const & line)
{
    char const * buffer = line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const filename(field_start, field_end);

    assert(*field_end == '\t');
    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const id(field_start, field_end);

    region reg{};
    std::vector<size_t> bin_indices{};

    // read begin
    assert(*field_end == '\t');
    ++field_end;
    auto res = std::from_chars(field_end, buffer_end, reg.begin);
    field_end = res.ptr;

    // read end
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, reg.end);
    field_end = res.ptr;

    // read bin indices
    assert(*field_end == '\t');
    ++field_end;
    while (field_end < buffer_end && *field_end != ';' && *field_end != '\n')
    {
        res = std::from_chars(field_end, buffer_end, reg.bin_index);
        field_end = res.ptr;
        ++field_end;
        bin_indices.push_back(reg.bin_index);
    }

    return std::make_tuple(std::move(filename), std::move(id), std::move(bin_indices), std::move(reg));
}

// data needs to be passed from outside since the graph in data cannot be moved
void read_chopper_split_file(build_data<chopper_split_record> & data, std::string const & chopper_split_filename)
{
    std::ifstream chopper_pack_file{chopper_split_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_split_filename + " for reading"};

    // parse header
    // -------------------------------------------------------------------------
    // header in split file was copied over from pack file!
    parse_chopper_pack_header(data.ibf_graph, data.node_map, chopper_pack_file);

    // assign to ibf tree structure
    // -------------------------------------------------------------------------
    std::unordered_map<std::string, chopper_split_record> record_map;
    std::string current_line;
    while (std::getline(chopper_pack_file, current_line))
    {
        auto const && [filename, seq_id, bin_indices, region_record] = parse_chopper_split_line(current_line);

        // go down the tree until you find the matching parent
        lemon::ListDigraph::Node current_node = data.ibf_graph.nodeFromId(0); // start at root
        for (size_t i = 0; i < bin_indices.size() - 1; ++i)
        {
            size_t const bin = bin_indices[i];
            auto & current_data = data.node_map[current_node];

            // update number of technical bins in current_node-IBF
            current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + 1);

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

        size_t const bin = bin_indices.back();
        auto & current_data = data.node_map[current_node];

        bool found_record{false};
        for (auto & record : current_data.remaining_records)
        {
            auto filenames_it = std::find(record.filenames.begin(), record.filenames.end(), filename);
            auto bins_it = std::find(record.bin_indices.begin(), record.bin_indices.end(), bin);

            if (bins_it != record.bin_indices.end()) // there is a record with that bin already
            {
                if (filenames_it == record.filenames.end())
                    record.filenames.push_back(filename);

                record.region_map[filename + seq_id].push_back(region_record);
                found_record = true;
            }
            else if (filenames_it != record.filenames.end()) // there is a record that already contains this file
            {
                if (bin == current_data.max_bin_index)
                    record.bin_indices.insert(record.bin_indices.begin(), bin);
                else
                    record.bin_indices.push_back(bin);

                record.region_map[filename + seq_id].push_back(region_record);
                found_record = true;
            }
        }

        if (!found_record) // then insert new one
        {
            chopper_split_record record{{filename}, {bin}, {1}, {{filename + seq_id, {region_record}}}};

            if (bin == current_data.max_bin_index)
            {
                current_data.remaining_records.insert(current_data.remaining_records.begin(), record);
            }
            else
            {
                current_data.remaining_records.push_back(record);
            }
        }

        // update number of technical bins in current_node-IBF
        current_data.number_of_technical_bins = std::max(current_data.number_of_technical_bins, bin + 1);
    }
};
