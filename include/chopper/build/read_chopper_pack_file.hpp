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

auto read_chopper_pack_file(std::string const & chopper_pack_filename)
{
    build_data data;
    std::vector<std::vector<chopper_pack_record>> records_per_hibf_bin{};

    std::ifstream chopper_pack_file{chopper_pack_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_pack_filename + " for reading"};

    // parse header
    // -------------------------------------------------------------------------
    lemon::ListDigraph & ibf_initialiser_graph{};
    parse_chopper_pack_header(ibf_initialiser_graph, chopper_pack_file);

    // parse lines
    // -------------------------------------------------------------------------
    do
    {
        chopper_pack_record const && record = parse_chopper_pack_line(current_line);
        data.num_libfs += (record.lidx != -1);

        if (records_per_hibf_bin.size() <= record.hidx)
            records_per_hibf_bin.resize(record.hidx + ((record.lidx == -1) ? record.bins : 1));

        records_per_hibf_bin[record.hidx].push_back(record);

    } while (std::getline(chopper_pack_file, current_line));

    for (auto & records : records_per_hibf_bin)
        std::ranges::sort(records, [] (auto const & rec1, auto const & rec2) { return rec1.lidx < rec2.lidx; });

    data.hibf_num_technical_bins = records_per_hibf_bin.size();

    return std::make_pair(std::move(data), std::move(records_per_hibf_bin));
};
