#pragma once

#include <vector>

#include <lemon/list_graph.h>

#include <chopper/detail_parse_chopper_pack_line.hpp> // for chopper_pack_record

template <typename record_type>
struct node_data // rename:ibf_data? or ibf_node_data
{
    size_t parent_bin_index{};
    size_t max_bin_index{};
    size_t number_of_technical_bins{};
    lemon::ListDigraph::Node favourite_child{lemon::INVALID};
    std::vector<record_type> remaining_records{}; // non-merged bins (either split or single)
};

std::ostream & operator<<(std::ostream & s, node_data<chopper_pack_record> const & nd)
{
    s << "{parent:"
      << nd.parent_bin_index
      << ", max:"
      << nd.max_bin_index
      << ", num:"
      << nd.number_of_technical_bins
      << ", fav:"
      << ((nd.favourite_child == lemon::INVALID) ? "NONE" : "some")
      << ", recs:<";
    for (auto const & rec : nd.remaining_records)
    {
        for (auto const & f : rec.filenames)
            s << f << ";";
        s << "|";
        for (auto const & f : rec.bin_indices)
            s << f << ";";
        s << "|";
        for (auto const & f : rec.number_of_bins)
            s << f << ";";
    }
    s  << ">}";

    return s;
}
