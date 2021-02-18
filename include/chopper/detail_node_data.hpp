#pragma once

#include <vector>

#include <lemon/list_graph.h>

#include <chopper/detail_parse_chopper_pack_line.hpp> // for chopper_pack_record
#include <chopper/build/chopper_split_record.hpp>

template <typename record_type>
struct node_data // rename:ibf_data? or ibf_node_data
{
    size_t parent_bin_index{};
    size_t max_bin_index{};
    size_t number_of_technical_bins{};
    lemon::ListDigraph::Node favourite_child{lemon::INVALID};
    std::vector<record_type> remaining_records{}; // non-merged bins (either split or single)

    bool operator==(node_data const & rhs) const
    {
        bool res = std::tie(parent_bin_index, max_bin_index, number_of_technical_bins, favourite_child) ==
                   std::tie(rhs.parent_bin_index, rhs.max_bin_index, rhs.number_of_technical_bins, rhs.favourite_child);

        if (remaining_records.size() != rhs.remaining_records.size())
          return false;

        for (size_t i = 0; i < remaining_records.size(); ++i)
          res &= (remaining_records[i] == rhs.remaining_records[i]);

        return res;
    }

    bool operator!=(node_data const & rhs) const
    {
        return !(*this == rhs);
    }
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

std::ostream & operator<<(std::ostream & s, node_data<chopper_split_record> const & nd)
{
    s << "{parent:"
      << nd.parent_bin_index
      << ", max:"
      << nd.max_bin_index
      << ", num:"
      << nd.number_of_technical_bins
      << ", fav:"
      << ((nd.favourite_child == lemon::INVALID) ? "NONE" : "some")
      << ", recs:[";
    for (auto const & rec : nd.remaining_records)
        s << rec << ',';
    s  << "]}";

    return s;
}
