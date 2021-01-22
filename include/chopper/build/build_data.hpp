#pragma once

#include <lemon/list_graph.h>

#include <chopper/build/batch.hpp>
#include <chopper/build/region.hpp>
#include <chopper/detail_node_data.hpp>
#include <chopper/detail_parse_chopper_pack_line.hpp>

struct build_data
{
    size_t hibf_num_technical_bins{};
    size_t num_libfs{};
    std::string hibf_max_bin_id{};
    size_t hibf_max_bin{};
    chopper_pack_record * hibf_max_record{nullptr};
    batch * hibf_max_batch_record{nullptr};

    std::unordered_map<size_t, size_t> merged_max_bin_map{};
    std::unordered_map<std::string, size_t> merged_bin_map{};
    std::unordered_map<std::string, std::vector<region>> region_map{};

    // new stuff
    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data> node_map{ibf_graph};
};
