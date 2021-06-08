#pragma once

#include <lemon/list_graph.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <chopper/build/batch.hpp>
#include <chopper/build/region.hpp>
#include <chopper/detail_node_data.hpp>
#include <chopper/detail_parse_chopper_pack_line.hpp>
#include <chopper/detail_hibf_user_bins.hpp>

template <typename record_type>
struct build_data
{
    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data<record_type>> node_map{ibf_graph};

    std::vector<seqan3::interleaved_bloom_filter<>> hibf;

    // maps for each ibf in hibf, each bin to the next ibf postition in hibf (if it is a merged bin)
    // or to the same ibf (if it is not a merged bin). You can thereby check if you need to query another
    // lower level IBF.
    std::vector<std::vector<int64_t>> hibf_bin_levels;

    hibf_user_bins user_bins;
};
