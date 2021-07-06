#pragma once

#include <atomic>

#include <lemon/list_graph.h>

#include <chopper/build/batch.hpp>
#include <chopper/build/region.hpp>
#include <chopper/detail_node_data.hpp>
#include <chopper/detail_parse_chopper_pack_line.hpp>
#include <chopper/hierarchical_interleaved_bloom_filter.hpp>

template <typename record_type>
struct build_data
{
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> ibf_number{};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> user_bin_number{};

    size_t number_of_user_bins{};
    size_t number_of_ibfs{};

    size_t request_ibf_idx()
    {
        return std::atomic_fetch_add(&ibf_number, 1u);
    }

    size_t request_user_bin_idx()
    {
        return std::atomic_fetch_add(&user_bin_number, 1u);
    }

    void resize()
    {
        hibf.hibf.resize(number_of_ibfs);
        hibf.user_bins.resize_bins(number_of_ibfs);
        hibf.user_bins.resize_filename(number_of_user_bins);
        hibf.hibf_bin_levels.resize(number_of_ibfs);
    }

    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data<record_type>> node_map{ibf_graph};

    hierarchical_interleaved_bloom_filter<> hibf{};
};
