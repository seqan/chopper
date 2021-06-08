#pragma once

#include <fstream>
#include <future>
#include <random>
#include <seqan3/std/ranges>

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <chopper/build/build_config.hpp>
#include <chopper/build/compute_bin_size.hpp>
#include <chopper/build/read_chopper_pack_file.hpp>

struct file_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_t = seqan3::sequence_file_input<file_traits,
                                                    seqan3::fields<seqan3::field::seq>,
                                                    seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

inline robin_hood::unordered_flat_set<size_t> compute_kmers(build_config const & config, chopper_pack_record const & record)
{
    robin_hood::unordered_flat_set<size_t> kmers{};

    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);

    return kmers;
}

inline void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                          build_config const & config,
                          chopper_pack_record const & record)
{
    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);
}

// automatically does naive splitting if number_of_bins > 1
template <typename parent_kmers_type>
inline void insert_into_ibf(parent_kmers_type & parent_kmers,
                            robin_hood::unordered_flat_set<size_t> const & kmers,
                            size_t const number_of_bins,
                            size_t const bin_index,
                            seqan3::interleaved_bloom_filter<> & ibf)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            ibf.emplace(value, bin_idx);
            if constexpr(std::same_as<parent_kmers_type, robin_hood::unordered_flat_set<size_t>>)
                parent_kmers.insert(value);
        }
    }
}

inline void insert_into_ibf(build_config const & config,
                            chopper_pack_record const & record,
                            seqan3::interleaved_bloom_filter<> & hibf)
{
    assert(record.bin_indices.back() >= 0);
    auto const bin_index = seqan3::bin_index{static_cast<size_t>(record.bin_indices.back())};
    auto hash_view = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | hash_view)
                hibf.emplace(hash, bin_index);
}

// forward declaration
template <typename parent_kmers_type>
inline size_t build(parent_kmers_type & parent_kmers,
                    lemon::ListDigraph::Node const & current_node,
                    build_data<chopper_pack_record> & data,
                    build_config const & config);

inline void update_user_bins(build_data<chopper_pack_record> & data, std::vector<int64_t> & filename_indices, chopper_pack_record const & record)
{
    size_t const idx = data.request_user_bin_idx();
    data.user_bins.filename_at(idx) = record.filenames | seqan3::views::join_with(std::string{";"}) | seqan3::views::to<std::string>;
    std::fill_n(filename_indices.begin() + record.bin_indices.back(), record.number_of_bins.back(), idx);
}

inline size_t initialise_max_bin_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                                       std::vector<int64_t> & ibf_positions,
                                       std::vector<int64_t> & filename_indices,
                                       lemon::ListDigraph::Node const & node,
                                       build_data<chopper_pack_record> & data,
                                       build_config const & config)
{
    auto & node_data = data.node_map[node];

    if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
    {
        // recursively initialize favourite child first
        ibf_positions[node_data.max_bin_index] = build(kmers, node_data.favourite_child, data, config);
        return 1;
    }
    else // max bin is not a merged bin
    {
        // we assume that the max record is at the beginning of the list of remaining records.
        auto const & record = node_data.remaining_records[0];
        compute_kmers(kmers, config, record);
        update_user_bins(data, filename_indices, record);

        return record.number_of_bins.back();
    }
}

class MyNumPunct : public std::numpunct<char>
{
protected:
    virtual char do_thousands_sep() const { return ','; }
    virtual std::string do_grouping() const { return "\03"; }
};

template <typename parent_kmers_type>
inline auto construct_ibf(parent_kmers_type & parent_kmers,
                          robin_hood::unordered_flat_set<size_t> & kmers,
                          size_t const number_of_bins,
                          lemon::ListDigraph::Node const & node,
                          build_data<chopper_pack_record> & data,
                          build_config const & config)
{
    auto & node_data = data.node_map[node];

    seqan3::bin_size const bin_size{compute_bin_size(config, kmers.size() / number_of_bins)};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter ibf{bin_count, bin_size, seqan3::hash_function_count{config.hash_funs}};

    if (config.verbose)
    {
        std::cout.imbue( std::locale( std::locale::classic(), new MyNumPunct ) );
        std::cout << "  > Initialised IBF with bin size " << bin_size.get() << std::endl;
    }

    insert_into_ibf(parent_kmers, kmers, number_of_bins, node_data.max_bin_index, ibf);

    return ibf;
}

template <typename parent_kmers_type>
void loop_over_children(parent_kmers_type & parent_kmers,
                        seqan3::interleaved_bloom_filter<> & ibf,
                        std::vector<int64_t> & ibf_positions,
                        lemon::ListDigraph::Node const & current_node,
                        build_data<chopper_pack_record> & data,
                        build_config const & config)
{
    auto & current_node_data = data.node_map[current_node];
    std::vector<lemon::ListDigraph::Node> children{};

    for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
        children.emplace_back(data.ibf_graph.target(arc_it));

    if (children.empty())
        return;

    size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
    std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

    auto worker = [&] (auto && index, auto &&)
    {
        auto & child = children[index];

        if (child != current_node_data.favourite_child)
        {
            robin_hood::unordered_flat_set<size_t> kmers{};
            size_t const ibf_pos = build(kmers, child, data, config);
            auto parent_bin_index = data.node_map[child].parent_bin_index;
            {
                size_t const mutex_id{parent_bin_index / 64};
                std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                ibf_positions[parent_bin_index] = ibf_pos;
                insert_into_ibf(parent_kmers, kmers, 1, parent_bin_index, ibf);
            }
        }
    };

    size_t number_of_threads{};
    auto indices_view = std::views::iota(0u, children.size()) | std::views::common;
    std::vector<size_t> indices{indices_view.begin(), indices_view.end()};

    if constexpr(!std::same_as<parent_kmers_type, robin_hood::unordered_flat_set<size_t>>)
    {
        // Shuffle indices: More likely to not block each other. Optimal: Interleave
        std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
        number_of_threads = config.threads;
    }
    else
    {
        number_of_threads = 1u;
    }

    seqan3::detail::execution_handler_parallel executioner{number_of_threads};
    executioner.bulk_execute(std::move(worker), std::move(indices), [] () {});
}

template <typename parent_kmers_type>
inline size_t build(parent_kmers_type & parent_kmers,
                    lemon::ListDigraph::Node const & current_node,
                    build_data<chopper_pack_record> & data,
                    build_config const & config)
{
    constexpr bool is_root = !std::same_as<parent_kmers_type, robin_hood::unordered_flat_set<size_t>>;
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<size_t> kmers{};

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers(kmers, ibf_positions, filename_indices, current_node, data, config);
    auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, config);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, config);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        if constexpr (is_root)
        {
            if (record.number_of_bins.back() == 1) // no splitting needed
            {
                insert_into_ibf(config, record, ibf);
            }
            else // same as non-root
            {
                compute_kmers(kmers, config, record);
                insert_into_ibf(parent_kmers, kmers, record.number_of_bins.back(), record.bin_indices.back(), ibf);
            }
        }
        else
        {
            compute_kmers(kmers, config, record);
            insert_into_ibf(parent_kmers, kmers, record.number_of_bins.back(), record.bin_indices.back(), ibf);
        }

        update_user_bins(data, filename_indices, record);
        kmers.clear();
    }

    data.hibf[ibf_pos] = std::move(ibf);
    data.hibf_bin_levels[ibf_pos] = std::move(ibf_positions);
    data.user_bins.bin_at(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}

inline void create_ibfs_from_chopper_pack(build_data<chopper_pack_record> & data, build_config const & config)
{
    read_chopper_pack_file(data, config.chopper_pack_filename);
    lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
    size_t dummy{};

    build(dummy, root, data, config);
}
