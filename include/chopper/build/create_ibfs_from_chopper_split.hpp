#pragma once

#if 0 // split is not implemented for multi-level

#include <fstream>
#include <seqan3/std/ranges>
#include <unordered_map>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <chopper/build/batch.hpp>
#include <chopper/build/build_config.hpp>
#include <chopper/build/compute_bin_size.hpp>
#include <chopper/build/read_chopper_split_file.hpp>
#include <chopper/detail_bin_prefixes.hpp>

struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

inline auto hash_infix(build_config const & config, auto const & seq, auto const begin, auto const end)
{
    return seq | seqan3::views::slice(begin, end + config.overlap)
               | seqan3::views::kmer_hash(seqan3::ungapped{config.k});
};

inline void compute_kmers(std::unordered_set<size_t> & kmers,
                          build_config const & config,
                          chopper_split_record const & record)
{
    for (auto const & [combined_id, seq] : record.info)
        for (auto const & reg : record.region_map.at(combined_id))
            if (reg.bin_index == record.bin_indices.front())
                for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
                    kmers.insert(hash);
}

inline void insert_into_ibf_and_set(build_config const & config,
                                    chopper_split_record const & record,
                                    std::unordered_set<size_t> & kmers,
                                    seqan3::interleaved_bloom_filter<> & ibf)
{
    for (auto const & [combined_id, seq] : record.info)
    {
        for (auto const & reg : record.region_map.at(combined_id))
        {
            // TODO: this will probably be more efficient if regions are sorted by hidx
            for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
            {
                ibf.emplace(hash, seqan3::bin_index{reg.bin_index});
                kmers.insert(hash);
            }
        }
    }
}

inline void insert_into_ibf_from_set(std::unordered_set<size_t> const & kmers,
                                     size_t const bin_index,
                                     seqan3::interleaved_bloom_filter<> & ibf)
{
    for (auto const kmer : kmers)
        ibf.emplace(kmer, seqan3::bin_index{bin_index});
}

inline void insert_into_ibf(build_config const & config,
                            chopper_split_record const & record,
                            seqan3::interleaved_bloom_filter<> & ibf)
{
    using sequence_file_t = seqan3::sequence_file_input<file_type_traits,
                                                        seqan3::fields<seqan3::field::seq>,
                                                        seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

    for (auto const & filename : record.filenames)
    {
        for (auto && [seq] : sequence_file_t{filename})
        {
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
            {
                assert(record.bin_indices.back() >= 0);
                ibf.emplace(hash, seqan3::bin_index{static_cast<size_t>(record.bin_indices.back())});
            }
        }
    }
}

inline void insert_regions_into_ibf(build_config const & config,
                                    chopper_split_record const & record,
                                    seqan3::interleaved_bloom_filter<> & ibf)
{
    for (auto const & [combined_id, seq] : record.info)
    {
        for (auto const & reg : record.region_map.at(combined_id))
        {
            // TODO: this will probably be more efficient if regions are sorted by hidx
            for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
            {
                ibf.emplace(hash, seqan3::bin_index{reg.bin_index});
            }
        }
    }
}

inline void update_user_bins(build_data<chopper_split_record> & data,
                             std::vector<int64_t> & ibf_filenames,
                             chopper_split_record const & record)
{
    auto const user_bin_pos = data.user_bins.add_user_bin(record.filenames);
    for (size_t i = 0; i < record.number_of_bins.back(); ++i)
        ibf_filenames.at(record.bin_indices.back() + i) = user_bin_pos;
}

inline void build(std::unordered_set<size_t> & parent_kmers,
                  lemon::ListDigraph::Node const & current_node,
                  build_data<chopper_split_record> & data,
                  build_config const & config);

inline size_t initialise_max_bin_kmers(std::unordered_set<size_t> & kmers,
                                       std::vector<int64_t> & ibf_positions,
                                       std::vector<int64_t> & ibf_filenames,
                                       lemon::ListDigraph::Node const & node,
                                       build_data<chopper_split_record> & data,
                                       build_config const & config)
{
    auto & node_data = data.node_map[node];

    if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
    {
        build(kmers, node_data.favourite_child, data, config); // recursively initialize favourite child first
        ibf_positions[node_data.max_bin_index] = data.hibf.size();
        return 1;
    }
    else // there a max bin, that is ot a merged bin
    {
        // we assume that the max record is at the beginning of the list of remaining records.
        auto & record = node_data.remaining_records[0];
        record.initialise_info();
        assert(node_data.max_bin_index == record.bin_indices.front());
        compute_kmers(kmers, config, record);
        update_user_bins(data, ibf_filenames, record);

        return record.bin_indices.size();
    }
}

inline auto construct_ibf(std::unordered_set<size_t> & kmers,
                          size_t const number_of_bins,
                          lemon::ListDigraph::Node const & node,
                          build_data<chopper_split_record> & data,
                          build_config const & config)
{
    auto & node_data = data.node_map[node];

    seqan3::bin_size const bin_size{compute_bin_size(config, kmers.size() / number_of_bins)};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter ibf{bin_count, bin_size, seqan3::hash_function_count{config.hash_funs}};

    if (config.verbose)
        std::cout << "  > Initialised IBF with bin size " << bin_size.get() << std::endl;

    insert_into_ibf_from_set(kmers, node_data.max_bin_index, ibf);

    return ibf;
}

template <typename parent_kmers_type>
void loop_over_children(parent_kmers_type & parent_kmers,
                        seqan3::interleaved_bloom_filter<> & ibf,
                        std::vector<int64_t> & ibf_positions,
                        lemon::ListDigraph::Node const & current_node,
                        build_data<chopper_split_record> & data,
                        build_config const & config)
{
    auto & current_node_data = data.node_map[current_node];

    for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
    {
        auto child = data.ibf_graph.target(arc_it);
        auto & child_data = data.node_map[child];
        if (child != current_node_data.favourite_child)
        {
            std::unordered_set<size_t> kmers{}; // todo: maybe it is more efficient if this is declared outside and cleared every iteration
            build(kmers, child, data, config); // also appends that childs counts to 'kmers'
            insert_into_ibf_from_set(kmers, child_data.parent_bin_index, ibf);
            ibf_positions[child_data.parent_bin_index] = data.hibf.size();

            if constexpr (std::same_as<parent_kmers_type, std::unordered_set<size_t>>)
                parent_kmers.merge(kmers);
        }
    }
}

inline void build(std::unordered_set<size_t> & parent_kmers,
                  lemon::ListDigraph::Node const & current_node,
                  build_data<chopper_split_record> & data,
                  build_config const & config)
{
    auto & current_node_data = data.node_map[current_node];

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, -1);
    std::vector<int64_t> ibf_filenames(current_node_data.number_of_technical_bins, -1);
    std::unordered_set<size_t> max_bin_kmers{};

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers(max_bin_kmers, ibf_positions, ibf_filenames, current_node, data, config);
    auto && ibf = construct_ibf(max_bin_kmers, max_bin_tbs, current_node, data, config);
    parent_kmers.merge(max_bin_kmers);
    max_bin_kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, config);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];
        std::unordered_set<size_t> kmers{};
        insert_into_ibf_and_set(config, record, kmers, ibf);
        parent_kmers.merge(kmers);
        update_user_bins(data, ibf_filenames, record);
    }

    data.hibf.push_back(std::move(ibf));

    for (auto & pos : ibf_positions)
        if (pos == -1)
            pos = data.hibf.size();

    data.hibf_bin_levels.push_back(std::move(ibf_positions));
    data.user_bins.add_user_bin_positions(std::move(ibf_filenames));
}

inline void create_ibfs_from_chopper_split(build_data<chopper_split_record> & data, build_config const & config)
{
    read_chopper_split_file(data, config.chopper_split_filename);

    lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
    auto & root_node_data = data.node_map[root];

    std::vector<int64_t> ibf_positions(root_node_data.number_of_technical_bins, 0);
    std::vector<int64_t> ibf_filenames(root_node_data.number_of_technical_bins, -1);
    std::unordered_set<size_t> max_bin_kmers{};

    // initialize high level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers(max_bin_kmers, ibf_positions, ibf_filenames, root, data, config);
    auto && high_level_ibf = construct_ibf(max_bin_kmers, max_bin_tbs, root, data, config);
    max_bin_kmers.clear(); // clear memory allocation

    // parse all other children (merged bins) of the high level-ibf (will build the whole HIBF)
    size_t dummy;
    loop_over_children(dummy, high_level_ibf, ibf_positions, root, data, config);

    // parse remaining (split) bins
    size_t const start{(root_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < root_node_data.remaining_records.size(); ++i)
    {
        auto const & record = root_node_data.remaining_records[i];

        if (record.number_of_bins.size() == 1)
        {
            insert_into_ibf(config, record, high_level_ibf);
        }
        else
        {
            max_bin_kmers.clear(); // reuse allocated space
            insert_regions_into_ibf(config, record, high_level_ibf);
        }

        update_user_bins(data, ibf_filenames, record);
    }

    data.hibf.insert(data.hibf.begin(), std::move(high_level_ibf)); // insert High level at the beginning
    data.hibf_bin_levels.insert(data.hibf_bin_levels.begin(), std::move(ibf_positions));
    data.user_bins.prepend_user_bin_positions(std::move(ibf_filenames));
}

#endif
