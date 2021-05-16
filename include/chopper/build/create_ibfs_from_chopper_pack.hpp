#pragma once

#include <fstream>
#include <unordered_set>
#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
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

inline std::unordered_set<size_t> compute_kmers(build_config const & config, chopper_pack_record const & record)
{
    std::unordered_set<size_t> kmers{};

    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);

    return kmers;
}

inline void compute_kmers(std::unordered_set<size_t> & kmers,
                          build_config const & config,
                          chopper_pack_record const & record)
{
    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);
}

// automatically does naive splitting if number_of_bins > 1
inline void insert_into_ibf(std::unordered_set<size_t> const & kmers,
                            size_t const number_of_bins,
                            size_t const bin_index,
                            seqan3::interleaved_bloom_filter<> & ibf)
{
    size_t const kmers_per_chunk = (kmers.size() / number_of_bins) + 1;
    auto it = kmers.begin();
    for (size_t chunk = 0; chunk < number_of_bins; ++chunk)
        for (size_t i = chunk * kmers_per_chunk; i < std::min((chunk + 1) * kmers_per_chunk, kmers.size()); ++i, ++it)
            ibf.emplace(*it, seqan3::bin_index{bin_index + chunk});
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
inline void build(std::unordered_set<size_t> & parent_kmers,
                  lemon::ListDigraph::Node const & current_node,
                  build_data<chopper_pack_record> & data,
                  build_config const & config);

inline void update_user_bins(build_data<chopper_pack_record> & data, std::vector<int64_t> & ibf_filenames, chopper_pack_record const & record)
{
    auto const user_bin_pos = data.user_bins.add_user_bin(record.filenames);
    for (size_t i = 0; i < record.number_of_bins.back(); ++i)
        ibf_filenames[record.bin_indices.back() + i] = user_bin_pos;
}

inline size_t initialise_max_bin_kmers(std::unordered_set<size_t> & kmers,
                                       std::vector<int64_t> & ibf_positions,
                                       std::vector<int64_t> & ibf_filenames,
                                       lemon::ListDigraph::Node const & node,
                                       build_data<chopper_pack_record> & data,
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
        auto const & record = node_data.remaining_records[0];
        compute_kmers(kmers, config, record);
        update_user_bins(data, ibf_filenames, record);

        return record.number_of_bins.back();
    }
}

class MyNumPunct : public std::numpunct<char>
{
protected:
    virtual char do_thousands_sep() const { return ','; }
    virtual std::string do_grouping() const { return "\03"; }
};

inline auto construct_ibf(std::unordered_set<size_t> & kmers,
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

    insert_into_ibf(kmers, number_of_bins, node_data.max_bin_index, ibf);

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
    std::unordered_set<size_t> kmers{};

    for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
    {
        kmers.clear();
        auto child = data.ibf_graph.target(arc_it);
        auto & child_data = data.node_map[child];
        if (child != current_node_data.favourite_child)
        {
            build(kmers, child, data, config); // also appends that childs counts to 'kmers'
            insert_into_ibf(kmers, 1, child_data.parent_bin_index, ibf);
            ibf_positions[child_data.parent_bin_index] = data.hibf.size();

            if constexpr (std::same_as<parent_kmers_type, std::unordered_set<size_t>>)
                parent_kmers.merge(kmers);
        }
    }
}

inline void build(std::unordered_set<size_t> & parent_kmers,
                  lemon::ListDigraph::Node const & current_node,
                  build_data<chopper_pack_record> & data,
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
    std::unordered_set<size_t> kmers{};
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];
        compute_kmers(kmers, config, record);
        insert_into_ibf(kmers, record.number_of_bins.back(), record.bin_indices.back(), ibf);
        parent_kmers.merge(kmers);
        update_user_bins(data, ibf_filenames, record);
        kmers.clear();
    }

    data.hibf.push_back(std::move(ibf));

    for (auto & pos : ibf_positions)
        if (pos == -1)
            pos = data.hibf.size();

    data.hibf_bin_levels.push_back(std::move(ibf_positions));
    data.user_bins.add_user_bin_positions(std::move(ibf_filenames));
}

inline void create_ibfs_from_chopper_pack(build_data<chopper_pack_record> & data, build_config const & config)
{
    read_chopper_pack_file(data, config.chopper_pack_filename);

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

        if (record.number_of_bins.back() == 1) // no splitting needed
        {
            insert_into_ibf(config, record, high_level_ibf);
        }
        else
        {
            max_bin_kmers.clear(); // reuse allocated space
            compute_kmers(max_bin_kmers, config, record);
            insert_into_ibf(max_bin_kmers, record.number_of_bins.back(), record.bin_indices.back(), high_level_ibf);
        }

        update_user_bins(data, ibf_filenames, record);
    }

    data.hibf.insert(data.hibf.begin(), std::move(high_level_ibf)); // insert High level at the beginning
    data.hibf_bin_levels.insert(data.hibf_bin_levels.begin(), std::move(ibf_positions));
    data.user_bins.prepend_user_bin_positions(std::move(ibf_filenames));
}
