#pragma once

#include <fstream>
#include <unordered_set>
#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
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

std::unordered_set<size_t> compute_kmers(build_config const & config, chopper_pack_record const & record)
{
    std::unordered_set<size_t> kmers{};

    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);

    return kmers;
}

void compute_kmers(std::unordered_set<size_t> & kmers,
                   build_config const & config,
                   chopper_pack_record const & record)
{
    for (auto const & filename : record.filenames)
        for (auto && [seq] : sequence_file_t{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                kmers.insert(hash);
}

void insert_into_ibf(std::unordered_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_bloom_filter<> & ibf)
{
    size_t const kmers_per_chunk = kmers.size() / number_of_bins + 1;
    auto it = kmers.begin();
    for (size_t chunk = 0; chunk < number_of_bins; ++chunk)
        for (size_t i = chunk * kmers_per_chunk; i < std::min((chunk + 1) * kmers_per_chunk, kmers.size()); ++i, ++it)
            ibf.emplace(*it, seqan3::bin_index{bin_index + chunk});
}

template <typename low_level_ibf_type>
void insert_into_ibf(build_config const & config,
                     chopper_pack_record const & record,
                     seqan3::interleaved_bloom_filter<> & hibf,
                     low_level_ibf_type & libf)
{
    for (auto const & filename : record.filenames)
    {
        for (auto && [seq] : sequence_file_t{filename})
        {
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
            {
                assert(record.hidx >= 0);
                hibf.emplace(hash, seqan3::bin_index{static_cast<size_t>(record.hidx)});

                if constexpr (std::same_as<low_level_ibf_type, seqan3::interleaved_bloom_filter<>>)
                {
                    assert(record.lidx >= 0);
                    libf.emplace(hash, seqan3::bin_index{static_cast<size_t>(record.lidx)});
                }
            }
        }
    }
}

auto initialise_hibf(build_config const & config,
                     build_data const & data,
                     std::vector<std::vector<chopper_pack_record>> const & records_per_hibf_bin,
                     std::vector<seqan3::interleaved_bloom_filter<>> & low_level_ibfs)
{
    assert(data.hibf_num_technical_bins != 0);
    std::unordered_set<size_t> hibf_kmers{};
    size_t hibf_highest_record_bins{1};

    // since the packing is supposedly done with minimizers and we probably have a different setting now
    // we need to calculate the maximum bin size. Since the relative number of minimizers should always correlate

    // we can estimate it by calculating the kmer content of the "highest bin".
    size_t const hidx = data.hibf_max_bin;
    auto const & highest_records = records_per_hibf_bin[hidx];
    assert(highest_records.size() > 0);

    if (config.verbose)
        seqan3::debug_stream << ">>> Initialising HIBF with bin " << hidx << std::endl;

    if (highest_records.size() == 1) // split bin
    {
        hibf_kmers = compute_kmers(config, highest_records[0]);
        hibf_highest_record_bins = highest_records[0].bins;
    }
    else // merged bin
    {
        // in this case we also initialize the respective low level IBF.
        size_t const max_id = data.merged_max_bin_map.at(hidx);
        auto max_it = std::ranges::find_if(highest_records, [max_id] (auto const & rec) {return rec.lidx == max_id;});

        std::unordered_set<size_t> libf_kmers{compute_kmers(config, *max_it)};

        // construct LIBF
        seqan3::bin_size const l_bin_size{compute_bin_size(config, libf_kmers.size() / max_it->bins)};
        seqan3::bin_count const l_bin_count{highest_records.back().lidx + highest_records.back().bins};
        seqan3::interleaved_bloom_filter libf{l_bin_count, l_bin_size, seqan3::hash_function_count{config.hash_funs}};

        if (config.verbose)
            seqan3::debug_stream << "  > Initialised LIBF with bin size " << l_bin_size.get() << std::endl;

        insert_into_ibf(libf_kmers, max_it->bins, max_it->lidx, libf);

        hibf_kmers.merge(libf_kmers);

        for (auto it = highest_records.begin(); it != highest_records.end(); ++it)
        {
            if (it != max_it) // already processed
            {
                libf_kmers = compute_kmers(config, *it);
                insert_into_ibf(libf_kmers, it->bins, it->lidx, libf);
                hibf_kmers.merge(libf_kmers);
            }
        }

        low_level_ibfs[hidx] = std::move(libf);
        hibf_highest_record_bins = 1;
    }

    seqan3::bin_size const hibf_bin_size{compute_bin_size(config, hibf_kmers.size() / hibf_highest_record_bins)};

    // construct HIBF
    seqan3::interleaved_bloom_filter high_level_ibf{seqan3::bin_count{data.hibf_num_technical_bins},
                                                    hibf_bin_size,
                                                    seqan3::hash_function_count{config.hash_funs}};

    if (config.verbose)
        seqan3::debug_stream << "  > Initialised HIBF with bin size " << hibf_bin_size.get() << std::endl;

    insert_into_ibf(hibf_kmers, hibf_highest_record_bins, hidx, high_level_ibf);

    return high_level_ibf;
}

auto process_merged_bin(build_config const & config,
                        std::vector<chopper_pack_record> const & records,
                        build_data const & data,
                        seqan3::interleaved_bloom_filter<> & high_level_ibf,
                        std::vector<seqan3::interleaved_bloom_filter<>> & low_level_ibfs)
{
    // we need to construct a low level ibf AND insert all kmers into the high level bin
    assert(records.size() > 1);
    size_t const max_id = data.merged_max_bin_map.at(records[0].hidx);
    auto max_it = std::ranges::find_if(records, [max_id] (auto const & rec) {return rec.lidx == max_id;});

    std::unordered_set<size_t> kmers{compute_kmers(config, *max_it)};

    // construct LIBF
    seqan3::bin_size const l_bin_size{compute_bin_size(config, kmers.size() / max_it->bins)};
    seqan3::bin_count const l_bin_count{records.back().lidx + records.back().bins};
    seqan3::interleaved_bloom_filter libf{l_bin_count, l_bin_size, seqan3::hash_function_count{config.hash_funs}};

    if (config.verbose)
        seqan3::debug_stream << "  > Initialised LIBF with bin size " << l_bin_size.get() << std::endl;

    insert_into_ibf(kmers, max_it->bins, max_it->lidx, libf);
    insert_into_ibf(kmers, 1, max_it->hidx, high_level_ibf);

    for (auto it = records.begin(); it != records.end(); ++it)
    {
        if (it != max_it) // already processed
        {
            if (it->bins == 1) // no splitting
            {
                insert_into_ibf(config, *it, high_level_ibf, libf);
            }
            else
            {
                kmers = compute_kmers(config, *it);
                insert_into_ibf(kmers, it->bins, it->lidx, libf);
                insert_into_ibf(kmers, 1, it->hidx, high_level_ibf);
            }
        }
    }

    low_level_ibfs[records[0].hidx] = std::move(libf);
}

auto process_split_bin(build_config const & config,
                       std::vector<chopper_pack_record> const & records,
                       seqan3::interleaved_bloom_filter<> & high_level_ibf)
{
    assert(records.size() == 1);
    auto const & record = records[0];
    if (record.bins == 1) // no splitting
    {
        insert_into_ibf(config, record, high_level_ibf, /*something that is not an ibf:*/record);
    }
    else
    {
        std::unordered_set<size_t> kmers{compute_kmers(config, record)};
        insert_into_ibf(kmers, record.bins, record.hidx, high_level_ibf);
    }
}

// node muss speichern:
// - IBF (referenz oder position)
// - max bin id
// - user bins
// - all the other user bins in this IBF

void build(std::unordered_set<size_t> & parent_kmers,
           const & current_node,
           const & tree,
           build_config const & config)
{
    std::unordered_set<size_t> current_node_kmers{};

    if (there is a favourite child, I am not a leaf) // favourite child -> max bin is a merged bin
    {
        initialize_favourite_child();
    }
    else // there a max bin, that is ot a merged bin
    {
        compute_kmers(current_node_kmers, config, max_record/*one line in file*/);
    }

    initialize IBF with respective bin size

    insert current_node_kmers into IBF

    while (there is another child that needs to be initialized beforehand) // (can be more than one child)
    {
        std::unordered_set<size_t> kmers{};
        initialize(kmers, child, tree); // also appends that childs counts to 'current_node_kmers'
        insert kmers into bin in IBF
        parent_kmers.merge(kmers);
    }

    hash and insert all remaining user bins of this IBF and into current_node_kmers;

    parent_kmers.merge(current_node_kmers);
}

void start_build(const & tree, build_config const & config)
{
    current_node = root; // high level

    std::unordered_set<size_t> current_node_kmers{};

    if (there is a favourite child, I am not a leaf) // favourite child -> max bin is a merged bin
    {
        initialize_favourite_child();
    }
    else // there a max bin, that is ot a merged bin
    {
        compute_kmers(current_node_kmers, config, max_record/*one line in file*/);
    }

    initialize IBF with respective bin size

    insert current_node_kmers into IBF

    while (there is another child that needs to be initialized beforehand) // (can be more than one child)
    {
        std::unordered_set<size_t> kmers{};
        initialize(kmers, child, tree); // also appends that childs counts to 'current_node_kmers'
        insert kmers into bin in IBF
    }

    hash and insert all remaining user bins of this IBF;
}

auto create_ibfs_from_chopper_pack(build_config const & config)
{
    auto const [data, records_per_hibf_bin] = read_chopper_pack_file(config.chopper_pack_filename);

    // fill libfs with dummy ibfs to resize the vector already
    seqan3::interleaved_bloom_filter dummy{seqan3::bin_count{1}, seqan3::bin_size{1}, seqan3::hash_function_count{1}};
    std::vector<seqan3::interleaved_bloom_filter<>> low_level_ibfs(data.hibf_num_technical_bins, dummy);

    auto high_level_ibf = initialise_hibf(config, data, records_per_hibf_bin, low_level_ibfs);

    for (auto const & records : records_per_hibf_bin)
    {
        if (records.size() == 0 || records[0].hidx == data.hibf_max_bin) // do not process empty bins or the highest
        {
            continue;
        }
        else if (records.size() == 1) // split bin
        {
            if (config.verbose)
            {
                seqan3::debug_stream << ">>> Processing split bin(s) "
                                     << records[0].hidx
                                     << ((records[0].bins > 1) ? "-" : "")
                                     << ((records[0].bins > 1) ? std::to_string(records[0].hidx + records[0].bins) : "")
                                     << " (out of " << data.hibf_num_technical_bins << ")" << std::endl;
            }

            process_split_bin(config, records, high_level_ibf);
        }
        else
        {
            if (config.verbose)
            {
                seqan3::debug_stream << ">>> Processing merged bin " << records[0].hidx
                                     << " (out of " << data.hibf_num_technical_bins << ")" << std::endl;
            }

            process_merged_bin(config, records, data, high_level_ibf, low_level_ibfs);
        }
    }

    return std::make_tuple(std::move(high_level_ibf), std::move(low_level_ibfs));
}
