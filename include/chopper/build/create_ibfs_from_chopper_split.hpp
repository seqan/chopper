#pragma once

#include <unordered_map>
#include <fstream>

#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <chopper/build/batch.hpp>
#include <chopper/build/build_config.hpp>
#include <chopper/build/compute_bin_size.hpp>
#include <chopper/build/read_data_file_and_set_high_level_bins.hpp>
#include <chopper/build/read_chopper_split_file.hpp>
#include <chopper/detail_bin_prefixes.hpp>

struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using seq_file_type = seqan3::sequence_file_input<file_type_traits,
                                                  seqan3::fields<seqan3::field::seq, seqan3::field::id>,
                                                  seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

using seq_file_type2 = seqan3::sequence_file_input<file_type_traits,
                                                   seqan3::fields<seqan3::field::seq>,
                                                   seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

auto read_sequences(std::vector<std::string> const & filenames)
{
    std::unordered_map<std::string, seqan3::dna4_vector> info;

    for (auto const & filename : filenames)
        for (auto && [seq, id] : seq_file_type{filename})
            info.emplace(filename + id, std::move(seq));

    return info;
}

auto hash_infix(build_config const & config, auto const & seq, auto const begin, auto const end)
{
    return seq | seqan3::views::drop(begin)
               | seqan3::views::take(end + config.overlap - begin) // views::take never goes over the end
               | seqan3::views::kmer_hash(seqan3::ungapped{config.k});
};

// TODO this will miss the low level IBF of a merged bin currently!
// But it's no bug yet as the highest record is parsed again in the create_ibfs loop below.
auto initialise_hibf(build_config const & config, build_data const & data)
{
    assert(data.hibf_num_technical_bins != 0);
    std::unordered_set<size_t> kmers{};

    // since the packing is supposedly done with minimizers and we probably have a different setting now
    // we need to calculate the maximum bin size. Since the relative number of minimizers should always correlate
    // wen can estimate it by calculating the kmer content of the "highest bin".
    auto const & highest_batch = *data.hibf_max_batch_record;

    std::vector<std::string> ids{};
    std::vector<seqan3::dna4_vector> seqs{};
    if (highest_batch.hibf_bins.size() != 1) // split bin
    {
        for (auto const & filename : highest_batch.filenames)
        {
            for (auto && [seq, id] : seq_file_type{filename})
            {
                for (auto const & reg : data.region_map.at(filename + id))
                {
                    if (reg.hidx == data.hibf_max_bin)
                        for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
                            kmers.insert(hash);

                    ids.push_back(filename + id);
                    seqs.push_back(seq);
                }
            }
        }
    }
    else
    {
        for (auto const & filename : highest_batch.filenames)
            for (auto && [seq] : seq_file_type2{filename})
                for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                    kmers.insert(hash);
    }

    seqan3::bin_size hibf_bin_size{compute_bin_size(config, kmers.size())};

    // construct HIBF
    seqan3::interleaved_bloom_filter high_level_ibf{seqan3::bin_count{data.hibf_num_technical_bins},
                                                    hibf_bin_size,
                                                    seqan3::hash_function_count{config.hash_funs}};

    // insert already hashed kmers
    for (auto const hash : kmers)
        high_level_ibf.emplace(hash, seqan3::bin_index{data.hibf_max_bin});

    if (highest_batch.hibf_bins.size() != 1) // it was a split bin and we need to hash the rest
    {
        for (size_t i = 0; i < ids.size(); ++i)
            for (auto const & reg : data.region_map.at(ids[i]))
                if (reg.hidx != data.hibf_max_bin)
                    for (auto hash : hash_infix(config, seqs[i], reg.begin, reg.end))
                        kmers.insert(hash);
    }

    if (config.verbose)
        seqan3::debug_stream << ">>> Initialised HIBF with bin size " << hibf_bin_size.get() << std::endl;

    return high_level_ibf;
}

auto process_bin(build_config const & config,
                 batch const & batch_record,
                 build_data const & data,
                 seqan3::interleaved_bloom_filter<> & high_level_ibf,
                 std::vector<seqan3::interleaved_bloom_filter<>> & low_level_ibfs)
{
    // we need to construct a low level ibf AND insert all kmers into the high level bin
    auto && info = read_sequences(batch_record.filenames);

    if (batch_record.libf_num_bins != 0) // this is a merged bin -> Initialise its IBF
    {
        assert(batch_record.hibf_bins.size() == 1); // a merged record only corresponds to a single bin in the HIBF
        size_t const hidx = batch_record.hibf_bins[0];

        std::unordered_set<size_t> kmers{};
        size_t const highest_libf_bin_idx = data.merged_max_bin_map.at(hidx);

        for (auto const & [combined_id, seq] : info)
            for (auto const & reg : data.region_map.at(combined_id))
                if (reg.lidx == highest_libf_bin_idx)
                    for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
                        kmers.insert(hash);

        seqan3::bin_size libf_bin_size{compute_bin_size(config, kmers.size())};

        // construct LIBF
        seqan3::interleaved_bloom_filter libf{seqan3::bin_count{batch_record.libf_num_bins},
                                              libf_bin_size,
                                              seqan3::hash_function_count{config.hash_funs}};

        if (config.verbose)
            seqan3::debug_stream << "  > Initialised LIBF with bin size " << libf_bin_size.get() << std::endl;

        low_level_ibfs[hidx] = std::move(libf);
    }

    for (auto const & [combined_id, seq] : info)
    {
        for (auto const & reg : data.region_map.at(combined_id))
        {
            // TODO: this will probably be more efficient if regions are sorted by hidx
            for (auto hash : hash_infix(config, seq, reg.begin, reg.end))
            {
                high_level_ibf.emplace(hash, seqan3::bin_index{reg.hidx});

                if (reg.lidx != -1)
                    low_level_ibfs[reg.hidx].emplace(hash, seqan3::bin_index{static_cast<size_t>(reg.lidx)});
            }
        }
    }
}

auto create_ibfs_from_chopper_split(build_config const & config)
{
    auto const [data, batches] = read_chopper_split_file(config.chopper_split_filename);

    auto high_level_ibf = initialise_hibf(config, data);

    std::unordered_map<size_t, size_t> libfs_pos_map{};
    // fill libfs with dummy ibfs to resize the vector already
    seqan3::interleaved_bloom_filter dummy{seqan3::bin_count{1}, seqan3::bin_size{1}, seqan3::hash_function_count{1}};
    std::vector<seqan3::interleaved_bloom_filter<>> low_level_ibfs(data.hibf_num_technical_bins, dummy);

    size_t bin_idx{};
    // todo: leave out highest record but fix initialise_ibf before!
    for (auto const & batch_record : batches)
    {
        if (config.verbose)
        {
            seqan3::debug_stream << ">>> Processing "
                                 << ((batch_record.libf_num_bins > 0) ? "merged bin " : " bin(s) ")
                                 << batch_record.hibf_bins
                                 << " (out of " << data.hibf_num_technical_bins << ")" << std::endl;
        }

        process_bin(config, batch_record, data, high_level_ibf, low_level_ibfs);

        ++bin_idx;
    }

    return std::make_tuple(std::move(high_level_ibf), std::move(low_level_ibfs));
}
