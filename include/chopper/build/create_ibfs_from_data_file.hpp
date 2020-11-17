#pragma once

#include <unordered_map>
#include <fstream>

#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

constexpr std::string_view merged_bin_prefix{"COLORFUL_MERGED_BIN"};
constexpr std::string_view split_bin_prefix{"SPLIT_BIN"};

#include <chopper/build/build_config.hpp>
#include <chopper/build/parse_traversal_file_line.hpp>
#include <chopper/build/read_data_file_and_set_high_level_bins.hpp>

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq, seqan3::field::id>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

auto read_sequences(std::vector<std::string> const & filenames)
{
    std::unordered_map<std::string, std::unordered_map<std::string, seqan3::dna4_vector>> info;

    for (auto const & filename : filenames)
        for (auto && [seq, id] : sequence_file_type{filename})
            info[filename][id] = seq;

    return info;
}

auto hash_infix(build_config const & config, auto const & seq, auto const begin, auto const end)
{
    return seq | seqan3::views::drop(begin)
               | seqan3::views::take(end + config.overlap - begin) // views::take never goes over the end
               | seqan3::views::kmer_hash(seqan3::ungapped{config.k});
};

size_t compute_bin_size(build_config const & config, size_t const number_of_kmers_to_be_stored)
{
    return (size_t)std::ceil((-(double)(config.hash_funs * number_of_kmers_to_be_stored)) /
                   (std::log(1 - std::pow(10.0, std::log10(config.FPR) / config.hash_funs))));
}

auto compute_maximum_technical_bin_size(build_config const & config, data_file_record const & record)
{
    auto && info = read_sequences(record.filenames); // into random access containers

    std::set<size_t> kmer_count{};

    if (starts_with(record.bin_name, split_bin_prefix) && record.bins != 1)
    {
        // choose bin index 0 as a representative
        std::ifstream fin{config.traversal_path_prefix + record.bin_name + ".out"};

        if (!fin.good() || !fin.is_open())
            throw std::logic_error{"Could not open file '" + config.traversal_path_prefix +
                                   record.bin_name + ".out' for reading."};

        std::string line;
        std::getline(fin, line); // skip header

        while (std::getline(fin, line))
        {
            auto && [filename, id, begin, end, idx] = parse_traversal_file_line(line);

            if (idx == 0)
                for (auto hash : hash_infix(config, info[filename][id], begin, end))
                    kmer_count.insert(hash);
        }
    }
    else // merged bin or record.bin = 1
    {
        for (auto const & filename : record.filenames)
            for (auto && [seq, id] : sequence_file_type{filename})
                for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                    kmer_count.insert(hash);
    }

    return compute_bin_size(config, kmer_count.size());
}

auto process_splitted_bin(build_config const & config,
                          data_file_record const & record,
                          seqan3::interleaved_bloom_filter<> & high_level_ibf,
                          size_t & bin_idx)
{
    auto && info = read_sequences(record.filenames); // into random access containers

    if (record.bins != 1)
    {
        std::ifstream fin{config.traversal_path_prefix + record.bin_name + ".out"};

        if (!fin.good() || !fin.is_open())
            throw std::logic_error{"Could not open file '" + config.traversal_path_prefix + record.bin_name + ".out' for reading."};

        std::string line;
        std::getline(fin, line); // skip header

        while (std::getline(fin, line))
        {
            auto && [filename, id, begin, end, idx] = parse_traversal_file_line(line);
            assert(bin_idx + idx < high_level_ibf.bin_count());

            // insert into high level
            for (auto hash : hash_infix(config, info[filename][id], begin, end))
                high_level_ibf.emplace(hash, seqan3::bin_index{bin_idx + idx});
        }
        bin_idx += record.bins - 1;
    }
    else
    {
        assert(bin_idx < high_level_ibf.bin_count());
        for (auto const & filename : record.filenames)
            for (auto && [seq, id] : sequence_file_type{filename})
                for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                    high_level_ibf.emplace(hash, seqan3::bin_index{bin_idx});
    }
}

auto process_merged_bin(build_config const & config,
                        data_file_record const & record,
                        seqan3::interleaved_bloom_filter<> & high_level_ibf,
                        std::vector<seqan3::interleaved_bloom_filter<>> & low_level_ibfs,
                        size_t & bin_idx)
{
    // we need to construct a low level ibf AND insert all kmers into the high level bin
    assert(record.bins != 0);
    seqan3::interleaved_bloom_filter low_level{seqan3::bin_count{record.bins},
                                               seqan3::bin_size{8192u/*todo*/},
                                               seqan3::hash_function_count{2}};

    auto && info = read_sequences(record.filenames); // into random access containers

    std::ifstream fin{config.traversal_path_prefix + record.bin_name + ".out"};
    std::string line;

    if (!fin.good() || !fin.is_open())
        throw std::logic_error{"Could not open file '" + config.traversal_path_prefix + record.bin_name + ".out' for reading"};

    std::getline(fin, line); // skip header
    while (std::getline(fin, line))
    {
        auto && [filename, id, begin, end, idx] = parse_traversal_file_line(line);
        assert(bin_idx < high_level_ibf.bin_count());
        assert(idx < low_level.bin_count());

        for (auto hash : hash_infix(config, info[filename][id], begin, end))
        {
            high_level_ibf.emplace(hash, seqan3::bin_index{bin_idx});
            low_level.emplace(hash, seqan3::bin_index{idx});
        }
    }

    low_level_ibfs.emplace_back(low_level);
}

auto create_ibfs_from_data_file(build_config const & config)
{
    // the data file records, e.g. {bin_name, filenames, number_of_technical_bins}
    auto records = read_data_file_and_set_high_level_bins(config);

    size_t high_level_ibf_num_technical_bins{};
    data_file_record highest_record{};

    for (auto const & record : records)
    {
        high_level_ibf_num_technical_bins += (starts_with(record.bin_name, split_bin_prefix)) ? record.bins : 1;
        if (highest_record.max_size < record.max_size)
            highest_record = record;
    }

    // since the packing is supposedly done with minimizers and we probably have a different setting now
    // we need to calculate the maximum bin size. Since the relative number of minimizers should always correlate
    // wen can estimate it by calculating the kmer content of the "highest bin".
    auto max_bin_size = compute_maximum_technical_bin_size(config, highest_record);

    assert(high_level_ibf_num_technical_bins != 0);
    seqan3::interleaved_bloom_filter high_level_ibf{seqan3::bin_count{high_level_ibf_num_technical_bins},
                                                    seqan3::bin_size{max_bin_size},
                                                    seqan3::hash_function_count{config.hash_funs}};

    std::vector<seqan3::interleaved_bloom_filter<>> low_level_ibfs;

    size_t bin_idx{};
    for (auto const & record : records)
    {
        if (starts_with(record.bin_name, split_bin_prefix))
        {
            process_splitted_bin(config, record, high_level_ibf, bin_idx);
        }
        else if (starts_with(record.bin_name, merged_bin_prefix))
        {
            process_merged_bin(config, record, high_level_ibf, low_level_ibfs, bin_idx);
        }

        ++bin_idx;
    }

    return std::make_pair(std::move(high_level_ibf), std::move(low_level_ibfs));
}
