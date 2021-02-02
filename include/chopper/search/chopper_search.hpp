#pragma once

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, search_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Read an HIBF on results from chopper-build and search queries in it.";
    parser.info.version = "1.0.0";

    // todo: k-mer size should be serialized with building the index to avoid inconsistencies
    parser.add_option(config.chopper_index_filename, 'i', "index", "Provide the HIBF index file produced by chopper build.");
    parser.add_option(config.chopper_index_map_filename, 'm', "index-map", "Provide the HIBF index mapping file produced by chopper build.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to build kmers.");
    parser.add_option(config.errors, 'e', "errors", "The errors to allow in the search.");
    parser.add_option(config.query_filename, 'q', "queries", "The query sequences to seach for in the index.");
    parser.add_flag(config.verbose, 'v', "verbose", "Output logging/progress information.");
}

struct file_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_t = seqan3::sequence_file_input<file_traits,
                                                    seqan3::fields<seqan3::field::id, seqan3::field::seq>,
                                                    seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

void search(std::unordered_set<std::pair<size_t, size_t>> & membership_result,
            std::vector<size_t> const & kmers,
            search_data const & data,
            size_t const ibf_idx)
{
    size_t const kmer_lemma = (kmers.size() > config.errors * config.k)
                              ? kmers.size() - config.errors * config.k
                              : 0;

    auto counting_agent = data.hibf[ibf_idx].template counting_agent<uint16_t>();

    auto && result = counter.count_hashes(kmers);

    for (size_t bin = 0; bin < result.size(); ++bin)
    {
        if (result[bin] > kmer_lemma)
        {
            size_t next_ibf_idx = data.ibf_mapping[ibf_idx];
            if (next_ibf != ibf_idx)
            {
                search(kmers, data, next_ibf_idx);
            }
            else
            {
                membership_result.emplace(ibf_idx, bin);
            }
        }
    }
}

std::vector<size_t> compute_kmers(dna4_vector const & query, search_config const & config)
{
    std::vector<size_t> kmers{};

    for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        kmers.insert(hash);

    return kmers;
}

void write_result(std::unordered_set<std::pair<size_t, size_t>> const & membership_result,
                  std::string const & id,
                  search_data const & data)
{
    std::cout << id << '\t';
    for (auto && [ibf_idx, bin_idx] : membership_result)
        std::cout << /*returns a number*/ data.user_bins[ibf_idx][bin_idx] << ','; // todo remove last comma
    std::cout << std::endl;
}

int chopper_search(seqan3::argument_parser & parser)
{
    search_config config{};
    initialize_argument_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[CHOPPER SEARCH ERROR] " << ext.what() << "\n";
        return -1;
    }

    search_data data;

    { // read ibf vector
        std::ifstream is{config.chopper_index_filename, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(data.hibf);
    }

    { // read ibf mapping
        std::ifstream is{config.chopper_index_map_filename, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(data.hibf_mapping);
    }

    write_header();

    std::vector<size_t> read_kmers; // allocate space once

    for (auto && [id, seq] : sequence_file_t{config.query_filename})
    {
        read_kmers = compute_kmers(seq, config);

        search(read_kmers, data, 0); // start at top level ibf

        write_result(membership_result, id, data);
    }

    return 0;
}

