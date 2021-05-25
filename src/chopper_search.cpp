#include <fstream>
#include <seqan3/std/ranges>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/views/interleave.hpp>

#include <chopper/search/search.hpp>
#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, search_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Read an HIBF on results from chopper-build and search queries in it.";
    parser.info.version = "1.0.0";

    // todo: k-mer size should be serialized with building the index to avoid inconsistencies
    parser.add_option(config.chopper_index_filename, 'i', "index", "Provide the HIBF index file produced by chopper build.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to build kmers.");
    parser.add_option(config.errors, 'e', "errors", "The errors to allow in the search.");
    parser.add_option(config.query_filename, 'q', "queries", "The query sequences to seach for in the index.");
    parser.add_flag(config.verbose, 'v', "verbose", "Output logging/progress information.");
}

struct search_file_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using search_sequence_file_t = seqan3::sequence_file_input<search_file_traits,
                                                    seqan3::fields<seqan3::field::id, seqan3::field::seq>,
                                                    seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

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

        if (!is.good() && !is.is_open())
            throw std::runtime_error{"File " + config.chopper_index_filename + " could not be opened"};

        cereal::BinaryInputArchive iarchive{is};
        iarchive(data.hibf);
        iarchive(data.hibf_bin_levels);
        iarchive(data.user_bins);
    }

    write_header(data, std::cout);

    std::vector<size_t> read_kmers;
    std::vector<std::pair<int32_t, uint32_t>> result{};

    for (auto && [id, seq] : search_sequence_file_t{config.query_filename})
    {
        clear_and_compute_kmers(read_kmers, seq, config);
        result.clear();

        search(result, read_kmers, data, config, 0); // start at top level ibf

        write_result(result, id, data, std::cout);
    }

    return 0;
}

