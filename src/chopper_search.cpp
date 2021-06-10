#include <fstream>
#include <seqan3/std/ranges>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/views/interleave.hpp>

#include <chopper/search/search.hpp>
#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>
#include <chopper/search/sync_out.hpp>

void initialize_argument_parser(seqan3::argument_parser & parser, search_config & config)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Read an HIBF on results from chopper-build and search queries in it.";
    parser.info.version = "1.0.0";

    parser.add_option(config.chopper_index_filename, 'i', "index", "Provide the HIBF index file produced by chopper build.");
    parser.add_option(config.errors, 'e', "errors", "The errors to allow in the search.");
    parser.add_option(config.query_filename, 'q', "queries", "The query sequences to seach for in the index.");
    parser.add_option(config.output_filename, 'o', "output", "The file to write results to.");
    parser.add_option(config.threads, 't', "threads", "The number of threads to use.");
    parser.add_flag(config.verbose, 'v', "verbose", "Output logging/progress information.");
    parser.add_flag(config.write_time, '\0', "time", "Write timing file.", seqan3::option_spec::advanced);
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

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

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
    sync_out sync_file{config.output_filename};

    auto cereal_worker = [&] ()
    {
        std::ifstream is{config.chopper_index_filename, std::ios::binary};

        if (!is.good() && !is.is_open())
            throw std::runtime_error{"File " + config.chopper_index_filename + " could not be opened"};

        cereal::BinaryInputArchive iarchive{is};

        auto start = std::chrono::high_resolution_clock::now();
        iarchive(data.hibf);
        iarchive(data.hibf_bin_levels);
        iarchive(data.user_bins);
        iarchive(config.k);
        auto end = std::chrono::high_resolution_clock::now();

        ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        write_header(data, sync_file);
    };

    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    search_sequence_file_t fin{config.query_filename};
    std::vector<search_sequence_file_t::record_type> records{};

    auto worker = [&] (auto && chunked_view, auto &&)
    {
        std::vector<size_t> read_kmers;
        std::vector<std::pair<int32_t, uint32_t>> result{};
        std::string buffer{};

        for (auto && [id, seq] : chunked_view)
        {
            clear_and_compute_kmers(read_kmers, seq, config);
            result.clear();

            search(result, read_kmers, data, config, 0); // start at top level ibf

            write_result(buffer, result, id, data, sync_file);
        }
    };

    for (auto && record_batch : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        records.clear();

        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(record_batch, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        size_t const records_per_thread = records.size() / config.threads;
        seqan3::detail::execution_handler_parallel executioner{config.threads};
        auto chunked_records = records | seqan3::views::chunk(records_per_thread);

        cereal_handle.wait();

        start = std::chrono::high_resolution_clock::now();
        executioner.bulk_execute(std::move(worker), std::move(chunked_records), [](){});
        end = std::chrono::high_resolution_clock::now();
        compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    }

    if (config.write_time)
    {
        std::filesystem::path file_path{config.output_filename};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }

    return 0;
}

