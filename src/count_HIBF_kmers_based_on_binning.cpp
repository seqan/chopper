#include <future>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/build/build_config.hpp>
#include <chopper/build/read_chopper_pack_file.hpp>
#include <chopper/pack/filenames_data_input.hpp>
#include <chopper/pack/configuration.hpp>
#include <chopper/pack/data_store.hpp>
#include <chopper/print_peak_memory_usage.hpp>

#include <robin_hood.h>

struct cmd_arguments
{
    std::filesystem::path binning_file{};
    std::filesystem::path count_file{};
    std::filesystem::path output_file{};
    uint8_t k{25};
    size_t threads{1u};
    bool verbose{false};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Count unique kmers in bins.";
    parser.info.version = "1.0.0";

    parser.add_option(args.binning_file, 'b', "binning", "Give me a binning file.", seqan3::option_spec::required);
    parser.add_option(args.count_file, 'c', "counts", "Give me a kmer counts file.", seqan3::option_spec::required);
    parser.add_option(args.k, 'k', "kmer-size", "The kmer size to count with.");
    parser.add_option(args.threads, 't', "threads", "The number of threads to use.");
    parser.add_option(args.output_file, 'o', "output", "Give me an output file.", seqan3::option_spec::required);
}

struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

std::mutex mutex;

auto print_kmer_content(chopper_pack_record const & record,
                        size_t const num_bins,
                        uint8_t const k,
                        robin_hood::unordered_map<std::string, size_t> const & counts,
                        std::ofstream & fout)
{
    using seq_file_type = seqan3::sequence_file_input<file_type_traits,
                                                      seqan3::fields<seqan3::field::seq>,
                                                      seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

    robin_hood::unordered_node_set<uint64_t> kmer_occurence{};

    size_t low_lvl_kmer_sum = 0;
    bool is_merged = starts_with(record.bin_name, merged_bin_prefix);

    for (auto const & filename : record.filenames)
    {
        if (is_merged)
            low_lvl_kmer_sum += counts.at(filename);

        for (auto && [seq] : seq_file_type{filename})
        {
            auto kmer_hash_view = seq | seqan3::views::kmer_hash(seqan3::ungapped{k});// | std::views::common;
            kmer_occurence.insert(kmer_hash_view.begin(), kmer_hash_view.end());
        }
    }

    std::string s;

    for (size_t i = 0; i < num_bins; ++i)
    {
        s += record.bin_name + '\t';
        s += std::to_string(kmer_occurence.size() / num_bins) + '\t';
        s += std::to_string(low_lvl_kmer_sum) + '\n';
    }

    {
        std::lock_guard<std::mutex> lock(mutex)
        fout << s;
    }
}

int main(int const argc, char const ** argv)
{
    seqan3::argument_parser myparser{"count_kmers_per_bin", argc, argv};
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[ERROR] " << ext.what() << "\n";
        return -1;
    }

    size_t const counting_threads = (args.threads <= 1) ? 1 : args.threads - 1;

    // output file
    std::ofstream fout{args.output_file};

    build_config config{};
    build_data data{};
    config.binning_filename = args.binning_file;

    // std::vector<chopper_pack_record> records{};
    // std::ifstream chopper_pack_file{config.binning_filename};
    // std::string line{};
    // std::getline(chopper_pack_file, line);
    // for (; line.substr(0, 6) != "#FILES"; std::getline(chopper_pack_file, line)) {}
    // while (std::getline(chopper_pack_file, line))
    //     records.push_back(parse_chopper_pack_line(line);

    read_chopper_pack_file(data, config.binning_filename);

    // Reuse the kmer counts file.
    configuration p_config;
    p_config.data_file = args.count_file;

    data_store data;
    read_filename_data_file(data, p_config);

    robin_hood::unordered_map<std::string, size_t> counts;

    // make a simple lookup table for kmer counts
    for (size_t i = 0; i < data.filenames.size(); ++i)
        counts[data.filenames[i]] = data.kmer_counts[i];

    // setup async execution
    auto async_buf_records = records | seqan3::views::async_input_buffer(counting_threads);

    auto worker = [&] ()
    {
        for (auto const & record : async_buf_records)
        {
            if (starts_with(record.bin_name, split_bin_prefix))
            {
                print_kmer_content(record, record.bins, args.k, counts, fout);
            }
            else if (starts_with(record.bin_name, merged_bin_prefix))
            {
                print_kmer_content(record, 1, args.k, counts, fout); // always one bin in high-level IBF, record.bins is for the low-level IBF
            }
        }
    };

    // launch threads with worker
    std::vector<decltype(std::async(std::launch::async, worker))> handles;

    for (size_t i = 0; i < counting_threads; ++i)
        handles.emplace_back(std::async(std::launch::async, worker));

    // wait for all threads to finish
    for (auto && handle : handles)
        handle.get();

    print_peak_memory_usage();

    return 0;
}
