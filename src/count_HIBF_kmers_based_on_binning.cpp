#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/build/read_data_file_and_set_high_level_bins.hpp>

struct cmd_arguments
{
    std::filesystem::path binning_file{};
    uint8_t k{25};
    bool verbose;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Count unique kmers in bins.";
    parser.info.version = "1.0.0";

    parser.add_option(args.binning_file, 'f', "files", "Give me a file.", seqan3::option_spec::required);
    parser.add_option(args.k, 'k', "kmer-size", "The kmer to count with.");
}

struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

auto print_kmer_content(chopper_pack_record const & record, size_t const num_bins, uint8_t const k)
{
    using seq_file_type = seqan3::sequence_file_input<file_type_traits,
                                                      seqan3::fields<seqan3::field::seq>,
                                                      seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

    std::unordered_set<uint64_t> kmer_occurence{};

    for (auto const & filename : record.filenames)
        for (auto && [seq] : seq_file_type{filename})
            for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{k}))
                kmer_occurence.insert(hash);

    for (size_t i = 0; i < num_bins; ++i)
        seqan3::debug_stream << record.bin_name << '\t' << kmer_occurence.size()/num_bins << std::endl;
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

    build_config config;
    config.binning_filename = args.binning_file;
    auto [header, records] = read_data_file_and_set_high_level_bins(config);

    for (auto const & record : records)
    {
        std::unordered_set<uint64_t> kmer_occurence{};

        if (starts_with(record.bin_name, split_bin_prefix))
        {
            print_kmer_content(record, record.bins, args.k);
        }
        else if (starts_with(record.bin_name, merged_bin_prefix))
        {
            print_kmer_content(record, 1, args.k); // always one bin in high-level IBF, record.bins is for the low-level IBF
        }
    }
}
