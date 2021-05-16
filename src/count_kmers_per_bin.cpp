#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/build/read_chopper_split_file.hpp>

struct cmd_arguments
{
    std::filesystem::path chopper_split_filename{};
    uint8_t k{25};
    size_t overlap{250};
    bool verbose;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Count unique kmers in bins.";
    parser.info.version = "1.0.0";

    parser.add_option(args.chopper_split_filename, 'f', "files", "Give me a file produced by chopper split.",
                      seqan3::option_spec::required);
    parser.add_option(args.k, 'k', "kmer-size", "The kmer to count with.");
    parser.add_option(args.overlap, 'l', "overlap", "The overlap between splitted bins.");
    parser.add_flag(args.verbose, 'v', "verbose", "Display more information.");
}

struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using seq_file_type = seqan3::sequence_file_input<file_type_traits,
                                                  seqan3::fields<seqan3::field::seq, seqan3::field::id>,
                                                  seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

auto read_sequences(std::vector<std::string> const & filenames)
{
    std::unordered_map<std::string, seqan3::dna4_vector> info;

    for (auto const & filename : filenames)
        for (auto && [seq, id] : seq_file_type{filename})
            info.emplace(filename + id, std::move(seq));

    return info;
}

auto hash_infix(cmd_arguments const & args, auto const & seq, auto const begin, auto const end)
{
    return seq | seqan3::views::drop(begin)
               | seqan3::views::take(end + args.overlap - begin) // views::take never goes over the end
               | seqan3::views::kmer_hash(seqan3::ungapped{args.k});
};

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

    auto [data, batches] = read_chopper_split_file(args.chopper_split_filename);

    std::vector<std::unordered_set<uint64_t>> kmer_counts;

    for (auto const & batch_record : batches)
    {
        auto && info = read_sequences(batch_record.filenames);

        if (batch_record.libf_num_bins > 0) // merged bin
        {
            kmer_counts.resize(batch_record.libf_num_bins);
            for (auto const & [combined_id, seq] : info)
                for (auto const & reg : data.region_map.at(combined_id))
                    for (auto hash : hash_infix(args, seq, reg.begin, reg.end))
                        kmer_counts[reg.lidx].insert(hash);

            for (size_t i = 0; i < batch_record.libf_num_bins; ++i)
                std::cout << batch_record.hibf_bins[0] << '_' << i << '\t' << kmer_counts[i].size() << '\n';

            // also output size of merged bin
            std::unordered_set<uint64_t> hibf_kmer_count{};
            for (auto & set : kmer_counts)
            {
                hibf_kmer_count.merge(set);
                set.clear();
            }

            std::cout << batch_record.hibf_bins[0] << '\t' << hibf_kmer_count.size();
        }
        else // split bin
        {
            kmer_counts.resize(batch_record.hibf_bins.size());
            for (auto const & [combined_id, seq] : info)
                for (auto const & reg : data.region_map.at(combined_id))
                    for (auto hash : hash_infix(args, seq, reg.begin, reg.end))
                        kmer_counts[reg.hidx].insert(hash);

            for (size_t i = 0; i < batch_record.hibf_bins.size(); ++i)
                std::cout << batch_record.hibf_bins[i] << '\t' << kmer_counts[i].size() << '\n';
        }

        kmer_counts.clear();
    }
}
