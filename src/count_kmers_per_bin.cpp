#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan3/argument_parser/all.hpp>

#include <chopper/build/parse_traversal_file_line.hpp>

auto read_traversal_file(std::filesystem::path const & traversal_filename)
{
    std::set<std::string> filenames;
    std::vector<std::vector<std::tuple<std::string, uint32_t, uint32_t>>> bin_ranges;

    std::ifstream fin{traversal_filename};

    if (!fin.good())
        throw std::runtime_error{"Could not open '" + traversal_filename.string() + "' for reading."};

    std::string line;
    std::getline(fin, line); // skip header line
    assert(line == std::string{"FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER"});

    while (std::getline(fin, line))
    {
        auto [file_name, seq_name, begin, end, bin] = parse_traversal_file_line(line);
        filenames.insert(file_name);
        if (bin >= bin_ranges.size())
            bin_ranges.resize(bin + 1);
        bin_ranges[bin].emplace_back(seq_name, begin, end);
    }

    return std::make_pair(std::move(filenames), std::move(bin_ranges));
}

struct cmd_arguments
{
    std::filesystem::path traversal_file{};
    size_t k{25};
    size_t overlap{250};
    bool verbose;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Avenja";
    parser.info.short_description = "Count unique kmers in bins.";
    parser.info.version = "1.0.0";

    parser.add_option(args.traversal_file, 'f', "files", "Give me a file.", seqan3::option_spec::REQUIRED);
    parser.add_option(args.k, 'k', "kmer-size", "The kmer to count with.");
    parser.add_option(args.overlap, 'l', "overlap", "The overlap between splitted bins.");
    parser.add_flag(args.verbose, 'v', "verbose", "Display more information.");
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

    auto [filenames, bin_ranges] = read_traversal_file(args.traversal_file);

    std::unordered_map<std::string, size_t> sequence_names_positions;

    seqan::String<seqan::String<seqan::Iupac>> seqs;

    for (auto filename : filenames)
    {
        seqan::SeqFileIn input_file(filename.c_str());

        while (!atEnd(input_file))
        {
            seqan::CharString id;
            seqan::String<seqan::Iupac> seq;
            seqan::readRecord(id, seq, input_file);
            appendValue(seqs, std::move(seq));
            sequence_names_positions[seqan::toCString(id)] = length(seqs) - 1;
        }
    }

    size_t total_sum_over_all_bins{};

    // count
    // -------------------------------------------------------------------------
    for (auto & current_bin_range : bin_ranges)
    {
        std::unordered_set<uint64_t> kmer_occurence{};

        for (auto const & [seq_id, begin, end] : current_bin_range)
        {
            auto const seq_pos = sequence_names_positions.at(seq_id);

            seqan::Shape<seqan::Iupac, seqan::SimpleShape> myShape;
            resize(myShape, args.k);
            seqan::hashInit(myShape, seqan::begin(seqs[seq_pos]));

            if (end < begin + seqan::length(myShape) + 1)
            {
                if (args.verbose)
                    std::cout << "[WARNING] invalid range: [" << begin << "," << end << "]" << std::endl;
                continue;
            }

            assert(end >= seqan::length(myShape) + 1);
            for (uint32_t i = begin; i < end - seqan::length(myShape) + 1; ++i)
            {
                auto hash = seqan::hashNext(myShape, seqan::begin(seqs[seq_pos]) + i);
                kmer_occurence.insert(hash);
            }
        }

        std::cout << kmer_occurence.size() << '\t';

        total_sum_over_all_bins += kmer_occurence.size();
    }

    if (args.verbose)
        std::cout << "--------------" << std::endl;

    std::cout << total_sum_over_all_bins << std::endl;
}
