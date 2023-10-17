// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cinttypes>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <sharg/parser.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <hibf/sketch/hyperloglog.hpp>

struct cli_args
{
    std::filesystem::path input_path{};
    std::filesystem::path output_path{};
    std::vector<uint8_t> hll_bits;
    uint8_t k{20};
};

struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

int main(int argc, char const * argv[])
{
    sharg::parser parser{"measure_hyperloglog", argc, argv, sharg::update_notifications::off};

    cli_args args{};

    parser.add_option(args.input_path,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input",
                                    .description = "Fasta formatted file with sequences.",
                                    .required = true});
    parser.add_option(args.output_path,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output-file",
                                    .description = "File where the output is written to in tsv format.",
                                    .required = true});
    parser.add_option(
        args.hll_bits,
        sharg::config{.short_id = 'b',
                      .long_id = "hll-bits",
                      .description = "Adds an integer value which is tested as HyperLogLog bits parameter."});
    parser.add_option(args.k,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer",
                                    .description = "The size of the k-mers of which the hash values are computed."});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[COMMAND LINE INPUT ERROR] " << ext.what() << std::endl;
        return -1;
    }

    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_file_type = seqan3::sequence_file_input<input_traits, fields>;

    sequence_file_type seq_file{args.input_path};
    std::ofstream fout{args.output_path};

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    std::unordered_set<uint64_t> control;

    // write metadata
    fout << "# input file: " << args.input_path << '\n';
    fout << "# k: " << static_cast<int>(args.k) << '\n';

    // write header
    fout << "sequence_id\tsequence_length\tsketch_register_size\testimated_cardinality\tactual_cardinality\t"
         << "expected_relative_error\tactual_relative_error\n";

    for (uint8_t bits : args.hll_bits)
        sketches.emplace_back(bits);

    for (auto && [seq, id] : seq_file)
    {
        for (uint64_t && hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{args.k}))
        {
            control.insert(hash);
            for (auto & sketch : sketches)
            {
                sketch.add(hash);
            }
        }

        for (auto & sketch : sketches)
        {
            double const expected_error = 1.04 / std::sqrt(sketch.data_size());
            double const actual_error =
                std::abs(1.0 - std::round(sketch.estimate()) / static_cast<double>(control.size()));

            fout << id << '\t' << seq.size() << '\t' << sketch.data_size() << '\t'
                 << static_cast<uint64_t>(sketch.estimate()) << '\t' << control.size() << '\t' << expected_error << '\t'
                 << actual_error << '\n';
        }

        // clear for the next sequence
        for (auto & sketch : sketches)
            sketch.reset();

        control.clear();
    }
}
