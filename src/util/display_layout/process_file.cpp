// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cinttypes>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <chopper/adjust_seed.hpp>

#include "shared.hpp"

struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<dna4_traits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

void process_file(std::string const & filename,
                  std::vector<uint64_t> & current_kmers,
                  seqan::hibf::sketch::hyperloglog & sketch,
                  bool const fill_current_kmers,
                  uint8_t const kmer_size)
{
    if (filename.ends_with(".minimiser"))
    {
        uint64_t hash{};
        char * const hash_data{reinterpret_cast<char *>(&hash)};
        std::streamsize const hash_bytes{sizeof(hash)};

        std::ifstream infile{filename, std::ios::binary};

        if (fill_current_kmers)
        {
            while (infile.read(hash_data, hash_bytes))
            {
                current_kmers.push_back(hash);
                sketch.add(hash);
            }
        }
        else
        {
            while (infile.read(hash_data, hash_bytes))
            {
                sketch.add(hash);
            }
        }
    }
    else
    {
        sequence_file_type fin{filename};

        seqan3::shape shape{seqan3::ungapped{kmer_size}};
        auto minimizer_view = seqan3::views::minimiser_hash(shape,
                                                            seqan3::window_size{kmer_size},
                                                            seqan3::seed{chopper::adjust_seed(shape.count())});
        if (fill_current_kmers)
        {
            for (auto && [seq] : fin)
            {
                for (uint64_t hash_value : seq | minimizer_view)
                {
                    current_kmers.push_back(hash_value);
                    sketch.add(hash_value);
                }
            }
        }
        else
        {
            for (auto && [seq] : fin)
            {
                for (uint64_t hash_value : seq | minimizer_view)
                {
                    sketch.add(hash_value);
                }
            }
        }
    }
}
