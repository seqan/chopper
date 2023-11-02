// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>
#include <cstddef>
#include <string>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>

#include <hibf/config.hpp>

namespace chopper
{

struct input_functor
{
    struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
    };

    using sequence_file_type =
        seqan3::sequence_file_input<dna4_traits,
                                    seqan3::fields<seqan3::field::seq>,
                                    seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

    std::vector<std::string> filenames;

    bool input_are_precomputed_files{false};

    uint8_t kmer_size{21};

    uint8_t window_size{21};

    void operator()(size_t const num, seqan::hibf::insert_iterator it);
};

} // namespace chopper
