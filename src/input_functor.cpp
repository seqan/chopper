// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <cinttypes>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <chopper/adjust_seed.hpp>
#include <chopper/input_functor.hpp>

namespace chopper
{

void input_functor::operator()(size_t const num, seqan::hibf::insert_iterator it)
{
    assert(filenames.size() > num);
    if (input_are_precomputed_files)
    {
        uint64_t hash{};
        char * const hash_data{reinterpret_cast<char *>(&hash)};
        std::streamsize const hash_bytes{sizeof(hash)};

        std::ifstream infile{filenames[num], std::ios::binary};

        while (infile.read(hash_data, hash_bytes))
            it = hash;
    }
    else
    {
        sequence_file_type fin{filenames[num]};

        seqan3::shape shape = seqan3::ungapped{kmer_size};
        auto minimizer_view = seqan3::views::minimiser_hash(shape,
                                                            seqan3::window_size{window_size},
                                                            seqan3::seed{adjust_seed(shape.count())});

        for (auto && [seq] : fin)
        {
            for (auto hash_value : seq | minimizer_view)
                it = hash_value;
        }
    }
}

} // namespace chopper
