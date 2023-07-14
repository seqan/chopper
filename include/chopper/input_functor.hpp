#pragma once

#include <string>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

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

    void operator()(size_t const num, hibf::insert_iterator it)
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

            for (auto && [seq] : fin)
            {
                for (auto hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{kmer_size}))
                    it = hash_value;
            }
        }
    }
};

} // namespace chopper
