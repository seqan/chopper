#include <gtest/gtest.h>

#include <cmath>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/sketch/compute_sketches.hpp>

#include "../api_test.hpp"

struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<dna4_traits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

TEST(compute_sketches_test, small_example)
{
    std::vector<std::vector<uint64_t>> hashes(1);

    sequence_file_type fin{data("small.fa")};
    for (auto && [seq] : fin)
    {
        for (auto hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{15}))
            hashes[0].push_back(hash_value);
    }

    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::compute_sketches_from_hashes(hashes, 12, 1, sketches);

    ASSERT_EQ(sketches.size(), 1);
    EXPECT_EQ(std::lround(sketches[0].estimate()), 571);
}

TEST(compute_sketches_test, small_example_parallel_2_threads)
{
    std::vector<std::vector<uint64_t>> hashes(4);

    sequence_file_type fin{data("small.fa")};
    for (auto && [seq] : fin)
    {
        for (auto hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{15}))
        {
            hashes[0].push_back(hash_value);
            hashes[1].push_back(hash_value);
            hashes[2].push_back(hash_value);
            hashes[3].push_back(hash_value);
        }
    }

    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::compute_sketches_from_hashes(hashes, 12, 2, sketches);

    ASSERT_EQ(sketches.size(), 4);
    EXPECT_EQ(std::lround(sketches[0].estimate()), 571);
    EXPECT_EQ(std::lround(sketches[1].estimate()), 571);
    EXPECT_EQ(std::lround(sketches[2].estimate()), 571);
    EXPECT_EQ(std::lround(sketches[3].estimate()), 571);
}
