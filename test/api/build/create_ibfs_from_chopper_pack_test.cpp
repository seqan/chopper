#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <chopper/build/create_ibfs_from_chopper_pack.hpp>

using seqan3::operator""_dna4;

struct create_ibfs_from_chopper_pack_test : public ::testing::Test
{
    auto count_kmers(typename seqan3::interleaved_bloom_filter<>::membership_agent & hibf_agent,
                     typename seqan3::interleaved_bloom_filter<>::membership_agent & libf_agent,
                     seqan3::dna4_vector const & query,
                     build_config const & config)
    {
        std::vector<size_t> hibf_counts(hibf_agent.result_buffer.size());
        std::vector<size_t> libf_counts(libf_agent.result_buffer.size());

        for (auto hash : query | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        {
            auto const & high_res_v = hibf_agent.bulk_contains(hash);
            auto const & low_res_v = libf_agent.bulk_contains(hash);

            for (size_t i = 0; i < high_res_v.size(); ++i)
                hibf_counts[i] += high_res_v[i];

            for (size_t i = 0; i < low_res_v.size(); ++i)
                libf_counts[i] += low_res_v[i];
        }

        return std::make_tuple(std::move(hibf_counts), std::move(libf_counts));
    }
};

TEST_F(create_ibfs_from_chopper_pack_test, small_example_2_levels)
{
    std::string seq1_filename = DATADIR"seq1.fa";
    std::string seq2_filename = DATADIR"seq2.fa";
    std::string seq3_filename = DATADIR"seq3.fa";

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:6\n"
             << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << seq1_filename << ";" << seq2_filename << "\t0\t1\t500\n"
             << seq3_filename << "\t1\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t2\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t3\t500\n"
             << seq1_filename << "\t6;0\t1;1\t500\n"
             << seq2_filename << "\t6;1\t1;1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t6;2\t1;3\t500\n"
             << seq3_filename << "\t6;5\t1;1\t500\n";
    }

    // HIGH LEVEL IBF
    // --------------
    // Bin 0: seq1, seq2
    // Bin 1: seq3
    // Bin 2: seq1, seq2, seq3
    // Bin 3: ? not easily to determine since kmers are split independently
    // Bin 4: ? not easily to determine since kmers are split independently
    // Bin 5: ? not easily to determine since kmers are split independently
    // Bin 6: seq1, seq2, seq3

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1
    // Bin 1: seq2
    // Bin 2: ? not easily to determine since kmers are split independently
    // Bin 3: ? not easily to determine since kmers are split independently
    // Bin 4: ? not easily to determine since kmers are split independently
    // Bin 5: seq3

    build_config config{};
    config.k = 15;
    config.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.ibfs.size(), 2);

    auto & high_level_ibf = data.ibfs[0];
    auto & low_level_ibf = data.ibfs[1];

    EXPECT_EQ(high_level_ibf.bin_size(), 114226);
    EXPECT_EQ(low_level_ibf.bin_size(), 76615);

    EXPECT_EQ(high_level_ibf.bin_count(), 7u);
    EXPECT_EQ(low_level_ibf.bin_count(), 6u);

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    auto hibf_agent = high_level_ibf.membership_agent();
    auto libf_agent = low_level_ibf.membership_agent();

    { // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, unspecific, config);

        size_t expected = std::ranges::distance(unspecific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], expected);
        EXPECT_EQ(hibf_counts[1], expected);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], expected);
        EXPECT_EQ(libf_counts[1], expected);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], expected);
    }

    { // SEQ2 SPECIFIC
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, seq2_specific, config);

        size_t expected = std::ranges::distance(seq2_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], expected);
        EXPECT_EQ(hibf_counts[1], 0);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], 0);
        EXPECT_EQ(libf_counts[1], expected);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], 0);
    }

    { // SEQ3 SPECIFIC
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, seq3_specific, config);

        size_t expected = std::ranges::distance(seq3_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], 0);
        EXPECT_EQ(hibf_counts[1], expected);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], 0);
        EXPECT_EQ(libf_counts[1], 0);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], expected);
    }
}

TEST_F(create_ibfs_from_chopper_pack_test, same_example_two_levels_but_split_bin_as_hibf_max_bin)
{
    std::string seq1_filename = DATADIR"seq1.fa";
    std::string seq2_filename = DATADIR"seq2.fa";
    std::string seq3_filename = DATADIR"seq3.fa";

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:2\n"
             << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << seq1_filename << ";" << seq2_filename << "\t0\t1\t500\n"
             << seq3_filename << "\t1\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t2\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t3\t500\n"
             << seq1_filename << "\t6;0\t1;1\t500\n"
             << seq2_filename << "\t6;1\t1;1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t6;2\t1;3\t500\n"
             << seq3_filename << "\t6;5\t1;1\t500\n";
    }

    // HIGH LEVEL IBF
    // --------------
    // Bin 0: seq1, seq2
    // Bin 1: seq3
    // Bin 2: seq1, seq2, seq3
    // Bin 3: ? not easily to determine since kmers are split independently
    // Bin 4: ? not easily to determine since kmers are split independently
    // Bin 5: ? not easily to determine since kmers are split independently
    // Bin 6: seq1, seq2, seq3

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1
    // Bin 1: seq2
    // Bin 2: ? not easily to determine since kmers are split independently
    // Bin 3: ? not easily to determine since kmers are split independently
    // Bin 4: ? not easily to determine since kmers are split independently
    // Bin 5: seq3

    build_config config{};
    config.k = 15;
    config.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.ibfs.size(), 2);

    auto & high_level_ibf = data.ibfs[0];
    auto & low_level_ibf = data.ibfs[1];

    EXPECT_EQ(high_level_ibf.bin_size(), 114226);
    EXPECT_EQ(low_level_ibf.bin_size(), 76615);

    EXPECT_EQ(high_level_ibf.bin_count(), 7u);
    EXPECT_EQ(low_level_ibf.bin_count(), 6u);

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    auto hibf_agent = high_level_ibf.membership_agent();
    auto libf_agent = low_level_ibf.membership_agent();

    { // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, unspecific, config);

        size_t expected = std::ranges::distance(unspecific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], expected);
        EXPECT_EQ(hibf_counts[1], expected);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], expected);
        EXPECT_EQ(libf_counts[1], expected);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], expected);
    }

    { // SEQ2 SPECIFIC
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, seq2_specific, config);

        size_t expected = std::ranges::distance(seq2_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], expected);
        EXPECT_EQ(hibf_counts[1], 0);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], 0);
        EXPECT_EQ(libf_counts[1], expected);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], 0);
    }

    { // SEQ3 SPECIFIC
        auto && [hibf_counts, libf_counts] = this->count_kmers(hibf_agent, libf_agent, seq3_specific, config);

        size_t expected = std::ranges::distance(seq3_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(hibf_counts[0], 0);
        EXPECT_EQ(hibf_counts[1], expected);
        EXPECT_EQ(hibf_counts[2], expected);
        EXPECT_EQ(hibf_counts[3] + hibf_counts[4] + hibf_counts[5], expected);
        EXPECT_EQ(hibf_counts[6], expected);

        EXPECT_EQ(libf_counts[0], 0);
        EXPECT_EQ(libf_counts[1], 0);
        EXPECT_EQ(libf_counts[2] + libf_counts[3] + libf_counts[4], expected);
        EXPECT_EQ(libf_counts[5], expected);
    }
}
