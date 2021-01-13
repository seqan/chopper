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

TEST_F(create_ibfs_from_chopper_pack_test, small_example)
{
    std::string seq1_filename = DATADIR"seq1.fa";
    std::string seq2_filename = DATADIR"seq2.fa";
    std::string seq3_filename = DATADIR"seq3.fa";

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_6\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             << "SPLIT_BIN_0\t" << seq1_filename << ";" << seq2_filename << "\t1\t500\n"
             << "SPLIT_BIN_1\t" << seq3_filename << "\t1\t500\n"
             << "SPLIT_BIN_2\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t1\t500\n"
             << "SPLIT_BIN_3\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t500\n"
             << "MERGED_BIN_6_0\t" << seq1_filename << "\t1\t500\n"
             << "MERGED_BIN_6_1\t" << seq2_filename << "\t1\t500\n"
             << "MERGED_BIN_6_2\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t500\n"
             << "MERGED_BIN_6_5\t" << seq3_filename << "\t1\t500\n";
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

    auto && [high_level_ibf, low_level_ibfs] = create_ibfs_from_chopper_pack(config);

    EXPECT_EQ(low_level_ibfs.size(), 7u);

    for (size_t i = 0; i < low_level_ibfs.size(); ++i)
        if (i != 6)
            EXPECT_EQ(low_level_ibfs[i].bin_size(), 1u); // dummy bin

    EXPECT_EQ(high_level_ibf.bin_size(), 114226);
    EXPECT_EQ(low_level_ibfs[6].bin_size(), 76615);

    EXPECT_EQ(high_level_ibf.bin_count(), 7u);
    EXPECT_EQ(low_level_ibfs[6].bin_count(), 6u);

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    auto hibf_agent = high_level_ibf.membership_agent();
    auto libf_agent = low_level_ibfs[6].membership_agent();

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

TEST_F(create_ibfs_from_chopper_pack_test, same_example_but_split_bin_as_hibf_max_bin)
{
    std::string seq1_filename = DATADIR"seq1.fa";
    std::string seq2_filename = DATADIR"seq2.fa";
    std::string seq3_filename = DATADIR"seq3.fa";

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_2\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             << "SPLIT_BIN_0\t" << seq1_filename << ";" << seq2_filename << "\t1\t500\n"
             << "SPLIT_BIN_1\t" << seq3_filename << "\t1\t500\n"
             << "SPLIT_BIN_2\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t1\t500\n"
             << "SPLIT_BIN_3\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t500\n"
             << "MERGED_BIN_6_0\t" << seq1_filename << "\t1\t500\n"
             << "MERGED_BIN_6_1\t" << seq2_filename << "\t1\t500\n"
             << "MERGED_BIN_6_2\t" << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t500\n"
             << "MERGED_BIN_6_5\t" << seq3_filename << "\t1\t500\n";
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

    auto && [high_level_ibf, low_level_ibfs] = create_ibfs_from_chopper_pack(config);

    EXPECT_EQ(low_level_ibfs.size(), 7u);

    for (size_t i = 0; i < low_level_ibfs.size(); ++i)
        if (i != 6)
            EXPECT_EQ(low_level_ibfs[i].bin_size(), 1u); // dummy bin

    EXPECT_EQ(high_level_ibf.bin_size(), 114226);
    EXPECT_EQ(low_level_ibfs[6].bin_size(), 76615);

    EXPECT_EQ(high_level_ibf.bin_count(), 7u);
    EXPECT_EQ(low_level_ibfs[6].bin_count(), 6u);

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    auto hibf_agent = high_level_ibf.membership_agent();
    auto libf_agent = low_level_ibfs[6].membership_agent();

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
