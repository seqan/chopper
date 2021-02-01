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
    auto count_kmers(typename seqan3::interleaved_bloom_filter<>::membership_agent & ibf_agent,
                     seqan3::dna4_vector const & query,
                     build_config const & config)
    {
        std::vector<size_t> ibf_counts(ibf_agent.result_buffer.size());

        for (auto hash : query | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        {
            auto const & result = ibf_agent.bulk_contains(hash);

            for (size_t i = 0; i < result.size(); ++i)
                ibf_counts[i] += result[i];
        }

        return ibf_counts;
    }

    auto compare_counts(typename seqan3::interleaved_bloom_filter<> & ibf,
                        seqan3::dna4_vector const & query,
                        std::vector<std::vector<size_t>> const & bins_with_expected_counts,
                        build_config const & config)
    {
        auto agent = ibf.membership_agent();
        auto && counts = this->count_kmers(agent, query, config);

        size_t expected = std::ranges::distance(query | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        std::vector<bool> without_expected_counts(counts.size(), true);

        for (auto const & bins : bins_with_expected_counts)
        {
            size_t sum{};
            for (auto const bin : bins)
            {
                without_expected_counts[bin] = false;
                sum += counts[bin];
            }
            EXPECT_EQ(sum, expected);
        }

        for (size_t i = 0; i < counts.size(); ++i)
            if (without_expected_counts[i])
                EXPECT_EQ(counts[i], 0) << " failed for i=" << i << std::endl;
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

    // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
    compare_counts(high_level_ibf, unspecific, {{0}, {1}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, unspecific, {{0}, {1}, {2, 3, 4}, {5}},config);

    // SEQ2 SPECIFIC
    compare_counts(high_level_ibf, seq2_specific, {{0}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, seq2_specific, {{1}, {2, 3, 4}},config);

    // SEQ3 SPECIFIC
    compare_counts(high_level_ibf, seq3_specific, {{1}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, seq3_specific, {{2, 3, 4}, {5}},config);
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

    // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
    compare_counts(high_level_ibf, unspecific, {{0}, {1}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, unspecific, {{0}, {1}, {2, 3, 4}, {5}},config);

    // SEQ2 SPECIFIC
    compare_counts(high_level_ibf, seq2_specific, {{0}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, seq2_specific, {{1}, {2, 3, 4}},config);

    // SEQ3 SPECIFIC
    compare_counts(high_level_ibf, seq3_specific, {{1}, {2}, {3,4,5}, {6}},config);
    compare_counts(low_level_ibf, seq3_specific, {{2, 3, 4}, {5}},config);
}

TEST_F(create_ibfs_from_chopper_pack_test, multi_level_ibf)
{
    std::string seq1_filename = DATADIR"seq1.fa";
    std::string seq2_filename = DATADIR"seq2.fa";
    std::string seq3_filename = DATADIR"seq3.fa";
    std::string seq12_filename = seq1_filename + ";" + seq2_filename;
    std::string seq13_filename = seq1_filename + ";" + seq3_filename;
    std::string seq23_filename = seq2_filename + ";" + seq3_filename;
    std::string seq123_filename = seq1_filename + ";" + seq2_filename + ";" + seq3_filename;

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:0\n"
             << "#MERGED_BIN_0;0;0 max_bin_id:3\n"
             << "#MERGED_BIN_0;0 max_bin_id:4\n"
             << "#MERGED_BIN_0;1 max_bin_id:2\n"
             << "#MERGED_BIN_0 max_bin_id:0\n"
             << "#MERGED_BIN_1 max_bin_id:2\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << seq1_filename   << "\t0;0;0;0\t1;1;1;3\t1650;350;80;1\n"
             << seq2_filename   << "\t0;0;0;3\t1;1;1;1\t1650;350;80;2\n"
             << seq3_filename   << "\t0;0;0;4\t1;1;1;2\t1650;350;80;2\n"
             << seq23_filename  << "\t0;0;0;6\t1;1;1;2\t1650;350;80;2\n"
             << seq123_filename << "\t0;0;0;8\t1;1;1;5\t1650;350;80;2\n"
             << seq1_filename   << "\t0;0;1\t1;1;1\t1650;350;40\n"
             << seq2_filename   << "\t0;0;2\t1;1;1\t1650;350;50\n"
             << seq3_filename   << "\t0;0;3\t1;1;1\t1650;350;80\n"
             << seq12_filename  << "\t0;0;4\t1;1;1\t1650;350;100\n"
             << seq1_filename   << "\t0;1;0\t1;1;1\t1650;400;6\n"
             << seq2_filename   << "\t0;1;1\t1;1;1\t1650;400;7\n"
             << seq3_filename   << "\t0;1;2\t1;1;1\t1650;400;7\n"
             << seq12_filename  << "\t0;2\t1;1\t1650;200\n"
             << seq13_filename  << "\t0;3\t1;1\t1650;300\n"
             << seq23_filename  << "\t0;4\t1;1\t1650;400\n"
             << seq1_filename   << "\t1;0\t1;1\t2000;31\n"
             << seq2_filename   << "\t1;1\t1;1\t2000;32\n"
             << seq12_filename  << "\t1;2\t1;2\t2000;32\n"
             << seq23_filename  << "\t2\t1\t1200\n"
             << seq123_filename << "\t3\t2\t1500\n";
    }

    build_config config{};
    config.k = 15;
    config.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.ibfs.size(), 6);
    // data.ibfs[0] is the high-level IBF
    // data.ibfs[1] is the LOW LEVEL IBF 0;0;0
    // data.ibfs[2] is the LOW LEVEL IBF 0;0
    // data.ibfs[3] is the LOW LEVEL IBF 0;1
    // data.ibfs[4] is the LOW LEVEL IBF 0
    // data.ibfs[5] is the LOW LEVEL IBF 1

    /* HIGH LEVEL IBF
     * --------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2 (merged)
     * Bin 2: seq2, seq3
     * Bin 3: split but together: seq1, seq2, seq3
     * Bin 4: --> belongs to bin 3
     */
    EXPECT_EQ(data.ibfs[0].bin_count(), 5u);
    EXPECT_EQ(data.ibfs[0].bin_size(), 114226);

    /* LOW LEVEL IBF 0;0;0
     * -------------------
     * Bin 0: split but together: seq1
     * Bin 1: --> belongs to bin 0
     * Bin 2: --> belongs to bin 0
     * Bin 3: seq2
     * Bin 4: split but together: seq3
     * Bin 5: --> belongs to bin 4
     * Bin 6: split but together: seq2, seq3
     * Bin 7: --> belongs to bin 6
     * Bin 8: split but together: seq1, seq2, seq3
     * Bin 9: --> belongs to bin 8
     * Bin 10: --> belongs to bin 8
     * Bin 11: --> belongs to bin 8
     * Bin 12: --> belongs to bin 8
     */
    EXPECT_EQ(data.ibfs[1].bin_count(), 13u);
    EXPECT_EQ(data.ibfs[1].bin_size(), 92535);

    /* LOW LEVEL IBF 0;0
     * -----------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1
     * Bin 2: seq2
     * Bin 3: seq3
     * Bin 4: seq1, seq2
     */
    EXPECT_EQ(data.ibfs[2].bin_count(), 5u);
    EXPECT_EQ(data.ibfs[2].bin_size(), 95321);

    /* LOW LEVEL IBF 0;1
     * -----------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: seq3
     */
    EXPECT_EQ(data.ibfs[3].bin_count(), 3u);
    EXPECT_EQ(data.ibfs[3].bin_size(), 92734);

    /* LOW LEVEL IBF 0
     * ---------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2, seq3 (merged)
     * Bin 2: seq1, seq2
     * Bin 3: seq1, seq3
     * Bin 4: seq2, seq3
     */
    EXPECT_EQ(data.ibfs[4].bin_count(), 5u);
    EXPECT_EQ(data.ibfs[4].bin_size(), 114226);

    /* LOW LEVEL IBF 1
     * ---------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: split but together: seq1, seq2
     * Bin 3: --> belongs to bin 2
     */
    EXPECT_EQ(data.ibfs[5].bin_count(), 4u);
    EXPECT_EQ(data.ibfs[5].bin_size(), 47561);

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
    compare_counts(data.ibfs[0], unspecific, {{0}, {1}, {2}, {3,4}}, config);
    compare_counts(data.ibfs[1], unspecific, {{0, 1, 2}, {3}, {4,5}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.ibfs[2], unspecific, {{0}, {1}, {2}, {3}, {4}}, config);
    compare_counts(data.ibfs[3], unspecific, {{0}, {1}, {2}}, config);
    compare_counts(data.ibfs[4], unspecific, {{0}, {1}, {2}, {3}, {4}}, config);
    compare_counts(data.ibfs[5], unspecific, {{0}, {1}, {2,3}}, config);

    // SEQ2 SPECIFIC
    compare_counts(data.ibfs[0], seq2_specific, {{0}, {1}, {2}, {3,4}}, config);
    compare_counts(data.ibfs[1], seq2_specific, {{3}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.ibfs[2], seq2_specific, {{0}, {2}, {4}}, config);
    compare_counts(data.ibfs[3], seq2_specific, {{1}}, config);
    compare_counts(data.ibfs[4], seq2_specific, {{0}, {1}, {2}, {4}}, config);
    compare_counts(data.ibfs[5], seq2_specific, {{1}, {2,3}}, config);

    // SEQ3 SPECIFIC
    compare_counts(data.ibfs[0], seq3_specific, {{0}, {2}, {3,4}}, config);
    compare_counts(data.ibfs[1], seq3_specific, {{4,5}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.ibfs[2], seq3_specific, {{0}, {3}}, config);
    compare_counts(data.ibfs[3], seq3_specific, {{2}}, config);
    compare_counts(data.ibfs[4], seq3_specific, {{0}, {1}, {3}, {4}}, config);
    compare_counts(data.ibfs[5], seq3_specific, {}, config);
}
