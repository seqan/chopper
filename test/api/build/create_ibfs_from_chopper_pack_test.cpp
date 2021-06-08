#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <chopper/build/create_ibfs_from_chopper_pack.hpp>

#include "../api_test.hpp"

using seqan3::operator""_dna4;

struct create_ibfs_from_chopper_pack_test : public ::testing::Test
{
    auto & count_kmers(typename seqan3::interleaved_bloom_filter<>::counting_agent_type<size_t> & ibf_agent,
                       seqan3::dna4_vector const & query,
                       build_config const & config)
    {
        return ibf_agent.bulk_count(query | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));
    }

    auto compare_counts(typename seqan3::interleaved_bloom_filter<> & ibf,
                        seqan3::dna4_vector const & query,
                        std::vector<std::vector<size_t>> const & bins_with_expected_counts,
                        build_config const & config)
    {
        auto agent = ibf.counting_agent<size_t>();
        auto const & counts = this->count_kmers(agent, query, config);

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
            EXPECT_EQ(sum, expected) << "failed for bin: " << bins[0] << std::endl;
        }

        for (size_t i = 0; i < counts.size(); ++i)
            if (without_expected_counts[i])
                EXPECT_EQ(counts[i], 0) << " failed for i=" << i << std::endl;
    }
};

TEST_F(create_ibfs_from_chopper_pack_test, small_example_2_levels)
{
    std::string seq1_filename = data("seq1.fa");
    std::string seq2_filename = data("seq2.fa");
    std::string seq3_filename = data("seq3.fa");

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

    build_data<chopper_pack_record> data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.hibf.size(), 2);

    EXPECT_EQ(data.hibf_bin_levels.size(), 2);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[0], (std::vector<int64_t>{0,0,0,0,0,0,1}));
    EXPECT_RANGE_EQ(data.hibf_bin_levels[1], (std::vector<int64_t>{1,1,1,1,1,1}));

    EXPECT_RANGE_EQ(data.user_bins[0],
                   (std::vector<std::string>
                   {
                        seq1_filename + ";" + seq2_filename,
                        seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        ""
                   }));

    EXPECT_RANGE_EQ(data.user_bins[1],
                   (std::vector<std::string>
                   {
                        seq1_filename,
                        seq2_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq3_filename
                   }));

    auto & high_level_ibf = data.hibf[0];
    auto & low_level_ibf = data.hibf[1];

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

TEST_F(create_ibfs_from_chopper_pack_test, uniform_splitting)
{
    std::string seq1_filename = data("seq1.fa");
    std::string seq2_filename = data("seq2.fa");
    std::string seq3_filename = data("seq3.fa");
    std::string all_seq_filename = data("small.fa");

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
    // Bin 3: split into 3 tb's containing seq1, seq2, seq3
    // Bin 4:  -> belongs to bin 3
    // Bin 5:  -> belongs to bin 3
    // Bin 6: seq1, seq2, seq3

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1
    // Bin 1: seq2
    // Bin 2: split into 3 tb's containing seq1, seq2, seq3
    // Bin 3:  -> belongs to bin 2
    // Bin 4:  -> belongs to bin 2
    // Bin 5: seq3

    build_config config{};
    config.k = 15;
    config.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data<chopper_pack_record> data{};

    create_ibfs_from_chopper_pack(data, config);

    ASSERT_EQ(data.hibf.size(), 2);

    auto & high_level_ibf = data.hibf[0];
    auto & low_level_ibf = data.hibf[1];


    {
        std::unordered_set<size_t> unique_kmers{};

        for (auto & [seq] : sequence_file_t{all_seq_filename})
            for (auto && hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                unique_kmers.insert(hash);

        auto agent = high_level_ibf.counting_agent<size_t>();
        auto const & all_counts = agent.bulk_count(unique_kmers);

        size_t const kmers_per_split_bin = (unique_kmers.size() / 3) + 1;
        EXPECT_EQ(all_counts[3], kmers_per_split_bin);
        EXPECT_EQ(all_counts[4], kmers_per_split_bin);
        EXPECT_EQ(all_counts[5], kmers_per_split_bin - 2);
    }

    {
        std::unordered_set<size_t> unique_kmers{};

        for (auto & [seq] : sequence_file_t{all_seq_filename})
            for (auto && hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                unique_kmers.insert(hash);

        auto agent = low_level_ibf.counting_agent<size_t>();
        auto const & all_counts = agent.bulk_count(unique_kmers);

        size_t const kmers_per_split_bin = (unique_kmers.size() / 3) + 1;
        EXPECT_EQ(all_counts[2], kmers_per_split_bin);
        EXPECT_EQ(all_counts[3], kmers_per_split_bin);
        EXPECT_EQ(all_counts[4], kmers_per_split_bin - 2);
    }
}

TEST_F(create_ibfs_from_chopper_pack_test, same_example_two_levels_but_split_bin_as_hibf_max_bin)
{
    std::string seq1_filename = data("seq1.fa");
    std::string seq2_filename = data("seq2.fa");
    std::string seq3_filename = data("seq3.fa");

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

    build_data<chopper_pack_record> data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.hibf.size(), 2);

    EXPECT_EQ(data.hibf_bin_levels.size(), 2);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[0], (std::vector<int64_t>{0,0,0,0,0,0,1}));
    EXPECT_RANGE_EQ(data.hibf_bin_levels[1], (std::vector<int64_t>{1,1,1,1,1,1}));

    EXPECT_RANGE_EQ(data.user_bins[0],
                   (std::vector<std::string>
                   {
                        seq1_filename + ";" + seq2_filename,
                        seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        ""
                   }));

    EXPECT_RANGE_EQ(data.user_bins[1],
                   (std::vector<std::string>
                   {
                        seq1_filename,
                        seq2_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq3_filename
                   }));

    auto & high_level_ibf = data.hibf[0];
    auto & low_level_ibf = data.hibf[1];

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
    std::string seq1_filename = data("seq1.fa");
    std::string seq2_filename = data("seq2.fa");
    std::string seq3_filename = data("seq3.fa");
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

    build_data<chopper_pack_record> data{};

    create_ibfs_from_chopper_pack(data, config);

    EXPECT_EQ(data.hibf.size(), 6);
    EXPECT_EQ(data.hibf_bin_levels.size(), 6);

    EXPECT_EQ(data.hibf.size(), 6);
    // data.hibf[0] is the high-level IBF
    // data.hibf[1] is the LOW LEVEL IBF 0
    // data.hibf[2] is the LOW LEVEL IBF 0;0
    // data.hibf[3] is the LOW LEVEL IBF 0;0;0
    // data.hibf[4] is the LOW LEVEL IBF 0;1
    // data.hibf[5] is the LOW LEVEL IBF 1

    /* HIGH LEVEL IBF
     * --------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2 (merged)
     * Bin 2: seq2, seq3
     * Bin 3: split but together: seq1, seq2, seq3
     * Bin 4: --> belongs to bin 3
     */
    EXPECT_EQ(data.hibf[0].bin_count(), 5u);
    EXPECT_EQ(data.hibf[0].bin_size(), 114226);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[0], (std::vector<int64_t>{1,5,0,0,0}));
    EXPECT_RANGE_EQ(data.user_bins[0],
                   (std::vector<std::string>
                   {
                        "",
                        "",
                        seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename
                   }));

    /* LOW LEVEL IBF 0
     * ---------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2, seq3 (merged)
     * Bin 2: seq1, seq2
     * Bin 3: seq1, seq3
     * Bin 4: seq2, seq3
     */
    EXPECT_EQ(data.hibf[1].bin_count(), 5u);
    EXPECT_EQ(data.hibf[1].bin_size(), 114226);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[1], (std::vector<int64_t>{2,4,1,1,1}));
    EXPECT_RANGE_EQ(data.user_bins[1],
                   (std::vector<std::string>
                   {
                        "",
                        "",
                        seq1_filename + ";" + seq2_filename,
                        seq1_filename + ";" + seq3_filename,
                        seq2_filename + ";" + seq3_filename
                   }));

    /* LOW LEVEL IBF 0;0
     * -----------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1
     * Bin 2: seq2
     * Bin 3: seq3
     * Bin 4: seq1, seq2
     */
    EXPECT_EQ(data.hibf[2].bin_count(), 5u);
    EXPECT_EQ(data.hibf[2].bin_size(), 95321);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[2], (std::vector<int64_t>{3,2,2,2,2}));
    EXPECT_RANGE_EQ(data.user_bins[2],
                   (std::vector<std::string>
                   {
                        "",
                        seq1_filename,
                        seq2_filename,
                        seq3_filename,
                        seq1_filename + ";" + seq2_filename
                   }));

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
    EXPECT_EQ(data.hibf[3].bin_count(), 13u);
    EXPECT_EQ(data.hibf[3].bin_size(), 92535);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[3], (std::vector<int64_t>(13, 3)));
    EXPECT_RANGE_EQ(data.user_bins[3],
                   (std::vector<std::string>
                   {
                        seq1_filename,
                        seq1_filename,
                        seq1_filename,
                        seq2_filename,
                        seq3_filename,
                        seq3_filename,
                        seq2_filename + ";" + seq3_filename,
                        seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename,
                        seq1_filename + ";" + seq2_filename + ";" + seq3_filename
                   }));

    /* LOW LEVEL IBF 0;1
     * -----------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: seq3
     */
    EXPECT_EQ(data.hibf[4].bin_count(), 3u);
    EXPECT_EQ(data.hibf[4].bin_size(), 92734);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[4], (std::vector<int64_t>{4,4,4}));
    EXPECT_RANGE_EQ(data.user_bins[4],
                   (std::vector<std::string>
                   {
                        seq1_filename,
                        seq2_filename,
                        seq3_filename
                   }));

    /* LOW LEVEL IBF 1
     * ---------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: split but together: seq1, seq2
     * Bin 3: --> belongs to bin 2
     */
    EXPECT_EQ(data.hibf[5].bin_count(), 4u);
    EXPECT_EQ(data.hibf[5].bin_size(), 47561);
    EXPECT_RANGE_EQ(data.hibf_bin_levels[5], (std::vector<int64_t>{5,5,5,5}));
    EXPECT_RANGE_EQ(data.user_bins[5],
                   (std::vector<std::string>
                   {
                        seq1_filename,
                        seq2_filename,
                        seq1_filename + ";" + seq2_filename,
                        seq1_filename + ";" + seq2_filename
                   }));

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    // UNSPECIFIC - unspecific region should be found in all bins that include a whole sequence
    compare_counts(data.hibf[0], unspecific, {{0}, {1}, {2}, {3,4}}, config);
    compare_counts(data.hibf[1], unspecific, {{0}, {1}, {2}, {3}, {4}}, config);
    compare_counts(data.hibf[2], unspecific, {{0}, {1}, {2}, {3}, {4}}, config);
    compare_counts(data.hibf[3], unspecific, {{0, 1, 2}, {3}, {4,5}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.hibf[4], unspecific, {{0}, {1}, {2}}, config);
    compare_counts(data.hibf[5], unspecific, {{0}, {1}, {2,3}}, config);

    // SEQ2 SPECIFIC
    compare_counts(data.hibf[0], seq2_specific, {{0}, {1}, {2}, {3,4}}, config);
    compare_counts(data.hibf[1], seq2_specific, {{0}, {1}, {2}, {4}}, config);
    compare_counts(data.hibf[2], seq2_specific, {{0}, {2}, {4}}, config);
    compare_counts(data.hibf[3], seq2_specific, {{3}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.hibf[4], seq2_specific, {{1}}, config);
    compare_counts(data.hibf[5], seq2_specific, {{1}, {2,3}}, config);

    // SEQ3 SPECIFIC
    compare_counts(data.hibf[0], seq3_specific, {{0}, {2}, {3,4}}, config);
    compare_counts(data.hibf[1], seq3_specific, {{0}, {1}, {3}, {4}}, config);
    compare_counts(data.hibf[2], seq3_specific, {{0}, {3}}, config);
    compare_counts(data.hibf[3], seq3_specific, {{4,5}, {6,7}, {8,9,10,11,12}}, config);
    compare_counts(data.hibf[4], seq3_specific, {{2}}, config);
    compare_counts(data.hibf[5], seq3_specific, {}, config);
}
