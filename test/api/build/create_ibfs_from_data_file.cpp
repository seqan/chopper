#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <chopper/build/create_ibfs_from_data_file.hpp>

using seqan3::operator""_dna4;

TEST(chopper_count_test, small_example_parallel_2_threads)
{
    std::string input_filename1 = DATADIR"small.fa";
    std::string input_filename2 = DATADIR"small2.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    seqan3::test::tmp_filename traversal_dir{""};
    std::string traversal_split_bin0{traversal_dir.get_path().string() + "/SPLIT_BIN_0.out"};
    std::string traversal_merged_bin2{traversal_dir.get_path().string() + "/COLORFUL_MERGED_BIN_2_1.out"};
    std::string traversal_split_bin3{traversal_dir.get_path().string() + "/SPLIT_BIN_3.out"};

    // generate data files
    {
        std::ofstream fout{data_filename.get_path()};
        fout << "BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
             << "SPLIT_BIN_0\t" << input_filename1 << "\t2\t500\n"
             << "SPLIT_BIN_1\t" << input_filename1 + "\t1\t500\n"
             << "COLORFUL_MERGED_BIN_2_0\t" << input_filename1 << "\t1\t2500\n"
             << "COLORFUL_MERGED_BIN_2_1\t" << input_filename1 << ";" << input_filename2 << "\t2\t2500\n"
             << "SPLIT_BIN_3\t" << input_filename2 + "\t3\t1000\n";
    }

    {
        std::ofstream fout{traversal_split_bin0};
        fout << "FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER\n"
             << input_filename1 << "\tseq1\t0\t400\t0\n"
             << input_filename1 << "\tseq2\t0\t480\t0\n"
             << input_filename1 << "\tseq3\t0\t481\t1\n";
    }
    {
        std::ofstream fout{traversal_merged_bin2};
        fout << "FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER\n"
             << input_filename1 << "\tseq1\t0\t400\t0\n"
             << input_filename2 << "\tseq10\t0\t400\t0\n"
             << input_filename1 << "\tseq2\t0\t480\t0\n"
             << input_filename2 << "\tseq20\t0\t480\t0\n"
             << input_filename1 << "\tseq3\t0\t481\t1\n"
             << input_filename2 << "\tseq30\t0\t481\t1\n";
    }
    {
        std::ofstream fout{traversal_split_bin3};
        fout << "FILE_ID\tSEQ_ID\tBEGIN\tEND\tBIN_NUMBER\n"
             << input_filename2 << "\tseq10\t0\t400\t0\n"
             << input_filename2 << "\tseq20\t0\t480\t1\n"
             << input_filename2 << "\tseq30\t0\t481\t2\n";
    }

    build_config config{};
    config.k = 15;
    config.traversal_path_prefix = traversal_dir.get_path().string() + "/";
    config.binning_filename = data_filename.get_path().string();

    auto && [high_level_ibf, low_level_ibfs] = create_ibfs_from_data_file(config);

    EXPECT_EQ(config.high_level_ibf_num_technical_bins, 7);
    EXPECT_EQ(low_level_ibfs.size(), 1u);
    EXPECT_EQ(low_level_ibfs[0].bin_count(), 3u);

    // The traversal files are made up like the following
    // HIGH LEVEL IBF
    // --------------
    // Bin 0: seq1, seq2
    // Bin 1: seq3
    // Bin 2: seq1, seq2, seq3
    // Bin 3: seq1, seq2, seq3
    // Bin 4: seq1
    // Bin 5: seq2
    // Bin 6: seq3

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1, seq2, seq3
    // Bin 1: seq1, seq2
    // Bin 2: seq3

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    auto high_level_agent = high_level_ibf.membership_agent();
    auto low_level_agent = low_level_ibfs[0].membership_agent();

    {
        std::vector<size_t> high_level_counts(high_level_agent.result_buffer.size());
        std::vector<size_t> low_level_counts(low_level_agent.result_buffer.size());

        for (auto hash : unspecific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        {
            auto const & high_res_v = high_level_agent.bulk_contains(hash);
            auto const & low_res_v = low_level_agent.bulk_contains(hash);

            for (size_t i = 0; i < high_res_v.size(); ++i)
                high_level_counts[i] += high_res_v[i];

            for (size_t i = 0; i < low_res_v.size(); ++i)
                low_level_counts[i] += low_res_v[i];
        }

        size_t expected = std::ranges::distance(unspecific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        for (size_t i = 0; i < high_level_counts.size(); ++i)
            EXPECT_EQ(high_level_counts[i], expected) << "[HIGH LEVEL] Failed for bin " << i << std::endl;

        for (size_t i = 0; i < low_level_counts.size(); ++i)
            EXPECT_EQ(low_level_counts[i], expected) << "[LOW LEVEL] Failed for bin " << i << std::endl;
    }

    {
        std::vector<size_t> high_level_counts(high_level_agent.result_buffer.size());
        std::vector<size_t> low_level_counts(low_level_agent.result_buffer.size());

        for (auto hash : seq2_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        {
            auto const & high_res_v = high_level_agent.bulk_contains(hash);
            auto const & low_res_v = low_level_agent.bulk_contains(hash);

            for (size_t i = 0; i < high_res_v.size(); ++i)
                high_level_counts[i] += high_res_v[i];

            for (size_t i = 0; i < low_res_v.size(); ++i)
                low_level_counts[i] += low_res_v[i];
        }

        size_t expected = std::ranges::distance(seq2_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(high_level_counts[0], expected);
        EXPECT_EQ(high_level_counts[1], 1); // 1 by change
        EXPECT_EQ(high_level_counts[2], expected);
        EXPECT_EQ(high_level_counts[3], expected);
        EXPECT_EQ(high_level_counts[4], 1); // 1 by change
        EXPECT_EQ(high_level_counts[5], expected);
        EXPECT_EQ(high_level_counts[6], 1); // 1 by change

        EXPECT_EQ(low_level_counts[0], expected);
        EXPECT_EQ(low_level_counts[1], expected);
        EXPECT_EQ(low_level_counts[2], 1); // 1 by change
    }

    {
        std::vector<size_t> high_level_counts(high_level_agent.result_buffer.size());
        std::vector<size_t> low_level_counts(low_level_agent.result_buffer.size());

        for (auto hash : seq3_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
        {
            auto const & high_res_v = high_level_agent.bulk_contains(hash);
            auto const & low_res_v = low_level_agent.bulk_contains(hash);

            for (size_t i = 0; i < high_res_v.size(); ++i)
                high_level_counts[i] += high_res_v[i];

            for (size_t i = 0; i < low_res_v.size(); ++i)
                low_level_counts[i] += low_res_v[i];
        }

        size_t expected = std::ranges::distance(seq3_specific | seqan3::views::kmer_hash(seqan3::ungapped{config.k}));

        EXPECT_EQ(high_level_counts[0], 0);
        EXPECT_EQ(high_level_counts[1], expected);
        EXPECT_EQ(high_level_counts[2], expected);
        EXPECT_EQ(high_level_counts[3], expected);
        EXPECT_EQ(high_level_counts[4], 0);
        EXPECT_EQ(high_level_counts[5], 0);
        EXPECT_EQ(high_level_counts[6], expected);

        EXPECT_EQ(low_level_counts[0], expected);
        EXPECT_EQ(low_level_counts[1], 0);
        EXPECT_EQ(low_level_counts[2], expected);
    }
}
