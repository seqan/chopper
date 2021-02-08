#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/build/create_ibfs_from_chopper_pack.hpp>
#include <chopper/search/chopper_search.hpp>

using seqan3::operator""_dna4;

TEST(chopper_search_test, first_example)
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
    // Bin 6: seq1, seq2, seq3 (merged)

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1
    // Bin 1: seq2
    // Bin 2: ? not easily to determine since kmers are split independently
    // Bin 3: ? not easily to determine since kmers are split independently
    // Bin 4: ? not easily to determine since kmers are split independently
    // Bin 5: seq3

    build_config bconfig{};
    bconfig.k = 15;
    bconfig.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data bdata{};

    create_ibfs_from_chopper_pack(bdata, bconfig);

    // move build data to search data
    search_data data{};
    for (auto & ibf : bdata.hibf)
        data.hibf.push_back(static_cast<seqan3::technical_binning_directory<> &&>(ibf));

    data.hibf_bin_levels = std::move(bdata.hibf_bin_levels);
    data.user_bins = std::move(bdata.user_bins);

    search_config config{};
    config.k = bconfig.k;

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    // auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    // auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    std::vector<size_t> unspecific_kmers = compute_kmers(unspecific, config);
    std::unordered_set<std::pair<int32_t, uint32_t>, pair_hash> unspecific_result{};

    search(unspecific_result, unspecific_kmers, data, config, 0); // start at top level ibf

    EXPECT_TRUE(unspecific_result.find({0,0}) != unspecific_result.end());
    EXPECT_TRUE(unspecific_result.find({0,1}) != unspecific_result.end());
    EXPECT_TRUE(unspecific_result.find({0,2}) != unspecific_result.end());

    EXPECT_TRUE(unspecific_result.find({1,0}) != unspecific_result.end());
    EXPECT_TRUE(unspecific_result.find({1,1}) != unspecific_result.end());
    EXPECT_TRUE(unspecific_result.find({1,5}) != unspecific_result.end());

    seqan3::debug_stream << unspecific_result << std::endl;
}
