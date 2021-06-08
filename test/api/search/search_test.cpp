#include <gtest/gtest.h>

#include <fstream>

#include <chopper/build/create_ibfs_from_chopper_pack.hpp>
#include <chopper/search/chopper_search.hpp>
#include <chopper/search/pair_hash.hpp>
#include <chopper/search/search.hpp>
#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

#include "../api_test.hpp"

using seqan3::operator""_dna4;

struct chopper_search_test : public ::testing::Test
{
    bool is_unique_range(std::vector<std::pair<int32_t, uint32_t>> & rng)
    {
        std::ranges::sort(rng);
        auto it = std::unique(rng.begin(), rng.end());
        return it == rng.end();
    }

    auto compare_result(std::vector<std::pair<int32_t, uint32_t>> & result_set,
                        std::vector<std::pair<int32_t, uint32_t>> && expected_set)
    {
        ASSERT_TRUE(is_unique_range(result_set)); // sanity check actual result
        ASSERT_TRUE(is_unique_range(expected_set)); // sanity check expected result

        EXPECT_RANGE_EQ(result_set, expected_set);
    }
};

TEST_F(chopper_search_test, first_example)
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
    // Bin 3: split bin with: seq1, seq2, seq3
    // Bin 4:   -> belongs to bin 3
    // Bin 5:   -> belongs to bin 3
    // Bin 6: seq1, seq2, seq3 (merged)

    // LOW LEVEL IBF (only one)
    // --------------
    // Bin 0: seq1
    // Bin 1: seq2
    // Bin 2: split bin with: seq1, seq2, seq3
    // Bin 3:   -> belongs to bin 2
    // Bin 4:   -> belongs to bin 2
    // Bin 5: seq3

    build_config bconfig{};
    bconfig.k = 15;
    bconfig.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data<chopper_pack_record> bdata{};

    create_ibfs_from_chopper_pack(bdata, bconfig);

    // move build data to search data
    search_data data{};

    data.hibf = std::move(bdata.hibf);
    data.hibf_bin_levels = std::move(bdata.hibf_bin_levels);
    data.user_bins = std::move(bdata.user_bins);

    search_config config{};
    config.k = bconfig.k;

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    std::vector<size_t> kmers{};

    { // unspecific
        clear_and_compute_kmers(kmers, unspecific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {{0,0},{0,1},{0,2},{0,5},{1,0},{1,1},{1,4},{1,5}});
    }

    { // seq2_specific
        clear_and_compute_kmers(kmers, seq2_specific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {{0,0},{0,2},{0,5},{1,1},{1,4}});
    }

    { // seq3_specific
        clear_and_compute_kmers(kmers, seq3_specific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {{0,1},{0,2},{0,5},{1,4},{1,5}});
    }
}

TEST_F(chopper_search_test, multi_level_example)
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

    build_config bconfig{};
    bconfig.k = 15;
    bconfig.chopper_pack_filename = chopper_pack_filename.get_path().string();

    build_data<chopper_pack_record> bdata{};

    create_ibfs_from_chopper_pack(bdata, bconfig);

    // move build data to search data
    search_data data{};

    data.hibf = std::move(bdata.hibf);
    data.hibf_bin_levels = std::move(bdata.hibf_bin_levels);
    data.user_bins = std::move(bdata.user_bins);

    search_config config{};
    config.k = bconfig.k;

    /* HIGH LEVEL IBF
     * --------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2 (merged)
     * Bin 2: seq2, seq3
     * Bin 3: split but together: seq1, seq2, seq3
     * Bin 4: --> belongs to bin 3
     */

    /* LOW LEVEL IBF 0
     * ---------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2, seq3 (merged)
     * Bin 2: seq1, seq2
     * Bin 3: seq1, seq3
     * Bin 4: seq2, seq3
     */

    /* LOW LEVEL IBF 0;0
     * -----------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1
     * Bin 2: seq2
     * Bin 3: seq3
     * Bin 4: seq1, seq2
     */

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

    /* LOW LEVEL IBF 0;1
     * -----------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: seq3
     */

    /* LOW LEVEL IBF 1
     * ---------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: split but together: seq1, seq2
     * Bin 3: --> belongs to bin 2
     */

    auto unspecific = "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"_dna4;
    auto seq2_specific = "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"_dna4;
    auto seq3_specific = "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"_dna4;

    std::vector<size_t> kmers{};

    { // unspecific
        clear_and_compute_kmers(kmers, unspecific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {
            /*high-level IBF      */ {0,2},{0,4},
            /*LOW LEVEL IBF 0     */ {1,2},{1,3},{1,4},
            /*LOW LEVEL IBF 0;0   */ {2,1},{2,2},{2,3},{2,4},
            /*LOW LEVEL IBF 0;0;0 */ {3,2},{3,3},{3,5},{3,7},{3,12},
            /*LOW LEVEL IBF 0;1   */ {4,0},{4,1},{4,2},
            /*LOW LEVEL IBF 1     */ {5,0},{5,1},{5,3}
        });
    }

    { // seq2_specific
        clear_and_compute_kmers(kmers, seq2_specific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {
            /*high-level IBF      */ {0,2},{0,4},
            /*LOW LEVEL IBF 0     */ {1,2},{1,4},
            /*LOW LEVEL IBF 0;0   */ {2,2},{2,4},
            /*LOW LEVEL IBF 0;0;0 */ {3,3},{3,7},{3,12},
            /*LOW LEVEL IBF 0;1   */ {4,1},
            /*LOW LEVEL IBF 1     */ {5,1},{5,3}
        });
    }

    { // seq3_specific
        clear_and_compute_kmers(kmers, seq3_specific, config);
        std::vector<std::pair<int32_t, uint32_t>> result{};

        search(result, kmers, data, config, 0); // start at top level ibf

        this->compare_result(result, {
            /*high-level IBF      */ {0,2},{0,4},
            /*LOW LEVEL IBF 0     */ {1,3},{1,4},
            /*LOW LEVEL IBF 0;0   */ {2,3},
            /*LOW LEVEL IBF 0;0;0 */ {3,5},{3,7},{3,12},
            /*LOW LEVEL IBF 0;1   */ {4,2}
            /*LOW LEVEL IBF 1     */
        });
    }
}
