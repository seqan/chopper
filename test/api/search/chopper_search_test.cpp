#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include <chopper/build/chopper_build.hpp>
#include <chopper/build/create_ibfs_from_chopper_pack.hpp>
#include <chopper/search/pair_hash.hpp>
#include <chopper/search/chopper_search.hpp>
#include <chopper/search/search.hpp>
#include <chopper/search/search_data.hpp>
#include <chopper/search/sync_out.hpp>

#include "../api_test.hpp"

using seqan3::operator""_dna4;

struct chopper_search_test : public ::testing::Test
{
    bool is_unique_range(std::vector<std::pair<int32_t, uint32_t>> v)
    {
        std::sort(v.begin(), v.end());
        auto it = std::unique(v.begin(), v.end());
        return it == v.end();
    }

    auto compare_result(std::unordered_set<std::pair<int32_t, uint32_t>, pair_hash> const & result_set,
                        std::vector<std::pair<int32_t, uint32_t>> const & bins_in_result_set)
    {
        ASSERT_TRUE(is_unique_range(bins_in_result_set)); // sanity check expected result

        EXPECT_EQ(result_set.size(), bins_in_result_set.size());

        for (auto const & bin : bins_in_result_set)
            EXPECT_TRUE(result_set.find(bin) != result_set.end());
    }

    static inline std::string const string_from_file(std::filesystem::path const & path,
                                                     std::ios_base::openmode const mode = std::ios_base::in)
    {
        std::ifstream file_stream(path, mode);
        if (!file_stream.is_open())
            throw std::logic_error{"[string_from_file] Cannot open " + path.string()};
        std::stringstream file_buffer;
        file_buffer << file_stream.rdbuf();
        return {file_buffer.str()};
    }
};

TEST_F(chopper_search_test, write_result)
{
    std::vector<std::pair<int32_t, uint32_t>> result;
    result.emplace_back(0, 0u);
    result.emplace_back(0, 1u);
    result.emplace_back(0, 2u);
    result.emplace_back(0, 5u);
    result.emplace_back(1, 0u);
    result.emplace_back(1, 2u);

    std::string query_id{"query1"};

    search_data data;

    data.user_bins.resize_filename(9);
    data.user_bins.filename_at(0) = "user_bin_0-0";
    data.user_bins.filename_at(1) = "user_bin_0-1";
    data.user_bins.filename_at(2) = "user_bin_0-2";
    data.user_bins.filename_at(3) = "user_bin_0-3";
    data.user_bins.filename_at(4) = "user_bin_0-4";
    data.user_bins.filename_at(5) = "user_bin_0-5";
    data.user_bins.filename_at(6) = "user_bin_1-0";
    data.user_bins.filename_at(7) = "user_bin_1-1";
    data.user_bins.filename_at(8) = "user_bin_1-2";

    data.user_bins.resize_bins(2);
    data.user_bins.bin_at(0) = std::vector<int64_t>{0, 1, 2, 3, 4, 5};
    data.user_bins.bin_at(1) = std::vector<int64_t>{6, 7, 8};

    seqan3::test::tmp_filename tmp_file{"chopper_search_test_write_result"};
    {
        sync_out out_file{tmp_file.get_path()};
        write_header(data, out_file);
        std::string buffer{};
        write_result(buffer, result, query_id, data, out_file);
    }

    std::string expected
    {
        "#0\tuser_bin_0-0\n"
        "#1\tuser_bin_0-1\n"
        "#2\tuser_bin_0-2\n"
        "#3\tuser_bin_0-3\n"
        "#4\tuser_bin_0-4\n"
        "#5\tuser_bin_0-5\n"
        "#6\tuser_bin_1-0\n"
        "#7\tuser_bin_1-1\n"
        "#8\tuser_bin_1-2\n"
        "#QUERY_NAME\tUSER_BINS\n"
        "query1\t0,1,2,5,6,8\n"
    };

    EXPECT_EQ(string_from_file(tmp_file.get_path()), expected);
}

TEST_F(chopper_search_test, first_example)
{
    std::string seq1_filename = data("seq1.fa");
    std::string seq2_filename = data("seq2.fa");
    std::string seq3_filename = data("seq3.fa");

    seqan3::test::tmp_filename chopper_pack_filename{"small.pack"};

    // generate data files (chopper pack)
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:6\n"
             << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << seq1_filename << ";" << seq2_filename                         << "\t0\t1\t500\n"
             << seq3_filename                                                 << "\t1\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t2\t1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t3\t3\t500\n"
             << seq1_filename                                                 << "\t6;0\t1;1\t500\n"
             << seq2_filename                                                 << "\t6;1\t1;1\t500\n"
             << seq1_filename << ";" << seq2_filename << ";" << seq3_filename << "\t6;2\t1;3\t500\n"
             << seq3_filename                                                 << "\t6;5\t1;1\t500\n";
    }

    seqan3::test::tmp_filename output_path{"chopper.test.index"};

    // first execute building to generate the needed input data
    {
        const char * argv[] = {"./chopper-build",
                               "-k", "15",
                               "-p", chopper_pack_filename.get_path().c_str(),
                               "-o", output_path.get_path().c_str()};
        int argc = 7;
        seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

        testing::internal::CaptureStderr();
        auto parse_res = chopper_build(build_parser);
        std::string std_cerr = testing::internal::GetCapturedStderr();
        ASSERT_EQ(parse_res, 0) << std_cerr;

        ASSERT_TRUE(std::filesystem::exists(output_path.get_path()));
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

    // generate query file
    seqan3::test::tmp_filename query_filename{"queries.fa"};
    {
        std::ofstream fout{query_filename.get_path()};
        fout << ">unspecific\n"
             << "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG\n"
             << ">seq2_specific\n"
             << "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG\n"
             << ">seq3_specific\n"
             << "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT\n";
    }

    seqan3::test::tmp_filename output_filename{"search.out"};
    std::filesystem::path timing_filename{output_filename.get_path()};
    timing_filename += ".time";

    const char * argv[] = {"./chopper-search",
                           "-i", output_path.get_path().c_str(),
                           "-q", query_filename.get_path().c_str(),
                           "-o", output_filename.get_path().c_str(),
                           "--time"};
    int argc = 8;
    seqan3::argument_parser search_parser{"chopper-search", argc, argv, seqan3::update_notifications::off};

    testing::internal::CaptureStderr();
    auto parse_res = chopper_search(search_parser);
    std::string std_cerr = testing::internal::GetCapturedStderr();
    ASSERT_EQ(parse_res, 0) << std_cerr;

    std::string expected
    {
        "#0\t" + seq1_filename + "\n"
        "#1\t" + seq2_filename + "\n"
        "#2\t" + seq1_filename + ";" + seq2_filename + ";" + seq3_filename + "\n"
        "#3\t" + seq3_filename + "\n"
        "#4\t" + seq1_filename + ";" + seq2_filename + "\n"
        "#5\t" + seq3_filename + "\n"
        "#6\t" + seq1_filename + ";" + seq2_filename + ";" + seq3_filename + "\n"
        "#7\t" + seq1_filename + ";" + seq2_filename + ";" + seq3_filename + "\n"
        "#QUERY_NAME\tUSER_BINS\n"
        "unspecific\t0,1,2,3,4,5,6,7\n"
        "seq2_specific\t1,2,4,6,7\n"
        "seq3_specific\t2,3,5,6,7\n"
    };

    EXPECT_TRUE(std::filesystem::exists(timing_filename));
    std::filesystem::remove(timing_filename);
    EXPECT_EQ(string_from_file(output_filename.get_path()), expected);
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

    seqan3::test::tmp_filename output_path{"chopper.test.index"};
    // first execute building to generate the needed input data
    {
        const char * argv[] = {"./chopper-build",
                               "-k", "15",
                               "-p", chopper_pack_filename.get_path().c_str(),
                               "-o", output_path.get_path().c_str()};
        int argc = 7;
        seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

        testing::internal::CaptureStderr();
        auto parse_res = chopper_build(build_parser);
        std::string std_cerr = testing::internal::GetCapturedStderr();
        ASSERT_EQ(parse_res, 0) << std_cerr;

        ASSERT_TRUE(std::filesystem::exists(output_path.get_path()));
    }

    /* HIGH LEVEL IBF
     * --------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2 (merged)
     * Bin 2: seq2, seq3
     * Bin 3: split but together: seq1, seq2, seq3
     * Bin 4: --> belongs to bin 3
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

    /* LOW LEVEL IBF 0;0
     * -----------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1
     * Bin 2: seq2
     * Bin 3: seq3
     * Bin 4: seq1, seq2
     */

    /* LOW LEVEL IBF 0;1
     * -----------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: seq3
     */

    /* LOW LEVEL IBF 0
     * ---------------
     * Bin 0: seq1, seq2, seq3 (merged)
     * Bin 1: seq1, seq2, seq3 (merged)
     * Bin 2: seq1, seq2
     * Bin 3: seq1, seq3
     * Bin 4: seq2, seq3
     */

    /* LOW LEVEL IBF 1
     * ---------------
     * Bin 0: seq1
     * Bin 1: seq2
     * Bin 2: split but together: seq1, seq2
     * Bin 3: --> belongs to bin 2
     */

    // generate query file
    seqan3::test::tmp_filename query_filename{"queries.fa"};
    {
        std::ofstream fout{query_filename.get_path()};
        fout << ">unspecific\n"
             << "ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG\n"
             << ">seq2_specific\n"
             << "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG\n"
             << ">seq3_specific\n"
             << "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT\n";
    }

    seqan3::test::tmp_filename output_filename{"search.out"};

    const char * argv[] = {"./chopper-search",
                           "-i", output_path.get_path().c_str(),
                           "-q", query_filename.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 7;
    seqan3::argument_parser search_parser{"chopper-search", argc, argv, seqan3::update_notifications::off};

    testing::internal::CaptureStderr();
    auto parse_res = chopper_search(search_parser);
    std::string std_cerr = testing::internal::GetCapturedStderr();
    ASSERT_EQ(parse_res, 0) << std_cerr;

    std::string expected
    {
        "#0\t" + seq12_filename + "\n"
        "#1\t" + seq2_filename + "\n"
        "#2\t" + seq1_filename + "\n"
        "#3\t" + seq3_filename + "\n"
        "#4\t" + seq23_filename + "\n"
        "#5\t" + seq123_filename + "\n"
        "#6\t" + seq1_filename + "\n"
        "#7\t" + seq2_filename + "\n"
        "#8\t" + seq3_filename + "\n"
        "#9\t" + seq3_filename + "\n"
        "#10\t" + seq1_filename + "\n"
        "#11\t" + seq2_filename + "\n"
        "#12\t" + seq12_filename + "\n"
        "#13\t" + seq13_filename + "\n"
        "#14\t" + seq23_filename + "\n"
        "#15\t" + seq12_filename + "\n"
        "#16\t" + seq1_filename + "\n"
        "#17\t" + seq2_filename + "\n"
        "#18\t" + seq23_filename + "\n"
        "#19\t" + seq123_filename + "\n"
        "#QUERY_NAME\tUSER_BINS\n"
        "unspecific\t0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19\n"
        "seq2_specific\t0,1,4,5,7,11,12,14,15,17,18,19\n"
        "seq3_specific\t3,4,5,8,9,13,14,18,19\n"
    };

    EXPECT_EQ(string_from_file(output_filename.get_path()), expected);
}
