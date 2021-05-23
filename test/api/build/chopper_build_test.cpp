#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <chopper/build/chopper_build.hpp>

#include "../api_test.hpp"

using seqan3::operator""_dna4;

// TEST(chopper_build_test, chopper_split_file)
// {
//     std::string input_filename1 = data("small.fa");
//     std::string input_filename2 = data("small2.fa");
//     seqan3::test::tmp_filename data_filename{"data.tsv"};

//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};

//     // generate data files
//     {
//         std::ofstream fout{chopper_split_filename.get_path()};
//         fout << "#MERGED_BIN_6 max_bin_id:0\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_6\n"
//              << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//              /*SPLIT_BIN_0*/
//              << input_filename1 << "\tseq1\t0\t400\t0\t-\n"
//              << input_filename1 << "\tseq2\t0\t480\t0\t-\n"
//              << input_filename1 << "\tseq3\t0\t481\t1\t-\n"
//              /*SPLIT_BIN_2*/
//              << input_filename2 << "\tseq10\t0\t400\t2\t-\n"
//              << input_filename2 << "\tseq20\t0\t480\t2\t-\n"
//              << input_filename2 << "\tseq30\t0\t481\t2\t-\n"
//              /*SPLIT_BIN_3*/
//              << input_filename2 << "\tseq10\t0\t400\t3\t-\n"
//              << input_filename2 << "\tseq20\t0\t480\t4\t-\n"
//              << input_filename2 << "\tseq30\t0\t481\t5\t-\n"
//              /*MERGED_BIN_6*/
//              << input_filename1 << "\tseq1\t0\t400\t6\t0\n"
//              << input_filename1 << "\tseq2\t0\t480\t6\t1\n"
//              << input_filename1 << "\tseq3\t0\t481\t6\t2\n"
//              << input_filename1 << "\tseq1\t0\t400\t6\t3\n"
//              << input_filename1 << "\tseq1\t0\t400\t6\t3\n"
//              << input_filename1 << "\tseq2\t0\t480\t6\t3\n"
//              << input_filename1 << "\tseq2\t0\t480\t6\t3\n"
//              << input_filename1 << "\tseq3\t0\t481\t6\t4\n"
//              << input_filename1 << "\tseq3\t0\t481\t6\t4\n";
//     }

//     seqan3::test::tmp_filename output_prefix{"TEST_"};
//     const char * argv[] = {"./chopper-build",
//                            "-k", "15",
//                            "-s", chopper_split_filename.get_path().c_str(),
//                            "-o", output_prefix.get_path().c_str()};
//     int argc = 7;
//     seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

//     chopper_build(build_parser);

//     EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "high_level.ibf"));
//     EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "low_level_6.ibf"));
// }

TEST(chopper_build_test, chopper_pack_file)
{
    std::string input_filename1 = data("small.fa");
    std::string input_filename2 = data("small2.fa");

    seqan3::test::tmp_filename chopper_pack_filename{"test.pack"};

    // generate data files
    {
        std::ofstream fout{chopper_pack_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:6\n"
             << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << input_filename1 << "\t0\t1\t500\n"
             << input_filename1 << "\t1\t1\t500\n"
             << input_filename1 << "\t2\t1\t500\n"
             << input_filename1 << "\t3\t3\t500\n"
             << input_filename1 << "\t6;0\t1;1\t500\n"
             << input_filename2 << "\t6;1\t1;1\t500\n"
             << input_filename1 << "\t6;2\t1;3\t500\n"
             << input_filename1 << "\t6;5\t1;1\t500\n";
    }

    seqan3::test::tmp_filename output_path{"chopper.test.index"};
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

    EXPECT_TRUE(std::filesystem::exists(output_path.get_path()));
}

// TEST(chopper_build_test, pack_and_split_file_given)
// {
//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};
//     seqan3::test::tmp_filename chopper_pack_filename{"test.pack"};

//     { // generate data files
//         std::ofstream fout{chopper_pack_filename.get_path()};
//         fout << "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_0\n"
//              << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
//              << "SPLIT_BIN_0\tfile1\t1\t300\n";

//         std::ofstream fout2{chopper_split_filename.get_path()};
//         fout2 << "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_0\n"
//               << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//               << "file1\tseq1\t0\t400\t0\t-\n";
//     }

//     seqan3::test::tmp_filename output_prefix{"TEST_"};
//     const char * argv[] = {"./chopper-build",
//                            "-k", "15",
//                            "-p", chopper_pack_filename.get_path().c_str(),
//                            "-s", chopper_split_filename.get_path().c_str(),
//                            "-o", output_prefix.get_path().c_str()};
//     int argc = 9;
//     seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

//     testing::internal::CaptureStderr();
//     EXPECT_EQ(chopper_build(build_parser), -1);
//     std::string std_cerr = testing::internal::GetCapturedStderr();
//     EXPECT_EQ(std_cerr, "[CHOPPER BIULD ERROR] Options -p/--pack-file and -s/--split_file are mututal exclusive.\n");

//     EXPECT_FALSE(std::filesystem::exists(output_prefix.get_path().string() + "high_level.ibf"));
// }

// TEST(chopper_build_test, neither_pack_nor_split_file_given)
// {
//     seqan3::test::tmp_filename output_prefix{"TEST_"};
//     const char * argv[] = {"./chopper-build",
//                            "-k", "15",
//                            "-o", output_prefix.get_path().c_str()};
//     int argc = 5;
//     seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

//     testing::internal::CaptureStderr();
//     EXPECT_EQ(chopper_build(build_parser), -1);
//     std::string std_cerr = testing::internal::GetCapturedStderr();
//     EXPECT_EQ(std_cerr, "[CHOPPER BIULD ERROR] Either option -p/--pack-file or -s/--split_file must be provided.\n");
// }

// TEST(chopper_build_test, create_output_dir_if_it_does_not_exist)
// {
//     std::string input_filename1 = data("small.fa");
//     std::string input_filename2 = data("small2.fa");

//     seqan3::test::tmp_filename chopper_pack_filename{"test.pack"};

//     // generate data files
//     {
//         std::ofstream fout{chopper_pack_filename.get_path()};
//         fout << "#MERGED_BIN_6 max_bin_id:0\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_6\n"
//              << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
//              << "SPLIT_BIN_0\t" << input_filename1 << "\t2\t300\n"
//              << "SPLIT_BIN_2\t" << input_filename2 << "\t1\t600\n"
//              << "SPLIT_BIN_3\t" << input_filename2 << "\t3\t200\n"
//              << "MERGED_BIN_6_0\t" << input_filename1 << "\t3\t200\n"
//              << "MERGED_BIN_6_3\t" << input_filename1 << "\t2\t300\n";
//     }

//     seqan3::test::tmp_filename output_prefix{"some/deep/directory/TEST_"};
//     const char * argv[] = {"./chopper-build",
//                            "-k", "15",
//                            "-p", chopper_pack_filename.get_path().c_str(),
//                            "-o", output_prefix.get_path().c_str()};
//     int argc = 7;
//     seqan3::argument_parser build_parser{"chopper-build", argc, argv, seqan3::update_notifications::off};

//     chopper_build(build_parser);

//     ASSERT_TRUE(std::filesystem::exists(output_prefix.get_path().parent_path()));
//     EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "high_level.ibf"));
//     EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "low_level_6.ibf"));
// }
