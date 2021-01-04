#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <chopper/build/chopper_build.hpp>

using seqan3::operator""_dna4;

TEST(chopper_count_test, small_example_parallel_2_threads)
{
    std::string input_filename1 = DATADIR"small.fa";
    std::string input_filename2 = DATADIR"small2.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    // generate data files
    {
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_6\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*SPLIT_BIN_0*/
             << input_filename1 << "\tseq1\t0\t400\t0\t-\n"
             << input_filename1 << "\tseq2\t0\t480\t0\t-\n"
             << input_filename1 << "\tseq3\t0\t481\t1\t-\n"
             /*SPLIT_BIN_2*/
             << input_filename2 << "\tseq10\t0\t400\t2\t-\n"
             << input_filename2 << "\tseq20\t0\t480\t2\t-\n"
             << input_filename2 << "\tseq30\t0\t481\t2\t-\n"
             /*SPLIT_BIN_3*/
             << input_filename2 << "\tseq10\t0\t400\t3\t-\n"
             << input_filename2 << "\tseq20\t0\t480\t4\t-\n"
             << input_filename2 << "\tseq30\t0\t481\t5\t-\n"
             /*MERGED_BIN_6*/
             << input_filename1 << "\tseq1\t0\t400\t6\t0\n"
             << input_filename1 << "\tseq2\t0\t480\t6\t1\n"
             << input_filename1 << "\tseq3\t0\t481\t6\t2\n"
             << input_filename1 << "\tseq1\t0\t400\t6\t3\n"
             << input_filename1 << "\tseq1\t0\t400\t6\t3\n"
             << input_filename1 << "\tseq2\t0\t480\t6\t3\n"
             << input_filename1 << "\tseq2\t0\t480\t6\t3\n"
             << input_filename1 << "\tseq3\t0\t481\t6\t4\n"
             << input_filename1 << "\tseq3\t0\t481\t6\t4\n";
    }

    seqan3::test::tmp_filename output_prefix{"TEST_"};
    const char * argv[] = {"./chopper-build",
                           "-k", "15",
                           "-p", traversal_filename.get_path().c_str(),
                           "-o", output_prefix.get_path().c_str()};
    int argc = 7;
    seqan3::argument_parser build_parser{"chopper-build", argc, argv, false};

    chopper_build(build_parser);

    EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "high_level.ibf"));
    EXPECT_TRUE(std::filesystem::exists(output_prefix.get_path().string() + "low_level_6.ibf"));
}
