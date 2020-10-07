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

    seqan3::test::tmp_filename traversal_dir{"/"};
    std::string traversal_split_bin0{traversal_dir.get_path().string() + "SPLIT_BIN_0.out"};
    std::string traversal_merged_bin2{traversal_dir.get_path().string() + "COLORFUL_MERGED_BIN_2_1.out"};
    std::string traversal_split_bin3{traversal_dir.get_path().string() + "SPLIT_BIN_3.out"};

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

    seqan3::test::tmp_filename output_filename{"/"};
    const char * argv[] = {"./chopper-build", "-k", "15",
                           "-f", data_filename.get_path().c_str(),
                           "-p", traversal_dir.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 9;
    seqan3::argument_parser build_parser{"chopper-build", argc, argv, false};

    chopper_build(build_parser);
}
