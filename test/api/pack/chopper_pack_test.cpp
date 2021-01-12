#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/pack/chopper_pack.hpp>

TEST(chopper_pack_test, small_example)
{
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    {
        std::ofstream fout{data_filename.get_path()};
        fout << "seq0\t500\n"
             << "seq1\t1000\n"
             << "seq2\t500\n"
             << "seq3\t500\n"
             << "seq4\t500\n"
             << "seq5\t500\n"
             << "seq6\t500\n"
             << "seq7\t500\n";
    }

    seqan3::test::tmp_filename output_filename{"output.binning"};
    const char * argv[] = {"./chopper-pack", "-b", "4",
                           "-f", data_filename.get_path().c_str(), "-o", output_filename.get_path().c_str()};
    int argc = 7;
    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};

    chopper_pack(pack_parser);

    std::string expected_file
    {
        "#MERGED_BIN_2 max_bin_id:16\n"
        "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2\n"
        "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
        "SPLIT_BIN_0\tseq7\t1\t500\n"
        "SPLIT_BIN_1\tseq6\t1\t500\n"
        "MERGED_BIN_2_0\tseq0\t16\t32\n"
        "MERGED_BIN_2_16\tseq2\t12\t42\n"
        "MERGED_BIN_2_28\tseq3\t12\t42\n"
        "MERGED_BIN_2_40\tseq4\t12\t42\n"
        "MERGED_BIN_2_52\tseq5\t12\t42\n"
        "SPLIT_BIN_3\tseq1\t1\t1000\n"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file);
}
