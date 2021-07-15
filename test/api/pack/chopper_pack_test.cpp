#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/chopper_pack.hpp>

#include "../api_test.hpp"

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
        "#HIGH_LEVEL_IBF max_bin_id:3\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#MERGED_BIN_2 max_bin_id:0\n"
        "#MERGED_BIN_3 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
        "seq7\t0\t1\n"
        "seq5\t1;0\t1;7\n"
        "seq6\t1;7\t1;57\n"
        "seq3\t2;0\t1;7\n"
        "seq4\t2;7\t1;57\n"
        "seq1\t3;0\t1;8\n"
        "seq0\t3;8\t1;4\n"
        "seq2\t3;12\t1;52\n"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file);
}

TEST(chopper_pack_test, small_example_debug)
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
    const char * argv[] = {"./chopper-pack", "-b", "4", "--debug",
                           "-f", data_filename.get_path().c_str(), "-o", output_filename.get_path().c_str()};
    int argc = 8;
    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};

    chopper_pack(pack_parser);

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:3\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#MERGED_BIN_2 max_bin_id:0\n"
        "#MERGED_BIN_3 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX\n"
        "seq7\t0\t1\t500\t1000\t1.00\t4\n"
        "seq5\t1;0\t1;7\t1000;72\t1000;124\t1.00;3.15\t4;64\n"
        "seq6\t1;7\t1;57\t1000;9\t1000;124\t1.00;14.24\t4;64\n"
        "seq3\t2;0\t1;7\t1000;72\t1000;124\t1.00;3.15\t4;64\n"
        "seq4\t2;7\t1;57\t1000;9\t1000;124\t1.00;14.24\t4;64\n"
        "seq1\t3;0\t1;8\t2000;125\t1000;126\t1.00;3.43\t4;64\n"
        "seq0\t3;8\t1;4\t2000;125\t1000;126\t1.00;2.23\t4;64\n"
        "seq2\t3;12\t1;52\t2000;10\t1000;126\t1.00;13.21\t4;64\n"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file);
}
