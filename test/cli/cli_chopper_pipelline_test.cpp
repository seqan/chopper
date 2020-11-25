#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/detail_bin_prefixes.hpp>

#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_pipeline)
{
    // CHOPPER COUNT
    // =========================================================================
    std::string seq_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq_filename << '\t' << "TAX1\n"
             << seq_filename << '\t' << "TAX2\n"
             << seq_filename << '\t' << "TAX2\n"
             << seq_filename << '\t' << "TAX3\n";
    }

    cli_test_result count_result = execute_app("chopper", "count",
                                               "-k", "15",
                                               "-w", "25",
                                               "-t", "2",
                                               "-c", "2",
                                               "-f", taxa_filename.get_path().c_str());

    std::string expected
    {
        seq_filename + "\t95\tTAX3\n" +
        seq_filename + ";" + seq_filename + "\t95\tTAX2\n" +
        seq_filename + "\t95\tTAX1\n"
    };

    // compare intermediate results
    EXPECT_EQ(expected, count_result.out);

    // Write count_result to output file (user would pipe the output)
    seqan3::test::tmp_filename const count_filename{"data.tsv"};

    {
        std::ofstream fout{count_filename.get_path()};
        fout << count_result.out;
    }

    // CHOPPER PACK
    // =========================================================================
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    cli_test_result pack_result = execute_app("chopper", "pack",
                                              "-b", "2",
                                              "-f", count_filename.get_path().c_str(),
                                              "-o", binning_filename.get_path().c_str());

    std::string expected_file
    {
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_1\n"
        "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n" +
        std::string{split_bin_prefix} + "_0\t" + seq_filename + "\t1\t95\n" +
        std::string{merged_bin_prefix} + "_1_0\t" + seq_filename + "\t32\t3\n" +
        std::string{merged_bin_prefix} + "_1_1\t" + seq_filename + ";" + seq_filename + "\t32\t3\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        EXPECT_EQ(output_file_str, expected_file);
    }

    // CHOPPER SPLIT
    // =========================================================================
    seqan3::test::tmp_filename const traversal_path_prefix{""};

    cli_test_result split_result = execute_app("chopper", "split",
                                               "-k", "15",
                                               "-w", "25",
                                               "-f", binning_filename.get_path().c_str(),
                                               "-o", traversal_path_prefix.get_path().c_str());

    std::vector<std::string> const traversal_filenames
    {
        "LOW_LEVEL_IBF_1.out"
    };

    // compare results
    for (size_t i = 0; i < traversal_filenames.size(); ++i)
    {
        std::string const output_filename{traversal_path_prefix.get_path().string() + traversal_filenames[i]};
        ASSERT_TRUE(std::filesystem::exists(output_filename));

        std::ifstream output_file{output_filename};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

        std::ifstream expected_file{(DATADIR + traversal_filenames[i])};
        std::string const expected_file_str((std::istreambuf_iterator<char>(expected_file)), std::istreambuf_iterator<char>());

        EXPECT_EQ(output_file_str, expected_file_str) << " failed at file " << traversal_filenames[i] << std::endl;
    }
}
