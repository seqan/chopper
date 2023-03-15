#include <fstream>
#include <ranges> // range comparisons
#include <string> // strings
#include <vector> // vectors

#include <cereal/archives/binary.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/join_with.hpp>

#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_layout)
{
    std::string const seq_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const taxa_filename{tmp_dir.path() / "data.tsv"};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename};
        fout << seq_filename << '\t' << "TAX1\n"
             << seq_filename << '\t'
             << "TAX2\n"
             /* << seq_filename << '\t' << "TAX2\n" */
             << seq_filename << '\t' << "TAX3\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--kmer-size",
                                         "15",
                                         "--threads",
                                         "2",
                                         "--input-file",
                                         taxa_filename.c_str(),
                                         "--tmax",
                                         "64",
                                         "--output-filename",
                                         binning_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string expected_file{"##CONFIG:\n"
                              "##{\n"
                              "##    \"config\": {\n"
                              "##        \"version\": 2,\n"
                              "##        \"data_file\": {\n"
                              "##            \"value0\": \""
                              + taxa_filename.string()
                              + "\"\n"
                                "##        },\n"
                                "##        \"debug\": false,\n"
                                "##        \"sketch_directory\": {\n"
                                "##            \"value0\": \"\"\n"
                                "##        },\n"
                                "##        \"k\": 15,\n"
                                "##        \"sketch_bits\": 12,\n"
                                "##        \"disable_sketch_output\": true,\n"
                                "##        \"precomputed_files\": false,\n"
                                "##        \"output_filename\": {\n"
                                "##            \"value0\": \""
                              + binning_filename.string()
                              + "\"\n"
                                "##        },\n"
                                "##        \"tmax\": 64,\n"
                                "##        \"num_hash_functions\": 2,\n"
                                "##        \"false_positive_rate\": 0.05,\n"
                                "##        \"alpha\": 1.2,\n"
                                "##        \"max_rearrangement_ratio\": 0.5,\n"
                                "##        \"threads\": 2,\n"
                                "##        \"disable_estimate_union\": false,\n"
                                "##        \"disable_rearrangement\": true,\n"
                                "##        \"determine_best_tmax\": false,\n"
                                "##        \"force_all_binnings\": false\n"
                                "##    }\n"
                                "##}\n"
                                "##ENDCONFIG\n"
                                "#HIGH_LEVEL_IBF max_bin_id:22\n"
                                "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
                              + seq_filename + "\t0\t22\n" + seq_filename + /* ";" + seq_filename + */ "\t22\t21\n"
                              + seq_filename + "\t43\t21\n"};

    ASSERT_TRUE(std::filesystem::exists(binning_filename));

    {
        std::ifstream output_file{binning_filename};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)),
                                          std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}

TEST_F(cli_test, chopper_layout2)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const taxa_filename{tmp_dir.path() / "data.tsv"};
    std::filesystem::path const binning_filename{tmp_dir.path() / "output.binning"};

    // we need to have filenames from the user
    {
        std::ofstream fout{taxa_filename};
        fout << seq1_filename << '\n' << seq2_filename << '\n' << seq3_filename << '\n' << seq4_filename << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--threads",
                                         "2",
                                         "--sketch-bits",
                                         "12",
                                         "--input-file",
                                         taxa_filename.c_str(),
                                         "--tmax",
                                         "64",
                                         "--output-filename",
                                         binning_filename.c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string expected_file{"##CONFIG:\n"
                              "##{\n"
                              "##    \"config\": {\n"
                              "##        \"version\": 2,\n"
                              "##        \"data_file\": {\n"
                              "##            \"value0\": \""
                              + taxa_filename.string()
                              + "\"\n"
                                "##        },\n"
                                "##        \"debug\": false,\n"
                                "##        \"sketch_directory\": {\n"
                                "##            \"value0\": \"\"\n"
                                "##        },\n"
                                "##        \"k\": 19,\n"
                                "##        \"sketch_bits\": 12,\n"
                                "##        \"disable_sketch_output\": true,\n"
                                "##        \"precomputed_files\": false,\n"
                                "##        \"output_filename\": {\n"
                                "##            \"value0\": \""
                              + binning_filename.string()
                              + "\"\n"
                                "##        },\n"
                                "##        \"tmax\": 64,\n"
                                "##        \"num_hash_functions\": 2,\n"
                                "##        \"false_positive_rate\": 0.05,\n"
                                "##        \"alpha\": 1.2,\n"
                                "##        \"max_rearrangement_ratio\": 0.5,\n"
                                "##        \"threads\": 2,\n"
                                "##        \"disable_estimate_union\": false,\n"
                                "##        \"disable_rearrangement\": false,\n"
                                "##        \"determine_best_tmax\": false,\n"
                                "##        \"force_all_binnings\": false\n"
                                "##    }\n"
                                "##}\n"
                                "##ENDCONFIG\n"
                                "#HIGH_LEVEL_IBF max_bin_id:54\n"
                                "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
                              + seq3_filename + "\t0\t15\n" + seq4_filename + "\t15\t24\n" + seq2_filename
                              + "\t39\t15\n" + seq1_filename + "\t54\t10\n"};

    ASSERT_TRUE(std::filesystem::exists(binning_filename));

    {
        std::ifstream output_file{binning_filename};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)),
                                          std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}
