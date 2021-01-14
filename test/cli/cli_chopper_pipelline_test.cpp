#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

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
        std::string{merged_bin_prefix} + "_1_32\t" + seq_filename + ";" + seq_filename + "\t32\t3\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        EXPECT_EQ(output_file_str, expected_file);
    }

    // CHOPPER SPLIT
    // =========================================================================
    seqan3::test::tmp_filename const chopper_split_filename{"small.split"};

    cli_test_result split_result = execute_app("chopper", "split",
                                               "-k", "15",
                                               "-w", "25",
                                               "-f", binning_filename.get_path().c_str(),
                                               "-o", chopper_split_filename.get_path().c_str());

    // compare results
    ASSERT_TRUE(std::filesystem::exists(chopper_split_filename.get_path()));
    ASSERT_TRUE(std::filesystem::exists(data("small.split")));

    std::ifstream output_file{chopper_split_filename.get_path()};
    std::filesystem::path expected_output_filename{data("small.split")};
    std::ifstream expected_chopper_split_file{expected_output_filename};

    // note that file small.split does not contain the full filename (with directory)
    // since the directory changes on every system
    std::string const directory{expected_output_filename.parent_path().string() + "/"};
    std::string output_line;
    std::string expected_line;
    while (std::getline(output_file, output_line) && std::getline(expected_chopper_split_file, expected_line))
    {
        std::string full_expected_line{(expected_line[0] == '#') ? expected_line : directory + expected_line};
        EXPECT_EQ(output_line, full_expected_line);
    }

    // CHOPPER BUILD from split
    // =========================================================================
    {
        seqan3::test::tmp_filename const chopper_build_dir{"build/"};

        cli_test_result build_result = execute_app("chopper", "build",
                                                   "--kmer-size", "15",
                                                   "--false-positive-rate", "0.01",
                                                   "--overlap", "20",
                                                   "-s", chopper_split_filename.get_path().c_str(),
                                                   "-o", chopper_build_dir.get_path().c_str());

        auto high_level_ibf_path = chopper_build_dir.get_path()/"high_level.ibf";
        auto low_level_ibf_path = chopper_build_dir.get_path()/"low_level_1.ibf";

        ASSERT_TRUE(std::filesystem::exists(high_level_ibf_path));
        ASSERT_TRUE(std::filesystem::exists(low_level_ibf_path));

        seqan3::interleaved_bloom_filter<> high_level_ibf;
        seqan3::interleaved_bloom_filter<> low_level_ibf;

        {   // Load High Level IBF
            std::ifstream is(high_level_ibf_path, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(high_level_ibf);
        }
        {   // Load Low Level IBF
            std::ifstream is(low_level_ibf_path, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(low_level_ibf);
        }

        EXPECT_EQ(high_level_ibf.bin_count(), 2);
        EXPECT_EQ(low_level_ibf.bin_count(), 64);
    }

    // CHOPPER BUILD from pack
    // =========================================================================
    {
        seqan3::test::tmp_filename const chopper_build_dir{"build/"};

        cli_test_result build_result = execute_app("chopper", "build",
                                                   "--kmer-size", "15",
                                                   "--false-positive-rate", "0.01",
                                                   "-p", binning_filename.get_path().c_str(),
                                                   "-o", chopper_build_dir.get_path().c_str());

        auto high_level_ibf_path = chopper_build_dir.get_path()/"high_level.ibf";
        auto low_level_ibf_path = chopper_build_dir.get_path()/"low_level_1.ibf";

        ASSERT_TRUE(std::filesystem::exists(high_level_ibf_path));
        ASSERT_TRUE(std::filesystem::exists(low_level_ibf_path));

        seqan3::interleaved_bloom_filter<> high_level_ibf;
        seqan3::interleaved_bloom_filter<> low_level_ibf;

        {   // Load High Level IBF
            std::ifstream is(high_level_ibf_path, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(high_level_ibf);
        }
        {   // Load Low Level IBF
            std::ifstream is(low_level_ibf_path, std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(low_level_ibf);
        }

        EXPECT_EQ(high_level_ibf.bin_count(), 2);
        EXPECT_EQ(low_level_ibf.bin_count(), 64);
    }
}
