#include <fstream>
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <cereal/archives/binary.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/to.hpp>

#include <chopper/detail_bin_prefixes.hpp>

#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_pipeline)
{
    // CHOPPER COUNT
    // =========================================================================
    std::string const seq_filename = data("small.fa");
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    seqan3::test::tmp_filename const count_filename{"kmer_counts.txt"};

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
                                               "-f", taxa_filename.get_path().c_str(),
                                               "-o", count_filename.get_path().c_str());

    EXPECT_EQ(count_result.exit_code, 0);
    EXPECT_EQ(count_result.out, std::string{});
    EXPECT_EQ(count_result.err, std::string{});

    std::vector<std::string> expected_components
    {
        seq_filename + "\t95\tTAX3",
        seq_filename + ";" + seq_filename + "\t95\tTAX2",
        seq_filename + "\t95\tTAX1"
    };

    std::ifstream count_file{count_filename.get_path()};
    std::string const count_file_str((std::istreambuf_iterator<char>(count_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : count_file_str | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end());
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);

    // Overwrite result file with expected order of elements.
    {
        std::ofstream fout{count_filename.get_path()};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::views::to<std::string>);
    }

    // CHOPPER PACK
    // =========================================================================
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    cli_test_result pack_result = execute_app("chopper", "pack",
                                              "-b", "2",
                                              "-f", count_filename.get_path().c_str(),
                                              "-o", binning_filename.get_path().c_str());

    EXPECT_EQ(pack_result.exit_code, 0);
    EXPECT_EQ(pack_result.out, std::string{});
    EXPECT_EQ(pack_result.err, std::string{});

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:1\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n" +
        seq_filename + "\t0\t1\t95\n" +
        seq_filename + "\t1;0\t1;32\t190;3\n" +
        seq_filename + ";" + seq_filename + "\t1;32\t1;32\t190;3\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }

    // CHOPPER SPLIT
    // =========================================================================
    seqan3::test::tmp_filename const chopper_split_filename{"small.split"};

    cli_test_result split_result = execute_app("chopper", "split",
                                               "-k", "15",
                                               "-w", "25",
                                               "-f", binning_filename.get_path().c_str(),
                                               "-o", chopper_split_filename.get_path().c_str());

    EXPECT_EQ(split_result.exit_code, 0);
    EXPECT_EQ(split_result.out, std::string{});
    EXPECT_EQ(split_result.err, std::string{});

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
        ASSERT_EQ(output_line, full_expected_line);
    }

    // CHOPPER BUILD from split
    // =========================================================================
    // {
    //     seqan3::test::tmp_filename output_path{"chopper.test.index"};

    //     cli_test_result build_result = execute_app("chopper", "build",
    //                                                "--kmer-size", "15",
    //                                                "--false-positive-rate", "0.01",
    //                                                "--overlap", "20",
    //                                                "-s", chopper_split_filename.get_path().c_str(),
    //                                                "-o", output_path.get_path().c_str());

    // EXPECT_EQ(build_result.exit_code, 0);
    // EXPECT_EQ(build_result.out, std::string{});
    // EXPECT_EQ(build_result.err, std::string{});

    //     ASSERT_TRUE(std::filesystem::exists(output_path.get_path()));

    //     std::vector<seqan3::interleaved_bloom_filter<>> hibf;

    //     {
    //         std::ifstream is(output_path.get_path(), std::ios::binary);
    //         cereal::BinaryInputArchive archive(is);
    //         archive(hibf);
    //     }

    //     ASSERT_EQ(hibf.size(), 2);
    //     EXPECT_EQ(hibf[0].bin_count(), 2);
    //     EXPECT_EQ(hibf[1].bin_count(), 64);
    // }

    // CHOPPER BUILD from pack
    // =========================================================================
    {
        seqan3::test::tmp_filename output_path{"chopper.test.index"};

        cli_test_result build_result = execute_app("chopper", "build",
                                                   "--kmer-size", "15",
                                                   "--false-positive-rate", "0.01",
                                                   "-p", binning_filename.get_path().c_str(),
                                                   "-o", output_path.get_path().c_str());

        EXPECT_EQ(build_result.exit_code, 0);
        EXPECT_EQ(build_result.out, std::string{});
        EXPECT_EQ(build_result.err, std::string{});

        ASSERT_TRUE(std::filesystem::exists(output_path.get_path()));

        std::vector<seqan3::interleaved_bloom_filter<>> hibf;

        {
            std::ifstream is(output_path.get_path(), std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(hibf);
        }

        ASSERT_EQ(hibf.size(), 2);
        EXPECT_EQ(hibf[0].bin_count(), 2);
        EXPECT_EQ(hibf[1].bin_count(), 64);
    }
}

class clear_directory
{
public:
    clear_directory() = delete;
    clear_directory(clear_directory const &) = delete;
    clear_directory & operator=(clear_directory const &) = delete;
    clear_directory(clear_directory &&) = default;
    clear_directory & operator=(clear_directory &&) = default;

    explicit clear_directory(std::filesystem::path const & path) : directory_path(path) {}

    ~clear_directory()
    {
        for (auto && path : std::filesystem::directory_iterator(directory_path))
            std::filesystem::remove_all(path);
    }

private:
    std::filesystem::path const directory_path{};
};

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_hll_pipeline)
{
    // CHOPPER COUNT
    // =========================================================================
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    seqan3::test::tmp_filename const count_filename{"kmer_counts.txt"};
    seqan3::test::tmp_directory const hll_dir{};
    clear_directory clear_hll{hll_dir.path()};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq1_filename << '\n'
             << seq2_filename << '\n'
             << seq3_filename << '\n'
             << seq4_filename << '\n';
    }

    cli_test_result count_result = execute_app("chopper", "count",
                                               "-e",
                                               "-t", "2",
                                               "-s", "12",
                                               "-d", hll_dir.path().c_str(),
                                               "-f", taxa_filename.get_path().c_str(),
                                               "-o", count_filename.get_path().c_str());

    EXPECT_EQ(count_result.exit_code, 0);
    EXPECT_EQ(count_result.out, std::string{});
    EXPECT_EQ(count_result.err, std::string{});

    std::vector<std::string> expected_components
    {
        seq4_filename + "\t2\t" + seq4_filename,
        seq1_filename + "\t1\t" + seq1_filename,
        seq2_filename + "\t1\t" + seq2_filename,
        seq3_filename + "\t1\t" + seq3_filename
    };

    std::ifstream count_file{count_filename.get_path()};
    std::string const count_file_str((std::istreambuf_iterator<char>(count_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : count_file_str | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end());
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);

    // Overwrite result file with expected order of elements.
    {
        std::ofstream fout{count_filename.get_path()};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::views::to<std::string>);
    }

    // CHOPPER PACK
    // =========================================================================
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    cli_test_result pack_result = execute_app("chopper", "pack",
                                              "-b", "2",
                                              "-t", "2",
                                              "-r",
                                              "-d", hll_dir.path().c_str(),
                                              "-f", count_filename.get_path().c_str(),
                                              "-o", binning_filename.get_path().c_str());

    EXPECT_EQ(pack_result.exit_code, 0);
    EXPECT_EQ(pack_result.out, std::string{});
    EXPECT_EQ(pack_result.err, std::string{});

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:0\n"
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n" +
        seq4_filename + "\t0;0\t1;62\t3;1\n" +
        seq3_filename + "\t0;62\t1;2\t3;1\n" +
        seq1_filename + "\t1;0\t1;62\t2;1\n" +
        seq2_filename + "\t1;62\t1;2\t2;1\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }

    // CHOPPER BUILD from pack
    // =========================================================================
    {
        seqan3::test::tmp_filename output_path{"chopper.test.index"};

        cli_test_result build_result = execute_app("chopper", "build",
                                                   "--kmer-size", "15",
                                                   "--false-positive-rate", "0.01",
                                                   "-p", binning_filename.get_path().c_str(),
                                                   "-o", output_path.get_path().c_str());

        EXPECT_EQ(build_result.exit_code, 0);
        EXPECT_EQ(build_result.out, std::string{});
        EXPECT_EQ(build_result.err, std::string{});

        ASSERT_TRUE(std::filesystem::exists(output_path.get_path()));

        std::vector<seqan3::interleaved_bloom_filter<>> hibf;

        {
            std::ifstream is(output_path.get_path(), std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(hibf);
        }

        ASSERT_EQ(hibf.size(), 3);
        EXPECT_EQ(hibf[0].bin_count(), 2);
        EXPECT_EQ(hibf[1].bin_count(), 64);
        EXPECT_EQ(hibf[2].bin_count(), 64);
    }
}
