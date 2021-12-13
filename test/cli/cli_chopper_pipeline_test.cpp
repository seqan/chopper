#include <fstream>
#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <cereal/archives/binary.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/to.hpp>

#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_pipeline)
{
    // CHOPPER COUNT
    // =========================================================================
    std::string const seq_filename = data("small.fa");
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    seqan3::test::tmp_filename const sketch_prefix{"small_sketch"};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq_filename << '\t' << "TAX1\n"
             << seq_filename << '\t' << "TAX2\n"
             /* << seq_filename << '\t' << "TAX2\n" */
             << seq_filename << '\t' << "TAX3\n";
    }

    cli_test_result count_result = execute_app("chopper", "count",
                                               "-k", "15",
                                               "-w", "25",
                                               "-t", "2",
                                               "-c", "2",
                                               "--disable-sketch-output",
                                               "-f", taxa_filename.get_path().c_str(),
                                               "-o", sketch_prefix.get_path().c_str());

    EXPECT_EQ(count_result.exit_code, 0);
    EXPECT_EQ(count_result.out, std::string{});
    EXPECT_EQ(count_result.err, std::string{});

    std::vector<std::string> expected_components
    {
        seq_filename + "\t86\tTAX3",
        seq_filename + /* ";" + seq_filename + */ "\t86\tTAX2",
        seq_filename + "\t86\tTAX1"
    };

    std::ifstream count_file{sketch_prefix.get_path().string() + ".count"};
    std::string const count_file_str((std::istreambuf_iterator<char>(count_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : count_file_str | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "missing:" << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);

    // Overwrite result file with expected order of elements.
    {
        std::ofstream fout{sketch_prefix.get_path().string() + ".count"};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::views::to<std::string>);
    }

    // CHOPPER LAYOUT
    // =========================================================================
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    cli_test_result layout_result = execute_app("chopper", "layout",
                                                "-b", "64",
                                                "-i", sketch_prefix.get_path().c_str(),
                                                "-o", binning_filename.get_path().c_str());

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, std::string{});
    EXPECT_EQ(layout_result.err, std::string{});

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:26\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n" +
        seq_filename + "\t0\t26\n" +
        seq_filename + /* ";" + seq_filename + */ "\t26\t19\n" +
        seq_filename + "\t45\t19\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_pipeline2)
{
    // CHOPPER COUNT
    // =========================================================================
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    seqan3::test::tmp_filename const sketch_prefix{"small_sketch"};

    // we need to have filenames from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq1_filename << '\n'
             << seq2_filename << '\n'
             << seq3_filename << '\n'
             << seq4_filename << '\n';
    }

    cli_test_result count_result = execute_app("chopper", "count",
                                               "-t", "2",
                                               "-s", "12",
                                               "-f", taxa_filename.get_path().c_str(),
                                               "-o", sketch_prefix.get_path().c_str());

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

    std::ifstream count_file{sketch_prefix.get_path().string() + ".count"};
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
        std::ofstream fout{sketch_prefix.get_path().string() + ".count"};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::views::to<std::string>);
    }

    // CHOPPER layout
    // =========================================================================
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    cli_test_result layout_result = execute_app("chopper", "layout",
                                                "-b", "64",
                                                "-t", "2",
                                                "-r",
                                                "-i", sketch_prefix.get_path().c_str(),
                                                "-o", binning_filename.get_path().c_str());

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, std::string{});
    EXPECT_EQ(layout_result.err, std::string{});

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n" +
        seq3_filename + "\t0\t54\n" +
        seq4_filename + "\t54\t6\n" +
        seq2_filename + "\t60\t2\n" +
        seq1_filename + "\t62\t2\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}
