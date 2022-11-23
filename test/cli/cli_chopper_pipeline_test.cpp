#include <fstream>
#include <ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <cereal/archives/binary.hpp>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/range/to.hpp>

#include "cli_test.hpp"

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_layout)
{
    std::string const seq_filename = data("small.fa");
    std::filesystem::path const sketch_prefix = "chopper_sketch";
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    // we need to have tax ids from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq_filename << '\t' << "TAX1\n"
             << seq_filename << '\t' << "TAX2\n"
             /* << seq_filename << '\t' << "TAX2\n" */
             << seq_filename << '\t' << "TAX3\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--kmer-size", "15",
                                         "--threads", "2",
                                         "--column-index", "2",
                                         "--disable-sketch-output",
                                         "--input-file", taxa_filename.get_path().c_str(),
                                         "--tmax", "64",
                                         "--output-filename", binning_filename.get_path().c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::vector<std::string> expected_components
    {
        seq_filename + "\t571\tTAX3",
        seq_filename + /* ";" + seq_filename + */ "\t571\tTAX2",
        seq_filename + "\t571\tTAX1"
    };

    std::ifstream count_file{sketch_prefix/* .get_path() */.string() + ".count"};
    std::string const count_file_str((std::istreambuf_iterator<char>(count_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : count_file_str | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "missing:" << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);

    // Overwrite result file with expected order of elements.
    {
        std::ofstream fout{sketch_prefix/* .get_path() */.string() + ".count"};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::ranges::to<std::string>());
    }

    std::string expected_file
    {
        "##CONFIG:\n"
        "##{\n"
        "##    \"config\": {\n"
        "##        \"version\": 1,\n"
        "##        \"input_prefix\": \"" + sketch_prefix/* .get_path() */.string() + "\",\n"
        "##        \"count_filename\": {\n"
        "##            \"value0\": \"" + sketch_prefix/* .get_path() */.string() + ".count\"\n"
        "##        },\n"
        "##        \"sketch_directory\": {\n"
        "##            \"value0\": \"" + sketch_prefix/* .get_path() */.string() + "_sketches\"\n"
        "##        },\n"
        "##        \"output_filename\": {\n"
        "##            \"value0\": \"" + binning_filename.get_path().string() + "\"\n"
        "##        },\n"
        "##        \"tmax\": 64,\n"
        "##        \"num_hash_functions\": 2,\n"
        "##        \"false_positive_rate\": 0.05,\n"
        "##        \"alpha\": 1.2,\n"
        "##        \"max_rearrangement_ratio\": 0.5,\n"
        "##        \"threads\": 2,\n"
        "##        \"estimate_union\": false,\n"
        "##        \"rearrange_user_bins\": false,\n"
        "##        \"determine_best_tmax\": false,\n"
        "##        \"force_all_binnings\": false\n"
        "##    }\n"
        "##}\n"
        "##ENDCONFIG\n"
        "#HIGH_LEVEL_IBF max_bin_id:22\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n" +
        seq_filename + "\t0\t22\n" +
        seq_filename + /* ";" + seq_filename + */ "\t22\t21\n" +
        seq_filename + "\t43\t21\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}

// check if each chopper submodule can work with the output of the other
TEST_F(cli_test, chopper_layout2)
{
    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const seq4_filename = data("small.fa");
    seqan3::test::tmp_filename const taxa_filename{"data.tsv"};
    // seqan3::test::tmp_filename const sketch_prefix{"small_sketch"};
    std::filesystem::path const sketch_prefix = "chopper_sketch";
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    // we need to have filenames from the user
    {
        std::ofstream fout{taxa_filename.get_path()};
        fout << seq1_filename << '\n'
             << seq2_filename << '\n'
             << seq3_filename << '\n'
             << seq4_filename << '\n';
    }

    cli_test_result result = execute_app("chopper",
                                         "--threads", "2",
                                         "--sketch-bits", "12",
                                         "--input-file", taxa_filename.get_path().c_str(),
                                         "--tmax", "64",
                                         "--rearrange-user-bins",
                                         "--output-filename", binning_filename.get_path().c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::vector<std::string> expected_components
    {
        seq4_filename + "\t588\t" + seq4_filename,
        seq1_filename + "\t383\t" + seq1_filename,
        seq2_filename + "\t467\t" + seq2_filename,
        seq3_filename + "\t468\t" + seq3_filename
    };

    std::ifstream count_file{sketch_prefix/* .get_path() */.string() + ".count"};
    std::string const count_file_str((std::istreambuf_iterator<char>(count_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : count_file_str | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "Missing:" << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);

    // Overwrite result file with expected order of elements.
    {
        std::ofstream fout{sketch_prefix/* .get_path() */.string() + ".count"};
        fout << (expected_components | seqan3::views::join_with(std::string{'\n'}) | seqan3::ranges::to<std::string>());
    }

    std::string expected_file
    {
        "##CONFIG:\n"
        "##{\n"
        "##    \"config\": {\n"
        "##        \"version\": 1,\n"
        "##        \"input_prefix\": \"" + sketch_prefix/* .get_path() */.string() + "\",\n"
        "##        \"count_filename\": {\n"
        "##            \"value0\": \"" + sketch_prefix/* .get_path() */.string() + ".count\"\n"
        "##        },\n"
        "##        \"sketch_directory\": {\n"
        "##            \"value0\": \"" + sketch_prefix/* .get_path() */.string() + "_sketches\"\n"
        "##        },\n"
        "##        \"output_filename\": {\n"
        "##            \"value0\": \"" + binning_filename.get_path().string() + "\"\n"
        "##        },\n"
        "##        \"tmax\": 64,\n"
        "##        \"num_hash_functions\": 2,\n"
        "##        \"false_positive_rate\": 0.05,\n"
        "##        \"alpha\": 1.2,\n"
        "##        \"max_rearrangement_ratio\": 0.5,\n"
        "##        \"threads\": 2,\n"
        "##        \"estimate_union\": true,\n"
        "##        \"rearrange_user_bins\": true,\n"
        "##        \"determine_best_tmax\": false,\n"
        "##        \"force_all_binnings\": false\n"
        "##    }\n"
        "##}\n"
        "##ENDCONFIG\n"
        "#HIGH_LEVEL_IBF max_bin_id:56\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n" +
        seq3_filename + "\t0\t14\n" +
        seq4_filename + "\t14\t29\n" +
        seq2_filename + "\t43\t13\n" +
        seq1_filename + "\t56\t8\n"
    };

    ASSERT_TRUE(std::filesystem::exists(binning_filename.get_path()));

    {
        std::ifstream output_file{binning_filename.get_path()};
        std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        ASSERT_EQ(output_file_str, expected_file);
    }
}
