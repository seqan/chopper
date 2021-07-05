#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, small_example)
{
    seqan3::test::tmp_filename binning_filename{"test.binning"};
    seqan3::test::tmp_filename counts_filename{"kmer_counts.txt"};
    seqan3::test::tmp_filename output_filename{"out.txt"};

    {
        std::ofstream fout{counts_filename.get_path()};
        fout << data("small.fa").string() << '\t' << "585" << '\t'
             << data("small.fa").string() << '\n';
    }

    {
        std::ofstream fout{binning_filename.get_path()};
        fout << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n" // header is ignored anyway
             << "SPLIT_BIN_0\t" << data("small.fa").string() << "\t2\t500\n"
             << "MERGED_BIN_2_0\t" << data("small.fa").string() << "\t2\t2500\n"
             << "MERGED_BIN_2_1\t" << data("small.fa").string() << "\t2\t2500\n"
             << "SPLIT_BIN_3\t" << data("small.fa").string() + "\t3\t1000\n";
    }

    cli_test_result result = execute_app("count_HIBF_kmers_based_on_binning",
                                         "-c", counts_filename.get_path().c_str(),
                                         "-k", "25",
                                         "-t", "1",
                                         "-b", binning_filename.get_path().c_str(),
                                         "-o", output_filename.get_path().c_str());

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});

    std::string expected_stdout
    {
        "SPLIT_BIN_0\t292\t0\n"
        "SPLIT_BIN_0\t292\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "MERGED_BIN_2\t585\t1170\n"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    // compare results
    EXPECT_EQ(expected_stdout, output_file_str);

    // peak memory usage
    // EXPECT_EQ(result.err, std::string{});
}
