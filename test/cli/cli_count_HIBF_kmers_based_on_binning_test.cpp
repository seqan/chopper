#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, small_example)
{
    seqan3::test::tmp_filename binning_filename{"test.binning"};

    {
        std::ofstream fout{binning_filename.get_path()};
        fout << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n" // header is ignored anyway
             << "SPLIT_BIN_0\t" << data("small.fa").string() << "\t2\t500\n"
             << "MERGED_BIN_2_0\t" << data("small.fa").string() << "\t2\t2500\n"
             << "MERGED_BIN_2_1\t" << data("small.fa").string() << "\t2\t2500\n"
             << "SPLIT_BIN_3\t" << data("small.fa").string() + "\t3\t1000\n";
    }

    cli_test_result result = execute_app("count_HIBF_kmers_based_on_binning",
                                         "-k", "25",
                                         "-f", binning_filename.get_path().c_str());

    std::string expected_stderr
    {
        "SPLIT_BIN_0\t292\n"
        "SPLIT_BIN_0\t292\n"
        "SPLIT_BIN_3\t195\n"
        "SPLIT_BIN_3\t195\n"
        "SPLIT_BIN_3\t195\n"
        "MERGED_BIN_2\t585\n"
    };

    // compare results
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_stderr);
}
