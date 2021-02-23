#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("chopper");
    std::string expected
    {
        "chopper\n"
        "=======\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_out_file)
{
    seqan3::test::tmp_filename output_filename{"small.split"};

    cli_test_result result = execute_app("chopper", "split",
                                         "-k", "15", "-w", "25", "-b", "3",
                                         "-s", data("small.fa"),
                                         "-s", data("small2.fa"),
                                         "-o", output_filename.get_path().c_str());

    // compare results
    std::string expected_file_str
    {
        "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tIBF_BIN_INDICES\n" +
        data("small.fa").string() + "\tseq1\t0\t163\t0\n" +
        data("small.fa").string() + "\tseq2\t0\t186\t0\n" +
        data("small.fa").string() + "\tseq3\t0\t163\t0\n" +
        data("small2.fa").string() + "\tseq10\t0\t163\t0\n" +
        data("small2.fa").string() + "\tseq20\t0\t186\t0\n" +
        data("small2.fa").string() + "\tseq30\t0\t163\t0\n" +
        data("small.fa").string() + "\tseq1\t163\t247\t1\n" +
        data("small.fa").string() + "\tseq2\t186\t327\t1\n" +
        data("small.fa").string() + "\tseq3\t163\t284\t1\n" +
        data("small2.fa").string() + "\tseq10\t163\t247\t1\n" +
        data("small2.fa").string() + "\tseq20\t186\t327\t1\n" +
        data("small2.fa").string() + "\tseq30\t163\t284\t1\n" +
        data("small.fa").string() + "\tseq1\t247\t400\t2\n" +
        data("small.fa").string() + "\tseq2\t327\t480\t2\n" +
        data("small.fa").string() + "\tseq3\t284\t481\t2\n" +
        data("small2.fa").string() + "\tseq10\t247\t400\t2\n" +
        data("small2.fa").string() + "\tseq20\t327\t480\t2\n" +
        data("small2.fa").string() + "\tseq30\t284\t481\t2\n"
    };

    // compare results
    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file_str);
}
