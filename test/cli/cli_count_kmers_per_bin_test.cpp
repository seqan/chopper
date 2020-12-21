#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, small_example_only_high_level_bins)
{
    seqan3::test::tmp_filename traversal_filename{"traversal.out"};

    {
        std::ofstream fout{traversal_filename.get_path()};
        fout << data("small.fa").string() << "\tseq1\t0\t163\t0\t-\n"
             << data("small.fa").string() << "\tseq2\t0\t186\t0\t-\n"
             << data("small.fa").string() << "\tseq3\t0\t163\t0\t-\n"
             << data("small2.fa").string() << "\tseq10\t0\t163\t0\t-\n"
             << data("small2.fa").string() << "\tseq20\t0\t186\t0\t-\n"
             << data("small2.fa").string() << "\tseq30\t0\t163\t0\t-\n"
             << data("small.fa").string() << "\tseq1\t163\t247\t1\t-\n"
             << data("small.fa").string() << "\tseq2\t186\t327\t1\t-\n"
             << data("small.fa").string() << "\tseq3\t163\t284\t1\t-\n"
             << data("small2.fa").string() << "\tseq10\t163\t247\t1\t-\n"
             << data("small2.fa").string() << "\tseq20\t186\t327\t1\t-\n"
             << data("small2.fa").string() << "\tseq30\t163\t284\t1\t-\n"
             << data("small.fa").string() << "\tseq1\t247\t400\t2\t-\n"
             << data("small.fa").string() << "\tseq2\t327\t480\t2\t-\n"
             << data("small.fa").string() << "\tseq3\t284\t481\t2\t-\n"
             << data("small2.fa").string() << "\tseq10\t247\t400\t2\t-\n"
             << data("small2.fa").string() << "\tseq20\t327\t480\t2\t-\n"
             << data("small2.fa").string() << "\tseq30\t284\t481\t2\t-\n";
    };

    seqan3::test::tmp_filename output_filename{"small_traverse.out"};

    cli_test_result result = execute_app("count_kmers_per_bin",
                                         "-l", "0",
                                         "-f", traversal_filename.get_path().c_str());

    std::string expected_stdout
    {
        "0\t165\n"
        "1\t161\n"
        "2\t173\n"
    };

    // compare results
    EXPECT_EQ(result.out, expected_stdout);
    EXPECT_EQ(result.err, std::string{""});
}
