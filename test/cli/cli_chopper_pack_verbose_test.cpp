#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, chopper_pack_summary)
{
    seqan3::test::tmp_filename const count_filename{"kmer_counts.txt"};
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    // Write count data to file.
    {
        std::ofstream fout{count_filename.get_path()};
        fout << "seq0\t10000\n"
                "seq1\t20000\n"
                "seq2\t30000\n"
                "seq3\t40000\n"
                "seq4\t50000\n"
                "seq5\t60000\n"
                "seq6\t70000\n"
                "seq7\t80000\n"
                "seq8\t90000\n"
                "seq9\t100000\n";
    }

    cli_test_result pack_result = execute_app("chopper", "pack",
                                              "-b", "128",
                                              "-f", count_filename.get_path().c_str(),
                                              "-o", binning_filename.get_path().c_str(),
                                              "--verbose",
                                              "--determine-num-bins",
                                              "--force-all-binnings");


    std::string expected_cout = 
        "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n"
        "64\t1.0000\t1.0000\n"
        "\n"
            "\tStatistics summary:\n"
            "\tlevel\tnum_ibfs\tlevel_size\tlevel_size_no_corr\ttotal_num_tbs\tavg_num_tbs\t"
                "split_tb_percentage\tmax_split_tb\tavg_split_tb\tmax_factor\tavg_factor\n"
            "\t0\t1\t1 MiB\t1 MiB\t64\t64\t100.0000\t20\t6.4000\t6.3777\t3.6157\n"
            "\tTotal HIBF size: 1 MiB\n"
            "\tTotal HIBF size no correction: 1 MiB\n"
        "\n"
        "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n"
        "128\t1.1000\t1.1000\n"
        "\n"
            "\tStatistics summary:\n"
            "\tlevel\tnum_ibfs\tlevel_size\tlevel_size_no_corr\ttotal_num_tbs\tavg_num_tbs\t"
            "split_tb_percentage\tmax_split_tb\tavg_split_tb\tmax_factor\tavg_factor\n"
            "\t0\t1\t3 MiB\t2 MiB\t128\t128\t100.0000\t47\t12.8000\t12.1721\t5.9559\n"
            "\tTotal HIBF size: 3 MiB\n"
            "\tTotal HIBF size no correction: 2 MiB\n"
        "\n"
        "Best t_max (regarding expected query runtime): 64\n";

    EXPECT_EQ(pack_result.exit_code, 0);
    EXPECT_EQ(pack_result.out, expected_cout);
    EXPECT_EQ(pack_result.err, std::string{});
}