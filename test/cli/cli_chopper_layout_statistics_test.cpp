#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, chopper_layout_statistics)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{count_file.get_path()};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    cli_test_result layout_result = execute_app("chopper", "layout",
                                              "-b", "64",
                                              "-f", count_file.get_path().c_str(),
                                              "-o", layout_file.get_path().c_str(),
                                              "--output-statistics");

    std::string expected_cout =
R"expected_cout(#number of user bins:96
#number of hash functions:2
#false positive rate:0.05

level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
0	1	37 KiB	37 KiB	64	64	81.25	1	1.00	1.00	1.00
1	12	48 KiB	8 KiB	768	64	100.00	32	17.45	9.02	6.50
#Total HIBF size:85 KiB
#Total HIBF size no correction:45 KiB

)expected_cout";

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, expected_cout) << layout_result.out;
    EXPECT_EQ(layout_result.err, std::string{});
}

TEST_F(cli_test, chopper_layout_statistics_determine_best_bins)
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

    cli_test_result layout_result = execute_app("chopper", "layout",
                                              "-b", "128",
                                              "-f", count_filename.get_path().c_str(),
                                              "-o", binning_filename.get_path().c_str(),
                                              "--output-statistics",
                                              "--determine-num-bins",
                                              "--force-all-binnings");


    std::string expected_cout =
R"expected_cout(#number of user bins:10
#number of hash functions:2
#false positive rate:0.05

#T_Max:64
#C_{T_Max}:1.00
#relative expected HIBF query time cost (l):1.00
#relative HIBF memory usage (m):1.00
#l*m:1.00
level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
0	1	1 MiB	1 MiB	64	64	100.00	20	6.40	6.38	3.62
#Total HIBF size:1 MiB
#Total HIBF size no correction:1 MiB

#T_Max:128
#C_{T_Max}:1.10
#relative expected HIBF query time cost (l):1.10
#relative HIBF memory usage (m):1.63
#l*m:1.79
level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
0	1	3 MiB	2 MiB	128	128	100.00	47	12.80	12.17	5.96
#Total HIBF size:3 MiB
#Total HIBF size no correction:2 MiB

#Best t_max (regarding expected query runtime):64
)expected_cout";

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, expected_cout) << layout_result.out;
    EXPECT_EQ(layout_result.err, std::string{});
}
