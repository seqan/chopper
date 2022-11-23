#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, chopper_layout_statistics)
{
    seqan3::test::tmp_filename const input_prefix{"test"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{input_prefix.get_path().string() + ".count"};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    cli_test_result layout_result = execute_app("chopper"
                                                "--tmax", "64",
                                                "--input-prefix", input_prefix.get_path().c_str(),
                                                "--output-filename", layout_file.get_path().c_str(),
                                                "--output-verbose-statistics");

    std::string expected_cout =
R"expected_cout(## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
## uncorr_size : The expected size of an tmax-HIBF without FPR correction
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size	uncorr_size	level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
64	1.00	1.26	1.00	1.26	85.2KiB	45.4KiB	:0:1	:1:12	:37.0KiB:48.2KiB	:37.0KiB:8.3KiB	:64:768	:64:64	:81.25:100.00	:1:32	:1.00:17.45	:1.00:9.02	:1.00:6.50
)expected_cout";

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, expected_cout) << layout_result.out;
    EXPECT_EQ(layout_result.err, std::string{});
}

TEST_F(cli_test, chopper_layout_statistics_determine_best_bins)
{
    seqan3::test::tmp_filename const input_prefixname{"test"};
    seqan3::test::tmp_filename const binning_filename{"output.binning"};

    // Write count data to file.
    {
        std::ofstream fout{input_prefixname.get_path().string() + ".count"};
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
                                              "--tmax", "128",
                                              "--input-prefix", input_prefixname.get_path().c_str(),
                                              "--output-filename", binning_filename.get_path().c_str(),
                                              "--output-verbose-statistics",
                                              "--determine-best-tmax",
                                              "--force-all-binnings");


    std::string expected_cout =
R"expected_cout(## ### Parameters ###
## number of user bins = 10
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
## uncorr_size : The expected size of an tmax-HIBF without FPR correction
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size	uncorr_size	level	num_ibfs	level_size	level_size_no_corr	total_num_tbs	avg_num_tbs	split_tb_percentage	max_split_tb	avg_split_tb	max_factor	avg_factor
64	1.00	1.00	1.00	1.00	1.9MiB	1.8MiB	:0	:1	:1.9MiB	:1.8MiB	:64	:64	:100.00	:20	:6.40	:6.38	:3.62
128	0.96	0.96	1.63	1.56	3.1MiB	2.4MiB	:0	:1	:3.1MiB	:2.4MiB	:128	:128	:100.00	:47	:12.80	:12.17	:5.96
# Best t_max (regarding expected query runtime): 128
)expected_cout";

    EXPECT_EQ(layout_result.exit_code, 0);
    EXPECT_EQ(layout_result.out, expected_cout) << layout_result.out;
    EXPECT_EQ(layout_result.err, std::string{});
}
