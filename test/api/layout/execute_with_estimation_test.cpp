#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/layout/execute.hpp>

#include "../api_test.hpp"

TEST(execute_estimation_test, few_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    {
        std::ofstream fout{count_file.get_path()};
        fout << "seq0\t500\n"
             << "seq1\t1000\n"
             << "seq2\t500\n"
             << "seq3\t500\n"
             << "seq4\t500\n"
             << "seq5\t500\n"
             << "seq6\t500\n"
             << "seq7\t500\n";
    }

    char const * const argv[] = {"./chopper-layout",
                                 "-b", "4",
                                 "--determine-num-bins",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", layout_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser layout_parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    chopper::layout::execute(layout_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.00"
                                                      "\t1.00\n#Best t_max (regarding expected query runtime):64\n");
    EXPECT_EQ(testing::internal::GetCapturedStderr(), "[CHOPPER LAYOUT WARNING]: Your requested number of technical "
                                                      "bins was not a multiple of 64. Due to the architecture of the "
                                                      "HIBF, it will use up space equal to the next multiple of 64 "
                                                      "anyway, so we increased your number of technical bins to 64.\n");
}

TEST(execute_estimation_test, many_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{count_file.get_path()};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    char const * const argv[] = {"./chopper-layout",
                                 "-b", "1024",
                                 "--determine-num-bins",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", layout_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser layout_parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    chopper::layout::execute(layout_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.00"
                                                      "\t1.26\n128\t1.10\t1.12\n256\t1.32\t1.32\n#Best t_max "
                                                      "(regarding expected query runtime):128\n");
}

TEST(execute_estimation_test, many_ubs_force_all)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{count_file.get_path()};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    char const * const argv[] = {"./chopper-layout",
                                 "-b", "256",
                                 "--determine-num-bins",
                                 "--force-all-binnings",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", layout_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser layout_parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    chopper::layout::execute(layout_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.00"
                                                      "\t1.26\n128\t1.10\t1.12\n256\t1.32\t1.32\n#Best t_max "
                                                      "(regarding expected query runtime):128\n");
}
