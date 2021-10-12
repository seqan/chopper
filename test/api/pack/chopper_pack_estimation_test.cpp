#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/chopper_pack.hpp>

#include "../api_test.hpp"

TEST(chopper_pack_estimation_test, few_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const pack_file{"pack.tsv"};

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

    char const * const argv[] = {"./chopper-pack",
                                 "-b", "4",
                                 "--determine-num-bins",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", pack_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    chopper_pack(pack_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.0000"
                                                      "\t1.0000\nBest t_max (regarding expected query runtime): 64\n");
}

TEST(chopper_pack_estimation_test, many_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const pack_file{"pack.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{count_file.get_path()};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    char const * const argv[] = {"./chopper-pack",
                                 "-b", "1024",
                                 "--determine-num-bins",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", pack_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    chopper_pack(pack_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.0000"
                                                      "\t1.2571\n128\t1.1000\t1.1214\n256\t1.3200\t1.3200\nBest t_max "
                                                      "(regarding expected query runtime): 128\n");
}

TEST(chopper_pack_estimation_test, many_ubs_force_all)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const pack_file{"pack.tsv"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{count_file.get_path()};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    char const * const argv[] = {"./chopper-pack",
                                 "-b", "256",
                                 "--determine-num-bins",
                                 "--force-all-binnings",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", pack_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};
    testing::internal::CaptureStdout();
    chopper_pack(pack_parser);

    EXPECT_EQ(testing::internal::GetCapturedStdout(), "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n64\t1.0000"
                                                      "\t1.2571\n128\t1.1000\t1.1214\n256\t1.3200\t1.3200\nBest t_max "
                                                      "(regarding expected query runtime): 128\n");
}
