#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/count/count_config.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/pack/chopper_pack.hpp>

#include "../api_test.hpp"
#include "print_debug_file.hpp"

TEST(chopper_pack_hll_test, few_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const pack_file{"pack.tsv"};
    seqan3::test::tmp_filename const hll_dir{"hll"};

    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    std::string const small2_filename = data("small2.fa");

    {
        count_config const config{.output_filename{count_file.get_path()},
                                  .hll_dir{hll_dir.get_path()},
                                  .exclusively_hlls{true}};

        std::vector<std::string> const files{seq1_filename,
                                             seq2_filename,
                                             seq3_filename,
                                             small_filename,
                                             small2_filename};

        robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters{};
        for (std::string const & filename : files)
            filename_clusters[filename].push_back(filename);

        count_kmers(filename_clusters, config);
    }

    char const * const argv[] = {"./chopper-pack",
                                 "-b", "4",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", pack_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};
    chopper_pack(pack_parser);

    std::vector<std::string> const expected_components
    {
        {"#HIGH_LEVEL_IBF max_bin_id:0"},
        {"#FILES\tBIN_INDICES\tNUMBER_OF_BINS"},
#if defined(__APPLE__)
        {seq1_filename + "\t48\t2"},
        {seq2_filename + "\t0\t48"},
        {seq3_filename + "\t50\t2"},
        {small_filename + "\t52\t6"},
        {small2_filename + "\t58\t6"}
#else
        {seq1_filename + "\t0\t48"},
        {seq2_filename + "\t50\t2"},
        {seq3_filename + "\t48\t2"},
        {small_filename + "\t58\t6"},
        {small2_filename + "\t52\t6"}
#endif
    };
    std::string const actual_file{string_from_file(pack_file.get_path())};

    size_t line_count{};
    for (auto && line : actual_file | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "Did not find " + line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);
}

TEST(chopper_pack_hll_test, many_ubs)
{
    seqan3::test::tmp_filename const count_file{"kmer_counts.tsv"};
    seqan3::test::tmp_filename const pack_file{"pack.tsv"};
    seqan3::test::tmp_filename const hll_dir{"hll"};

    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    std::string const small2_filename = data("small2.fa");

    {
        count_config const config{.output_filename{count_file.get_path()},
                                  .hll_dir{hll_dir.get_path()},
                                  .exclusively_hlls{true}};

        std::vector<std::string> const files{seq1_filename,
                                             seq2_filename,
                                             seq3_filename,
                                             small_filename,
                                             small2_filename};

        robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters{};
        std::string cluster_name{};
        for (size_t i{0}; i < 96u; ++i)
        {
            cluster_name = "cluster" + std::to_string(i);
            filename_clusters[cluster_name].push_back(files[i % 5]);
        }

        count_kmers(filename_clusters, config);
    }

    char const * const argv[] = {"./chopper-pack",
                                 "-b", "4",
                                 "-f", count_file.get_path().c_str(),
                                 "-o", pack_file.get_path().c_str()};
    int const argc = sizeof(argv) / sizeof(*argv);

    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};
    chopper_pack(pack_parser);

    std::vector<std::string> const expected_components
    {
        {"#HIGH_LEVEL_IBF max_bin_id:63"},
        {"#MERGED_BIN_11 max_bin_id:0"},
        {"#MERGED_BIN_12 max_bin_id:0"},
        {"#MERGED_BIN_13 max_bin_id:0"},
        {"#MERGED_BIN_14 max_bin_id:0"},
        {"#MERGED_BIN_15 max_bin_id:0"},
        {"#MERGED_BIN_16 max_bin_id:0"},
        {"#MERGED_BIN_17 max_bin_id:0"},
        {"#MERGED_BIN_18 max_bin_id:0"},
        {"#MERGED_BIN_19 max_bin_id:0"},
        {"#MERGED_BIN_20 max_bin_id:0"},
        {"#MERGED_BIN_21 max_bin_id:0"},
        {"#MERGED_BIN_22 max_bin_id:0"},
        {"#MERGED_BIN_23 max_bin_id:0"},
        {"#MERGED_BIN_24 max_bin_id:0"},
        {"#MERGED_BIN_25 max_bin_id:0"},
        {"#MERGED_BIN_26 max_bin_id:0"},
        {"#MERGED_BIN_63 max_bin_id:0"},
        {"#FILES\tBIN_INDICES\tNUMBER_OF_BINS"},
        {seq2_filename + "\t0\t1"},
        {seq3_filename + "\t1\t1"},
        {seq1_filename + "\t2\t1"},
        {seq1_filename + "\t3\t1"},
        {seq3_filename + "\t4\t1"},
        {seq2_filename + "\t5\t1"},
        {seq2_filename + "\t6\t1"},
        {seq1_filename + "\t7\t1"},
        {seq2_filename + "\t8\t1"},
        {seq3_filename + "\t9\t1"},
        {seq3_filename + "\t10\t1"},
        {seq2_filename + "\t11;0\t1;62"},
        {seq3_filename + "\t11;62\t1;2"},
        {seq1_filename + "\t12;0\t1;60"},
        {seq1_filename + "\t12;60\t1;2"},
        {seq2_filename + "\t12;62\t1;2"},
        {seq1_filename + "\t13;0\t1;60"},
        {seq2_filename + "\t13;60\t1;2"},
        {seq3_filename + "\t13;62\t1;2"},
        {seq3_filename + "\t14;0\t1;60"},
        {seq1_filename + "\t14;60\t1;2"},
        {seq2_filename + "\t14;62\t1;2"},
        {seq1_filename + "\t15;0\t1;60"},
        {seq2_filename + "\t15;60\t1;2"},
        {seq2_filename + "\t15;62\t1;2"},
        {seq3_filename + "\t16;0\t1;60"},
        {seq1_filename + "\t16;60\t1;2"},
        {seq1_filename + "\t16;62\t1;2"},
        {seq2_filename + "\t17;0\t1;60"},
        {seq3_filename + "\t17;60\t1;2"},
        {seq3_filename + "\t17;62\t1;2"},
        {seq1_filename + "\t18;0\t1;60"},
        {seq1_filename + "\t18;60\t1;2"},
        {seq2_filename + "\t18;62\t1;2"},
        {seq2_filename + "\t19;0\t1;60"},
        {seq1_filename + "\t19;60\t1;2"},
        {seq2_filename + "\t19;62\t1;2"},
        {seq3_filename + "\t20;0\t1;60"},
        {seq2_filename + "\t20;60\t1;2"},
        {seq1_filename + "\t20;62\t1;2"},
        {seq2_filename + "\t21;0\t1;60"},
        {seq3_filename + "\t21;60\t1;2"},
        {seq3_filename + "\t21;62\t1;2"},
        {seq3_filename + "\t22;0\t1;60"},
        {seq1_filename + "\t22;60\t1;2"},
        {seq2_filename + "\t22;62\t1;2"},
        {seq1_filename + "\t23;0\t1;60"},
        {seq3_filename + "\t23;60\t1;2"},
        {seq3_filename + "\t23;62\t1;2"},
        {seq3_filename + "\t24;0\t1;60"},
        {seq1_filename + "\t24;60\t1;2"},
        {seq1_filename + "\t24;62\t1;2"},
        {seq3_filename + "\t25;0\t1;60"},
        {seq1_filename + "\t25;60\t1;2"},
        {seq1_filename + "\t25;62\t1;2"},
        {seq2_filename + "\t26;0\t1;60"},
        {seq3_filename + "\t26;60\t1;2"},
        {seq2_filename + "\t26;62\t1;2"},
        {small2_filename + "\t27\t1"},
        {small2_filename + "\t28\t1"},
        {small2_filename + "\t29\t1"},
        {small_filename + "\t30\t1"},
        {small_filename + "\t31\t1"},
        {small2_filename + "\t32\t1"},
        {small2_filename + "\t33\t1"},
        {small_filename + "\t34\t1"},
        {small_filename + "\t35\t1"},
        {small_filename + "\t36\t1"},
        {small_filename + "\t37\t1"},
        {small_filename + "\t38\t1"},
        {small_filename + "\t39\t1"},
        {small2_filename + "\t40\t1"},
        {small2_filename + "\t41\t1"},
        {small2_filename + "\t42\t1"},
        {small2_filename + "\t43\t1"},
        {small2_filename + "\t44\t1"},
        {small2_filename + "\t45\t1"},
        {small_filename + "\t46\t1"},
        {small2_filename + "\t47\t1"},
        {small_filename + "\t48\t1"},
        {small2_filename + "\t49\t1"},
        {small_filename + "\t50\t1"},
        {small_filename + "\t51\t1"},
        {small_filename + "\t52\t1"},
        {small_filename + "\t53\t1"},
        {small2_filename + "\t54\t1"},
        {small2_filename + "\t55\t1"},
        {small_filename + "\t56\t1"},
        {small_filename + "\t57\t1"},
        {small2_filename + "\t58\t1"},
        {small_filename + "\t59\t1"},
        {small_filename + "\t60\t1"},
        {small2_filename + "\t61\t1"},
        {small_filename + "\t62\t1"},
        {small2_filename + "\t63;0\t1;58"},
        {small2_filename + "\t63;58\t1;6"}
    };
    std::string const actual_file{string_from_file(pack_file.get_path())};

    size_t line_count{};
    for (auto && line : actual_file | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "Did not find " + line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);
}
