#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>

#include "../api_test.hpp"
#include "print_debug_file.hpp"

TEST(execute_hll_test, few_ubs)
{
    seqan3::test::tmp_filename const data_file{"test.tsv"};
    seqan3::test::tmp_filename const io_prefix{"hll"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};

    std::string const seq1_filename = data("seq1.fa");
    std::string const seq2_filename = data("seq2.fa");
    std::string const seq3_filename = data("seq3.fa");
    std::string const small_filename = data("small.fa");
    std::string const small2_filename = data("small2.fa");

    {
        chopper::configuration config{.data_file{data_file.get_path()},
                                             .output_prefix{io_prefix.get_path().string()}};
        chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

        std::vector<std::string> const files{seq1_filename,
                                             seq2_filename,
                                             seq3_filename,
                                             small_filename,
                                             small2_filename};

        robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters{};
        for (std::string const & filename : files)
            filename_clusters[filename].push_back(filename);

        chopper::count::count_kmers(filename_clusters, config);
    }

    {
        std::ofstream fout{io_prefix.get_path().string() + ".count"};
        fout << seq2_filename << "\t1\t" << seq2_filename << '\n'
             << small_filename << "\t2\t" << small_filename << '\n'
             << small2_filename << "\t2\t" << small2_filename << '\n'
             << seq3_filename << "\t1\t" << seq3_filename << '\n'
             << seq1_filename << "\t1\t" << seq1_filename << '\n';
    }

    chopper::configuration config{};
    config.tmax = 64;
    config.input_prefix = io_prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = layout_file.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    // char const * const argv[] = {"./chopper-layout",
    //                              "--tmax", "64",
    //                              "--input-prefix", io_prefix.get_path().c_str(),
    //                              "--output-filename", layout_file.get_path().c_str()};

    chopper::layout::execute(config);

    std::vector<std::string> const expected_components
    {
        {"##CONFIG:"},
        {"##{"},
        {"##    \"config\": {"},
        {"##        \"version\": 1,"},
        {"##        \"input_prefix\": \"" + io_prefix.get_path().string() + "\","},
        {"##        \"count_filename\": {"},
        {"##            \"value0\": \"" + io_prefix.get_path().string() + ".count\""},
        {"##        },"},
        {"##        \"sketch_directory\": {"},
        {"##            \"value0\": \"" + io_prefix.get_path().string() + "_sketches\""},
        {"##        },"},
        {"##        \"output_filename\": {"},
        {"##            \"value0\": \"" + layout_file.get_path().string() + "\""},
        {"##        },"},
        {"##        \"tmax\": 64,"},
        {"##        \"num_hash_functions\": 2,"},
        {"##        \"false_positive_rate\": 0.05,"},
        {"##        \"alpha\": 1.2,"},
        {"##        \"max_rearrangement_ratio\": 0.5,"},
        {"##        \"threads\": 1,"},
        {"##        \"estimate_union\": false,"},
        {"##        \"rearrange_user_bins\": false,"},
        {"##        \"determine_best_tmax\": false,"},
        {"##        \"force_all_binnings\": false"},
        {"##    }"},
        {"##}"},
        {"##ENDCONFIG"},
        {"#HIGH_LEVEL_IBF max_bin_id:0"},
        {"#FILES\tBIN_INDICES\tNUMBER_OF_BINS"},
        {seq1_filename + "\t0\t48"},
        {seq2_filename + "\t50\t2"},
        {seq3_filename + "\t48\t2"},
        {small_filename + "\t58\t6"},
        {small2_filename + "\t52\t6"}
    };

    std::string const actual_file{string_from_file(layout_file.get_path())};

    size_t line_count{};
    for (auto && line : actual_file | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "Did not find " + line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);
}

// TEST(execute_hll_test, many_ubs)
// {
//     seqan3::test::tmp_filename const data_file{"test.tsv"};
//     seqan3::test::tmp_filename const io_prefix{"hll"};
//     seqan3::test::tmp_filename const layout_file{"layout.tsv"};

//     std::string const seq1_filename = data("seq1.fa");
//     std::string const seq2_filename = data("seq2.fa");
//     std::string const seq3_filename = data("seq3.fa");
//     std::string const small_filename = data("small.fa");
//     std::string const small2_filename = data("small2.fa");

//     {
//         chopper::configuration config{.data_file{data_file.get_path()},
//                                              .output_prefix{io_prefix.get_path().string()}};
//         chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

//         std::vector<std::string> const files{seq1_filename,
//                                              seq2_filename,
//                                              seq3_filename,
//                                              small_filename,
//                                              small2_filename};

//         robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters{};
//         std::string cluster_name{};
//         for (size_t i{0}; i < 96u; ++i)
//         {
//             cluster_name = "cluster" + std::to_string(i);
//             filename_clusters[cluster_name].push_back(files[i % 5]);
//         }

//         chopper::count::count_kmers(filename_clusters, config);
//     }

//     {
//         std::ofstream fout{io_prefix.get_path().string() + ".count"};
//         fout << small_filename << "\t2\tcluster18\n"
//              << seq1_filename << "\t1\tcluster45\n"
//              << small2_filename << "\t2\tcluster44\n"
//              << small_filename << "\t2\tcluster23\n"
//              << small_filename << "\t2\tcluster8\n"
//              << small_filename << "\t2\tcluster63\n"
//              << seq2_filename << "\t1\tcluster61\n"
//              << seq1_filename << "\t1\tcluster15\n"
//              << small_filename << "\t2\tcluster13\n"
//              << small2_filename << "\t2\tcluster19\n"
//              << seq3_filename << "\t1\tcluster92\n"
//              << seq2_filename << "\t1\tcluster1\n"
//              << seq2_filename << "\t1\tcluster31\n"
//              << seq1_filename << "\t1\tcluster75\n"
//              << seq1_filename << "\t1\tcluster60\n"
//              << small2_filename << "\t2\tcluster59\n"
//              << seq1_filename << "\t1\tcluster10\n"
//              << seq3_filename << "\t1\tcluster42\n"
//              << seq3_filename << "\t1\tcluster62\n"
//              << small2_filename << "\t2\tcluster89\n"
//              << seq3_filename << "\t1\tcluster47\n"
//              << small_filename << "\t2\tcluster43\n"
//              << small2_filename << "\t2\tcluster64\n"
//              << seq2_filename << "\t1\tcluster76\n"
//              << seq2_filename << "\t1\tcluster26\n"
//              << seq1_filename << "\t1\tcluster25\n"
//              << seq1_filename << "\t1\tcluster65\n"
//              << seq2_filename << "\t1\tcluster16\n"
//              << seq1_filename << "\t1\tcluster35\n"
//              << seq2_filename << "\t1\tcluster36\n"
//              << small_filename << "\t2\tcluster53\n"
//              << small2_filename << "\t2\tcluster4\n"
//              << seq3_filename << "\t1\tcluster7\n"
//              << seq1_filename << "\t1\tcluster80\n"
//              << small_filename << "\t2\tcluster88\n"
//              << small_filename << "\t2\tcluster73\n"
//              << small2_filename << "\t2\tcluster39\n"
//              << small2_filename << "\t2\tcluster79\n"
//              << small2_filename << "\t2\tcluster74\n"
//              << seq1_filename << "\t1\tcluster90\n"
//              << seq3_filename << "\t1\tcluster12\n"
//              << seq2_filename << "\t1\tcluster91\n"
//              << seq2_filename << "\t1\tcluster86\n"
//              << seq1_filename << "\t1\tcluster55\n"
//              << seq2_filename << "\t1\tcluster41\n"
//              << seq3_filename << "\t1\tcluster37\n"
//              << seq3_filename << "\t1\tcluster22\n"
//              << seq3_filename << "\t1\tcluster82\n"
//              << seq1_filename << "\t1\tcluster30\n"
//              << small_filename << "\t2\tcluster93\n"
//              << small_filename << "\t2\tcluster33\n"
//              << seq2_filename << "\t1\tcluster66\n"
//              << seq2_filename << "\t1\tcluster81\n"
//              << seq1_filename << "\t1\tcluster40\n"
//              << seq1_filename << "\t1\tcluster70\n"
//              << small2_filename << "\t2\tcluster54\n"
//              << seq3_filename << "\t1\tcluster87\n"
//              << seq2_filename << "\t1\tcluster11\n"
//              << small2_filename << "\t2\tcluster49\n"
//              << small2_filename << "\t2\tcluster24\n"
//              << small2_filename << "\t2\tcluster69\n"
//              << small_filename << "\t2\tcluster83\n"
//              << seq1_filename << "\t1\tcluster0\n"
//              << seq1_filename << "\t1\tcluster50\n"
//              << small_filename << "\t2\tcluster38\n"
//              << small2_filename << "\t2\tcluster9\n"
//              << seq2_filename << "\t1\tcluster6\n"
//              << seq3_filename << "\t1\tcluster17\n"
//              << seq3_filename << "\t1\tcluster72\n"
//              << seq3_filename << "\t1\tcluster2\n"
//              << seq1_filename << "\t1\tcluster20\n"
//              << seq2_filename << "\t1\tcluster21\n"
//              << seq2_filename << "\t1\tcluster46\n"
//              << seq3_filename << "\t1\tcluster32\n"
//              << seq3_filename << "\t1\tcluster77\n"
//              << seq2_filename << "\t1\tcluster51\n"
//              << seq3_filename << "\t1\tcluster27\n"
//              << small_filename << "\t2\tcluster48\n"
//              << seq3_filename << "\t1\tcluster57\n"
//              << seq3_filename << "\t1\tcluster67\n"
//              << small_filename << "\t2\tcluster68\n"
//              << seq3_filename << "\t1\tcluster52\n"
//              << small2_filename << "\t2\tcluster94\n"
//              << small_filename << "\t2\tcluster58\n"
//              << small2_filename << "\t2\tcluster84\n"
//              << seq2_filename << "\t1\tcluster71\n"
//              << seq1_filename << "\t1\tcluster5\n"
//              << small2_filename << "\t2\tcluster14\n"
//              << small_filename << "\t2\tcluster78\n"
//              << seq1_filename << "\t1\tcluster85\n"
//              << seq2_filename << "\t1\tcluster56\n"
//              << small2_filename << "\t2\tcluster34\n"
//              << seq1_filename << "\t1\tcluster95\n"
//              << small_filename << "\t2\tcluster28\n"
//              << small_filename << "\t2\tcluster3\n"
//              << small2_filename << "\t2\tcluster29\n";
//     }

//     char const * const argv[] = {"./chopper-layout",
//                                  "--tmax", "4",
//                                  "--input-prefix", io_prefix.get_path().c_str(),
//                                  "--output-filename", layout_file.get_path().c_str()};
//     int const argc = sizeof(argv) / sizeof(*argv);

//     seqan3::argument_parser layout_parser{"chopper-layout", argc, argv, seqan3::update_notifications::off};
//     chopper::layout::execute(layout_parser);

//     std::vector<std::string> const expected_components
//     {
//         {"##CONFIG:"},
//         {"##{"},
//         {"##    \"config\": {"},
//         {"##        \"version\": 1,"},
//         {"##        \"input_prefix\": \"" + io_prefix.get_path().string() + "\","},
//         {"##        \"count_filename\": {"},
//         {"##            \"value0\": \"" + io_prefix.get_path().string() + ".count\""},
//         {"##        },"},
//         {"##        \"sketch_directory\": {"},
//         {"##            \"value0\": \"" + io_prefix.get_path().string() + "_sketches\""},
//         {"##        },"},
//         {"##        \"output_filename\": {"},
//         {"##            \"value0\": \"" + layout_file.get_path().string() + "\""},
//         {"##        },"},
//         {"##        \"tmax\": 64,"},
//         {"##        \"num_hash_functions\": 2,"},
//         {"##        \"false_positive_rate\": 0.05,"},
//         {"##        \"alpha\": 1.2,"},
//         {"##        \"max_rearrangement_ratio\": 0.5,"},
//         {"##        \"threads\": 1,"},
//         {"##        \"estimate_union\": false,"},
//         {"##        \"rearrange_user_bins\": false,"},
//         {"##        \"determine_best_tmax\": false,"},
//         {"##        \"force_all_binnings\": false"},
//         {"##    }"},
//         {"##}"},
//         {"##ENDCONFIG"},
//         {"#HIGH_LEVEL_IBF max_bin_id:10"},
//         {"#MERGED_BIN_10 max_bin_id:0"},
//         {"#MERGED_BIN_11 max_bin_id:0"},
//         {"#MERGED_BIN_12 max_bin_id:0"},
//         {"#MERGED_BIN_13 max_bin_id:0"},
//         {"#MERGED_BIN_14 max_bin_id:0"},
//         {"#MERGED_BIN_15 max_bin_id:0"},
//         {"#MERGED_BIN_16 max_bin_id:0"},
//         {"#MERGED_BIN_17 max_bin_id:0"},
//         {"#MERGED_BIN_18 max_bin_id:0"},
//         {"#MERGED_BIN_19 max_bin_id:0"},
//         {"#MERGED_BIN_20 max_bin_id:0"},
//         {"#MERGED_BIN_21 max_bin_id:0"},
//         {"#MERGED_BIN_22 max_bin_id:0"},
//         {"#MERGED_BIN_23 max_bin_id:0"},
//         {"#MERGED_BIN_24 max_bin_id:0"},
//         {"#MERGED_BIN_25 max_bin_id:0"},
//         {"#FILES\tBIN_INDICES\tNUMBER_OF_BINS"},
//         {seq2_filename + "\t0\t1"},
//         {seq3_filename + "\t1\t1"},
//         {seq1_filename + "\t2\t1"},
//         {seq1_filename + "\t3\t1"},
//         {seq3_filename + "\t4\t1"},
//         {seq2_filename + "\t5\t1"},
//         {seq2_filename + "\t6\t1"},
//         {seq1_filename + "\t7\t1"},
//         {seq2_filename + "\t8\t1"},
//         {seq3_filename + "\t9\t1"},
//         {seq2_filename + "\t10;0\t1;60"},
//         {seq3_filename + "\t10;60\t1;2"},
//         {seq3_filename + "\t10;62\t1;2"},
//         {seq1_filename + "\t11;0\t1;60"},
//         {seq1_filename + "\t11;60\t1;2"},
//         {seq2_filename + "\t11;62\t1;2"},
//         {seq1_filename + "\t12;0\t1;60"},
//         {seq2_filename + "\t12;60\t1;2"},
//         {seq3_filename + "\t12;62\t1;2"},
//         {seq3_filename + "\t13;0\t1;60"},
//         {seq1_filename + "\t13;60\t1;2"},
//         {seq2_filename + "\t13;62\t1;2"},
//         {seq1_filename + "\t14;0\t1;60"},
//         {seq2_filename + "\t14;60\t1;2"},
//         {seq2_filename + "\t14;62\t1;2"},
//         {seq3_filename + "\t15;0\t1;60"},
//         {seq1_filename + "\t15;60\t1;2"},
//         {seq1_filename + "\t15;62\t1;2"},
//         {seq2_filename + "\t16;0\t1;60"},
//         {seq3_filename + "\t16;60\t1;2"},
//         {seq3_filename + "\t16;62\t1;2"},
//         {seq1_filename + "\t17;0\t1;60"},
//         {seq1_filename + "\t17;60\t1;2"},
//         {seq2_filename + "\t17;62\t1;2"},
//         {seq2_filename + "\t18;0\t1;60"},
//         {seq1_filename + "\t18;60\t1;2"},
//         {seq2_filename + "\t18;62\t1;2"},
//         {seq3_filename + "\t19;0\t1;60"},
//         {seq2_filename + "\t19;60\t1;2"},
//         {seq1_filename + "\t19;62\t1;2"},
//         {seq2_filename + "\t20;0\t1;60"},
//         {seq3_filename + "\t20;60\t1;2"},
//         {seq3_filename + "\t20;62\t1;2"},
//         {seq3_filename + "\t21;0\t1;60"},
//         {seq1_filename + "\t21;60\t1;2"},
//         {seq2_filename + "\t21;62\t1;2"},
//         {seq1_filename + "\t22;0\t1;60"},
//         {seq3_filename + "\t22;60\t1;2"},
//         {seq3_filename + "\t22;62\t1;2"},
//         {seq3_filename + "\t23;0\t1;60"},
//         {seq1_filename + "\t23;60\t1;2"},
//         {seq1_filename + "\t23;62\t1;2"},
//         {seq3_filename + "\t24;0\t1;60"},
//         {seq1_filename + "\t24;60\t1;2"},
//         {seq1_filename + "\t24;62\t1;2"},
//         {seq2_filename + "\t25;0\t1;60"},
//         {seq3_filename + "\t25;60\t1;2"},
//         {seq2_filename + "\t25;62\t1;2"},
//         {small2_filename + "\t26\t1"},
//         {small2_filename + "\t27\t1"},
//         {small2_filename + "\t28\t1"},
//         {small_filename + "\t29\t1"},
//         {small_filename + "\t30\t1"},
//         {small2_filename + "\t31\t1"},
//         {small2_filename + "\t32\t1"},
//         {small_filename + "\t33\t1"},
//         {small_filename + "\t34\t1"},
//         {small_filename + "\t35\t1"},
//         {small_filename + "\t36\t1"},
//         {small_filename + "\t37\t1"},
//         {small_filename + "\t38\t1"},
//         {small2_filename + "\t39\t1"},
//         {small2_filename + "\t40\t1"},
//         {small2_filename + "\t41\t1"},
//         {small2_filename + "\t42\t1"},
//         {small2_filename + "\t43\t1"},
//         {small2_filename + "\t44\t1"},
//         {small_filename + "\t45\t1"},
//         {small2_filename + "\t46\t1"},
//         {small_filename + "\t47\t1"},
//         {small2_filename + "\t48\t1"},
//         {small_filename + "\t49\t1"},
//         {small_filename + "\t50\t1"},
//         {small_filename + "\t51\t1"},
//         {small_filename + "\t52\t1"},
//         {small2_filename + "\t53\t1"},
//         {small2_filename + "\t54\t1"},
//         {small_filename + "\t55\t1"},
//         {small_filename + "\t56\t1"},
//         {small2_filename + "\t57\t1"},
//         {small_filename + "\t58\t1"},
//         {small_filename + "\t59\t1"},
//         {small2_filename + "\t60\t1"},
//         {small_filename + "\t61\t1"},
//         {small2_filename + "\t62\t1"},
//         {small2_filename + "\t63\t1"}
//     };
//     std::string const actual_file{string_from_file(layout_file.get_path())};

//     size_t line_count{};
//     for (auto && line : actual_file | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
//     {
// #if defined(__GNUC__) && (__GNUC__ == 12)
//         if (line.empty())
//             continue;
// #endif
//         EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "Did not find " + line;
//         ++line_count;
//     }

//     EXPECT_EQ(expected_components.size(), line_count);
// }
