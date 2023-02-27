#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include "../api_test.hpp"
#include "print_debug_file.hpp"
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>

TEST(execute_test, few_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    chopper::configuration config{};
    config.tmax = 64;
    config.output_filename = layout_file;
    config.disable_sketch_output = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer,
                              .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                              .kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500}};

    chopper::layout::execute(config, store);

    std::string const expected_file{"##CONFIG:\n"
                                    "##{\n"
                                    "##    \"config\": {\n"
                                    "##        \"version\": 2,\n"
                                    "##        \"data_file\": {\n"
                                    "##            \"value0\": \"\"\n"
                                    "##        },\n"
                                    "##        \"debug\": false,\n"
                                    "##        \"sketch_directory\": {\n"
                                    "##            \"value0\": \"\"\n"
                                    "##        },\n"
                                    "##        \"k\": 19,\n"
                                    "##        \"sketch_bits\": 12,\n"
                                    "##        \"disable_sketch_output\": true,\n"
                                    "##        \"precomputed_files\": false,\n"
                                    "##        \"output_filename\": {\n"
                                    "##            \"value0\": \""
                                    + layout_file.string()
                                    + "\"\n"
                                      "##        },\n"
                                      "##        \"tmax\": 64,\n"
                                      "##        \"num_hash_functions\": 2,\n"
                                      "##        \"false_positive_rate\": 0.05,\n"
                                      "##        \"alpha\": 1.2,\n"
                                      "##        \"max_rearrangement_ratio\": 0.5,\n"
                                      "##        \"threads\": 1,\n"
                                      "##        \"estimate_union\": false,\n"
                                      "##        \"rearrange_user_bins\": false,\n"
                                      "##        \"determine_best_tmax\": false,\n"
                                      "##        \"force_all_binnings\": false\n"
                                      "##    }\n"
                                      "##}\n"
                                      "##ENDCONFIG\n"
                                      "#HIGH_LEVEL_IBF max_bin_id:0\n"
                                      "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\n"
                                      "seq7\t0\t6\n"
                                      "seq6\t6\t6\n"
                                      "seq5\t12\t6\n"
                                      "seq4\t18\t6\n"
                                      "seq3\t24\t6\n"
                                      "seq2\t30\t6\n"
                                      "seq0\t36\t6\n"
                                      "seq1\t42\t22\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file) << actual_file;
}

TEST(execute_test, set_default_tmax)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    chopper::configuration config{}; // tmax == 0 triggers to set default to the sqrt(#samples)
    config.output_filename = layout_file;
    config.disable_sketch_output = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer,
                              .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                              .kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500}};

    chopper::layout::execute(config, store);

    EXPECT_EQ(config.tmax, 64u);
}

TEST(execute_test, few_ubs_debug)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    chopper::configuration config{};
    config.tmax = 64;
    config.output_filename = layout_file;
    config.debug = true;
    config.disable_sketch_output = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer,
                              .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                              .kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500}};

    chopper::layout::execute(config, store);

    std::string const expected_file{

        "##CONFIG:\n"
        "##{\n"
        "##    \"config\": {\n"
        "##        \"version\": 2,\n"
        "##        \"data_file\": {\n"
        "##            \"value0\": \"\"\n"
        "##        },\n"
        "##        \"debug\": true,\n"
        "##        \"sketch_directory\": {\n"
        "##            \"value0\": \"\"\n"
        "##        },\n"
        "##        \"k\": 19,\n"
        "##        \"sketch_bits\": 12,\n"
        "##        \"disable_sketch_output\": true,\n"
        "##        \"precomputed_files\": false,\n"
        "##        \"output_filename\": {\n"
        "##            \"value0\": \""
        + layout_file.string()
        + "\"\n"
          "##        },\n"
          "##        \"tmax\": 64,\n"
          "##        \"num_hash_functions\": 2,\n"
          "##        \"false_positive_rate\": 0.05,\n"
          "##        \"alpha\": 1.2,\n"
          "##        \"max_rearrangement_ratio\": 0.5,\n"
          "##        \"threads\": 1,\n"
          "##        \"estimate_union\": false,\n"
          "##        \"rearrange_user_bins\": false,\n"
          "##        \"determine_best_tmax\": false,\n"
          "##        \"force_all_binnings\": false\n"
          "##    }\n"
          "##}\n"
          "##ENDCONFIG\n"
          "#HIGH_LEVEL_IBF max_bin_id:0\n"
          "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX\n"
          "seq7\t0\t6\t84\t232\t2.61\t64\n"
          "seq6\t6\t6\t84\t232\t2.61\t64\n"
          "seq5\t12\t6\t84\t232\t2.61\t64\n"
          "seq4\t18\t6\t84\t232\t2.61\t64\n"
          "seq3\t24\t6\t84\t232\t2.61\t64\n"
          "seq2\t30\t6\t84\t232\t2.61\t64\n"
          "seq0\t36\t6\t84\t232\t2.61\t64\n"
          "seq1\t42\t22\t46\t232\t5.12\t64\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file);
    // print_debug_file(layout_file); // Formatted output
}

TEST(execute_test, many_ubs_debug)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};

    std::vector<std::string> many_filenames;
    std::vector<size_t> many_kmer_counts;

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    for (size_t i{0}; i < 96u; ++i)
    {
        many_filenames.push_back(seqan3::detail::to_string("seq", i));
        many_kmer_counts.push_back(100 * ((i + 20) / 20));
    }

    chopper::configuration config{};
    config.tmax = 64;
    config.output_filename = layout_file;
    config.debug = true;
    config.disable_sketch_output = true;

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .output_buffer = &output_buffer,
                             .header_buffer = &header_buffer,
                             .filenames = many_filenames,
                             .kmer_counts = many_kmer_counts};

    chopper::layout::execute(config, data);

    std::string const expected_file{"##CONFIG:\n"
                                    "##{\n"
                                    "##    \"config\": {\n"
                                    "##        \"version\": 2,\n"
                                    "##        \"data_file\": {\n"
                                    "##            \"value0\": \"\"\n"
                                    "##        },\n"
                                    "##        \"debug\": true,\n"
                                    "##        \"sketch_directory\": {\n"
                                    "##            \"value0\": \"\"\n"
                                    "##        },\n"
                                    "##        \"k\": 19,\n"
                                    "##        \"sketch_bits\": 12,\n"
                                    "##        \"disable_sketch_output\": true,\n"
                                    "##        \"precomputed_files\": false,\n"
                                    "##        \"output_filename\": {\n"
                                    "##            \"value0\": \""
                                    + layout_file.string()
                                    + "\"\n"
                                      "##        },\n"
                                      "##        \"tmax\": 64,\n"
                                      "##        \"num_hash_functions\": 2,\n"
                                      "##        \"false_positive_rate\": 0.05,\n"
                                      "##        \"alpha\": 1.2,\n"
                                      "##        \"max_rearrangement_ratio\": 0.5,\n"
                                      "##        \"threads\": 1,\n"
                                      "##        \"estimate_union\": false,\n"
                                      "##        \"rearrange_user_bins\": false,\n"
                                      "##        \"determine_best_tmax\": false,\n"
                                      "##        \"force_all_binnings\": false\n"
                                      "##    }\n"
                                      "##}\n"
                                      "##ENDCONFIG\n"
                                      "#HIGH_LEVEL_IBF max_bin_id:0\n"
                                      "#MERGED_BIN_0 max_bin_id:14\n"
                                      "#MERGED_BIN_1 max_bin_id:14\n"
                                      "#MERGED_BIN_2 max_bin_id:14\n"
                                      "#MERGED_BIN_3 max_bin_id:52\n"
                                      "#MERGED_BIN_4 max_bin_id:0\n"
                                      "#MERGED_BIN_5 max_bin_id:0\n"
                                      "#MERGED_BIN_6 max_bin_id:0\n"
                                      "#MERGED_BIN_7 max_bin_id:0\n"
                                      "#MERGED_BIN_8 max_bin_id:0\n"
                                      "#MERGED_BIN_9 max_bin_id:0\n"
                                      "#MERGED_BIN_26 max_bin_id:0\n"
                                      "#MERGED_BIN_27 max_bin_id:0\n"
                                      "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX\n"
                                      "seq15\t0;0\t1;14\t600;8\t600;34\t1.00;4.06\t64;64\n"
                                      "seq16\t0;14\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq17\t0;24\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq18\t0;34\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq19\t0;44\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq0\t0;54\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq9\t1;0\t1;14\t600;8\t600;34\t1.00;4.06\t64;64\n"
                                      "seq10\t1;14\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq11\t1;24\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq12\t1;34\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq13\t1;44\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq14\t1;54\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq3\t2;0\t1;14\t600;8\t600;34\t1.00;4.06\t64;64\n"
                                      "seq4\t2;14\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq5\t2;24\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq6\t2;34\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq7\t2;44\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq8\t2;54\t1;10\t600;10\t600;34\t1.00;3.41\t64;64\n"
                                      "seq33\t3;0\t1;27\t600;8\t600;43\t1.00;5.68\t64;64\n"
                                      "seq32\t3;27\t1;25\t600;8\t600;43\t1.00;5.46\t64;64\n"
                                      "seq1\t3;52\t1;6\t600;17\t600;43\t1.00;2.61\t64;64\n"
                                      "seq2\t3;58\t1;6\t600;17\t600;43\t1.00;2.61\t64;64\n"
                                      "seq30\t4;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq35\t4;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq34\t4;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq38\t5;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq37\t5;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq36\t5;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq20\t6;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq31\t6;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq39\t6;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq23\t7;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq22\t7;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq21\t7;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq27\t8;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq26\t8;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq25\t8;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq29\t9;0\t1;22\t600;10\t600;47\t1.00;5.12\t64;64\n"
                                      "seq24\t9;22\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq28\t9;43\t1;21\t600;10\t600;47\t1.00;5.00\t64;64\n"
                                      "seq40\t10\t1\t300\t600\t1.00\t64\n"
                                      "seq41\t11\t1\t300\t600\t1.00\t64\n"
                                      "seq42\t12\t1\t300\t600\t1.00\t64\n"
                                      "seq43\t13\t1\t300\t600\t1.00\t64\n"
                                      "seq44\t14\t1\t300\t600\t1.00\t64\n"
                                      "seq45\t15\t1\t300\t600\t1.00\t64\n"
                                      "seq46\t16\t1\t300\t600\t1.00\t64\n"
                                      "seq47\t17\t1\t300\t600\t1.00\t64\n"
                                      "seq48\t18\t1\t300\t600\t1.00\t64\n"
                                      "seq59\t19\t1\t300\t600\t1.00\t64\n"
                                      "seq58\t20\t1\t300\t600\t1.00\t64\n"
                                      "seq57\t21\t1\t300\t600\t1.00\t64\n"
                                      "seq56\t22\t1\t300\t600\t1.00\t64\n"
                                      "seq55\t23\t1\t300\t600\t1.00\t64\n"
                                      "seq54\t24\t1\t300\t600\t1.00\t64\n"
                                      "seq53\t25\t1\t300\t600\t1.00\t64\n"
                                      "seq51\t26;0\t1;32\t600;10\t600;58\t1.00;6.20\t64;64\n"
                                      "seq52\t26;32\t1;32\t600;10\t600;58\t1.00;6.20\t64;64\n"
                                      "seq49\t27;0\t1;32\t600;10\t600;58\t1.00;6.20\t64;64\n"
                                      "seq50\t27;32\t1;32\t600;10\t600;58\t1.00;6.20\t64;64\n"
                                      "seq79\t28\t1\t400\t600\t1.00\t64\n"
                                      "seq78\t29\t1\t400\t600\t1.00\t64\n"
                                      "seq77\t30\t1\t400\t600\t1.00\t64\n"
                                      "seq76\t31\t1\t400\t600\t1.00\t64\n"
                                      "seq75\t32\t1\t400\t600\t1.00\t64\n"
                                      "seq74\t33\t1\t400\t600\t1.00\t64\n"
                                      "seq73\t34\t1\t400\t600\t1.00\t64\n"
                                      "seq70\t35\t1\t400\t600\t1.00\t64\n"
                                      "seq71\t36\t1\t400\t600\t1.00\t64\n"
                                      "seq72\t37\t1\t400\t600\t1.00\t64\n"
                                      "seq60\t38\t1\t400\t600\t1.00\t64\n"
                                      "seq61\t39\t1\t400\t600\t1.00\t64\n"
                                      "seq62\t40\t1\t400\t600\t1.00\t64\n"
                                      "seq63\t41\t1\t400\t600\t1.00\t64\n"
                                      "seq64\t42\t1\t400\t600\t1.00\t64\n"
                                      "seq65\t43\t1\t400\t600\t1.00\t64\n"
                                      "seq66\t44\t1\t400\t600\t1.00\t64\n"
                                      "seq67\t45\t1\t400\t600\t1.00\t64\n"
                                      "seq68\t46\t1\t400\t600\t1.00\t64\n"
                                      "seq69\t47\t1\t400\t600\t1.00\t64\n"
                                      "seq88\t48\t1\t500\t600\t1.00\t64\n"
                                      "seq94\t49\t1\t500\t600\t1.00\t64\n"
                                      "seq93\t50\t1\t500\t600\t1.00\t64\n"
                                      "seq92\t51\t1\t500\t600\t1.00\t64\n"
                                      "seq91\t52\t1\t500\t600\t1.00\t64\n"
                                      "seq90\t53\t1\t500\t600\t1.00\t64\n"
                                      "seq89\t54\t1\t500\t600\t1.00\t64\n"
                                      "seq87\t55\t1\t500\t600\t1.00\t64\n"
                                      "seq86\t56\t1\t500\t600\t1.00\t64\n"
                                      "seq85\t57\t1\t500\t600\t1.00\t64\n"
                                      "seq84\t58\t1\t500\t600\t1.00\t64\n"
                                      "seq83\t59\t1\t500\t600\t1.00\t64\n"
                                      "seq82\t60\t1\t500\t600\t1.00\t64\n"
                                      "seq81\t61\t1\t500\t600\t1.00\t64\n"
                                      "seq80\t62\t1\t500\t600\t1.00\t64\n"
                                      "seq95\t63\t1\t500\t600\t1.00\t64\n"};
    std::string const actual_file{string_from_file(layout_file)};

    EXPECT_EQ(actual_file, expected_file) << actual_file << std::endl;
}
