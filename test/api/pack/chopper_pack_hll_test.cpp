#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/count/count_config.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/pack/chopper_pack.hpp>

#include "../api_test.hpp"
#include "print_debug_file.hpp"

TEST(chopper_pack_hll_test, small_example_hll)
{
    seqan3::test::tmp_filename count_file{"kmer_counts.txt"};
    seqan3::test::tmp_filename hll_dir{"hll"};

    // COUNT
    {
        count_config config;
        config.k = 15;
        config.w = 25;
        config.num_threads = 1;
        config.exclusively_hlls = true;
        config.hll_dir = hll_dir.get_path();
        config.output_filename = count_file.get_path();

        robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters;
        // std::string expected{}; // Needs check for components
        std::string short_names{};

        std::vector<std::pair<std::string, std::string>> files{{"seq2.fa", "77"},
                                                               {"small.fa", "96"},
                                                               {"small2.fa", "96"},
                                                               {"seq3.fa", "76"},
                                                               {"seq1.fa", "65"}};

        for (auto && [filename, size] : files)
        {
            std::string input_file = data(filename.c_str());

            filename_clusters[input_file].push_back(input_file);

            // expected += input_file + '\t' + size + '\t' + input_file + '\n';
            short_names += filename + '\t' + "10" + '\n'; // size=10 to actually make a difference for estimate_union
        }

        count_kmers(filename_clusters, config);

        // {
        //     std::ifstream output_file{count_file.get_path()};
        //     std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
        //     EXPECT_EQ(expected, output_file_str);
        // }

        {
            std::ofstream fout{count_file.get_path()};
            fout << short_names;
        }
    }

    // PACK
    seqan3::test::tmp_filename output_filename{"output.binning"};
    const char * argv[] = {"./chopper-pack", "-b", "4", "--debug", "-d", hll_dir.get_path().c_str(), "-u",
                           "-f", count_file.get_path().c_str(), "-o", output_filename.get_path().c_str()};
    int argc = 11;
    seqan3::argument_parser pack_parser{"chopper-pack", argc, argv, seqan3::update_notifications::off};

    chopper_pack(pack_parser);

    std::string expected_file
    {
        "#HIGH_LEVEL_IBF max_bin_id:0\n"
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\tSCORE\tCORR\tT_MAX\n"
        "seq3.fa\t0;0\t1;39\t20;1\t80;2\t1.00;10.50\t4;64\n"
        "seq1.fa\t0;39\t1;25\t20;1\t80;2\t1.00;7.50\t4;64\n"
        "small2.fa\t1\t1\t10\t80\t1.00\t4\n"
        "small.fa\t2\t1\t10\t80\t1.00\t4\n"
        "seq2.fa\t3\t1\t10\t80\t1.00\t4\n"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file);
    // print_debug_file(output_filename.get_path());
}
