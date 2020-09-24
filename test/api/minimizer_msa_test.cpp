#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

#include "sequence_input.hpp"
#include "minimizer_msa.hpp"

TEST(minimizer_msa_test, simple_example)
{
    // set up config
    segment_generation_config<int> seg_gen_config;
    seg_gen_config.output_graph_file = DATADIR"output_small_graph.dot";
    seg_gen_config.kmer_size = 15;
    seg_gen_config.window_size = 25;

    seqan::StringSet<seqan::String<minimizer>> seqs;
    seqan::StringSet<seqan::String<char>> ids;
    seqan::String<size_t> lengths;
    load_minimizer_sequences(seqs, ids, lengths, seg_gen_config, DATADIR"small.fa");
    ASSERT_EQ(seqan::length(seqs), 3); // sanity check
    ASSERT_EQ(seqan::length(ids), 3); // sanity check
    ASSERT_EQ(seqan::length(lengths), 3); // sanity check

    // run MSA
    minimizer_msa(seqs, ids, lengths, seg_gen_config); // writes graph.out as a result

    // compare results
    std::ifstream expected_file{DATADIR"small_graph.dot"};
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    std::ifstream output_file{seg_gen_config.output_graph_file};

    std::string expected_line;
    std::string output_line;

    while (std::getline(expected_file, expected_line) && std::getline(output_file, output_line))
        EXPECT_EQ(expected_line, output_line);

    EXPECT_FALSE(std::getline(expected_file, expected_line)); // both files are exhausted
    EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted
}
