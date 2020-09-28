#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

#include "split/split_data.hpp"
#include "split/sequence_input.hpp"
#include "split/minimizer_msa.hpp"

TEST(minimizer_msa_test, simple_example)
{
    // set up config
    split_config config;
    config.output_graph_file = DATADIR"output_small_graph.dot";
    config.kmer_size = 15;
    config.window_size = 25;

    split_data data;
    load_minimizer_sequences(data, config, DATADIR"small.fa");
    ASSERT_EQ(seqan::length(data.sequences), 3); // sanity check
    ASSERT_EQ(seqan::length(data.ids), 3); // sanity check
    ASSERT_EQ(seqan::length(data.lengths), 3); // sanity check

    // run MSA
    minimizer_msa(data, config); // writes graph.out as a result

    // compare results
    std::ifstream expected_file{DATADIR"small_graph.dot"};
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    std::ifstream output_file{config.output_graph_file};

    std::string expected_line;
    std::string output_line;

    while (std::getline(expected_file, expected_line) && std::getline(output_file, output_line))
        EXPECT_EQ(expected_line, output_line);

    EXPECT_FALSE(std::getline(expected_file, expected_line)); // both files are exhausted
    EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted
}
