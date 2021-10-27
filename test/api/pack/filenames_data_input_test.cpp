#include <gtest/gtest.h>

#include <seqan3/std/filesystem>

#include <chopper/pack/filenames_data_input.hpp>
#include <chopper/pack/configuration.hpp>
#include <chopper/pack/data_store.hpp>

#include "../api_test.hpp"

TEST(read_filename_data_file_test, file_open_error)
{
    chopper::pack::configuration config;
    chopper::pack::data_store data;
    config.data_file = "non_existing.file";

    EXPECT_THROW(read_filename_data_file(data, config), std::runtime_error);
}

TEST(read_filename_data_file_test, only_filenames)
{
    chopper::pack::configuration config;
    config.data_file = data("only_filenames.tsv");

    chopper::pack::data_store data;
    EXPECT_THROW(read_filename_data_file(data, config), std::runtime_error);
}

TEST(read_filename_data_file_test, filenames_and_counts)
{
    chopper::pack::configuration config;
    config.data_file = data("filenames_and_counts.tsv");

    chopper::pack::data_store data;
    read_filename_data_file(data, config);

    ASSERT_EQ(data.filenames.size(), 3u);
    EXPECT_EQ(data.filenames[0], "file1");
    EXPECT_EQ(data.filenames[1], "file2");
    EXPECT_EQ(data.filenames[2], "file3");

    ASSERT_EQ(data.kmer_counts.size(), 3u);
    EXPECT_EQ(data.kmer_counts[0], 1000u);
    EXPECT_EQ(data.kmer_counts[1], 2000u);
    EXPECT_EQ(data.kmer_counts[2], 3000u);

    EXPECT_EQ(data.extra_information.size(), 3u);
    for (auto const & info : data.extra_information)
        EXPECT_TRUE(info.empty());
}

TEST(read_filename_data_file_test, filenames_counts_and_extra_information)
{
    chopper::pack::configuration config;
    config.data_file = data("filenames_counts_and_extra_information.tsv");

    chopper::pack::data_store data;
    read_filename_data_file(data, config);

    ASSERT_EQ(data.filenames.size(), 3u);
    EXPECT_EQ(data.filenames[0], "file1");
    EXPECT_EQ(data.filenames[1], "file2");
    EXPECT_EQ(data.filenames[2], "file3");

    ASSERT_EQ(data.kmer_counts.size(), 3u);
    EXPECT_EQ(data.kmer_counts[0], 1000u);
    EXPECT_EQ(data.kmer_counts[1], 2000u);
    EXPECT_EQ(data.kmer_counts[2], 3000u);

    ASSERT_EQ(data.extra_information.size(), 3u);
    EXPECT_EQ(data.extra_information[0], (std::vector<std::string>{"info_a", "info_b"}));
    EXPECT_EQ(data.extra_information[1], (std::vector<std::string>{"123", "456"}));
    EXPECT_EQ(data.extra_information[2], (std::vector<std::string>{"hi ho", "mi mo"}));
}
