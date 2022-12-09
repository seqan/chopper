#include <gtest/gtest.h>

#include <filesystem>

#include "../api_test.hpp"
#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/filenames_data_input.hpp>

TEST(read_filename_data_file_test, file_open_error)
{
    chopper::configuration config;
    chopper::data_store data;
    config.count_filename = "non_existing.file";

    EXPECT_THROW(chopper::layout::read_filename_data_file(data, config), std::runtime_error);
}

TEST(read_filename_data_file_test, only_filenames)
{
    chopper::configuration config;
    config.count_filename = data("only_filenames.tsv");

    chopper::data_store data;
    EXPECT_THROW(chopper::layout::read_filename_data_file(data, config), std::runtime_error);
}

TEST(read_filename_data_file_test, filenames_and_counts)
{
    chopper::configuration config;
    config.count_filename = data("filenames_and_counts.tsv");

    chopper::data_store data;
    chopper::layout::read_filename_data_file(data, config);

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
    chopper::configuration config;
    config.count_filename = data("filenames_counts_and_extra_information.tsv");

    chopper::data_store data;
    chopper::layout::read_filename_data_file(data, config);

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
