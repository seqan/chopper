#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include "../api_test.hpp"

TEST(read_data_file_test, file_open_error)
{
    chopper::configuration config{};
    std::vector<std::string> filenames{};
    config.data_file = data("non_existing.file");
    EXPECT_THROW(chopper::sketch::read_data_file(config, filenames), std::runtime_error);
}

TEST(read_data_file_test, small_example)
{
    chopper::configuration config;
    std::vector<std::string> filenames{};
    config.data_file = data("seqinfo.tsv");

    chopper::sketch::read_data_file(config, filenames);

    std::vector<std::string> expected_filenames{"file1", "file2", "file3", "file4", "file5"};
    EXPECT_RANGE_EQ(filenames, expected_filenames);
}
