#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include "../api_test.hpp"

TEST(read_data_file_test, file_open_error)
{
    chopper::configuration config{};
    chopper::data_store store{};
    config.data_file = data("non_existing.file");
    EXPECT_THROW(chopper::sketch::read_data_file(config, store), std::runtime_error);
}

TEST(read_data_file_test, small_example)
{
    chopper::configuration config;
    chopper::data_store store{};
    config.data_file = data("seqinfo.tsv");

    chopper::sketch::read_data_file(config, store);

    std::vector<std::string> filenames{"file1", "file2", "file3", "file4", "file5"};
    std::vector<std::string> extra_information_strings{"1	foo", "2	foo", "1	moo", "2	moo", "3	moo"};
    EXPECT_RANGE_EQ(store.filenames, filenames);
    EXPECT_RANGE_EQ(store.extra_information_strings, extra_information_strings);
}
