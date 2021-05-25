#include <gtest/gtest.h>

#include <sstream>
#include <vector>
#include <chopper/count/count_config.hpp>
#include <chopper/count/read_data_file.hpp>

#include "../api_test.hpp"

TEST(read_data_file_test, small_example)
{
    count_config config;
    config.data_file = data("seqinfo.tsv");

    {
        config.column_index_to_cluster = 2;
        auto filename_clusters = read_data_file(config);

        EXPECT_RANGE_EQ(filename_clusters[std::string{"1"}], (std::vector<std::string>{"file1", "file3"}));
        EXPECT_RANGE_EQ(filename_clusters[std::string{"2"}], (std::vector<std::string>{"file2", "file4"}));
        EXPECT_RANGE_EQ(filename_clusters[std::string{"3"}], (std::vector<std::string>{"file5"}));
    }

    {
        config.column_index_to_cluster = 3;
        auto filename_clusters = read_data_file(config);

        EXPECT_RANGE_EQ(filename_clusters[std::string{"foo"}], (std::vector<std::string>{"file1", "file2"}));
        EXPECT_RANGE_EQ(filename_clusters[std::string{"moo"}], (std::vector<std::string>{"file3", "file4", "file5"}));
    }
}
