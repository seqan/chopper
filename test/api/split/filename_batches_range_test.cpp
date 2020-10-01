#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/split/split_config.hpp>
#include <chopper/split/filename_batches_range.hpp>

TEST(filename_batches_range_test, high_level_data_file)
{
    split_config config;
    config.data_filename = DATADIR"high_level_ibf.binning";

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::high_level);

    auto it = r.begin();
    EXPECT_RANGE_EQ((*it).seqfiles, std::vector<std::string>{"seq7"});
    EXPECT_EQ((*it).bins, 2);

    ++it;

    EXPECT_RANGE_EQ((*it).seqfiles, (std::vector<std::string>{"seq1.1", "seq1.2"}));
    EXPECT_EQ((*it).bins, 2);

    ++it;

    EXPECT_TRUE(it == r.end());
}

TEST(filename_batches_range_test, low_level_data_file)
{
    split_config config;
    config.data_filename = DATADIR"low_level_ibfs.binning";

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::low_level);

    std::vector<std::vector<std::string>> expected_seqfiles_range
    {
        {"seq0"}, {"seq2"}, {"seq3.1", "seq3.2", "seq3.3"}, {"seq4"}, {"seq5"}
    };

    std::vector<int> expected_bins_range{16, 12, 12, 12, 12};

    auto it = r.begin();
    for (size_t i = 0; i < expected_seqfiles_range.size(); ++i, ++it)
    {
        EXPECT_RANGE_EQ((*it).seqfiles, expected_seqfiles_range[i]);
        EXPECT_EQ((*it).bins, expected_bins_range[i]);
    }
    EXPECT_TRUE(it == r.end());
}

TEST(filename_batches_range_test, no_data_file)
{
    split_config config;
    config.seqfiles = {"seq3.1", "seq3.2", "seq3.3"};

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::seqfiles_given);

    auto it = r.begin();

    EXPECT_RANGE_EQ((*it).seqfiles, config.seqfiles);
    EXPECT_EQ((*it).bins, config.bins);

    ++it;

    EXPECT_TRUE(it == r.end());
}
