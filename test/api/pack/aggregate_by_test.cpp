#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <chopper/pack/aggregate_by.hpp>
#include <chopper/pack/pack_data.hpp>

#include "../api_test.hpp"

TEST(sort_by_test, small_example)
{
    pack_data data;
    data.filenames = std::vector<std::string>{"seq1", "seq2", "seq3", "seq4", "seq5"};
    data.kmer_counts = std::vector<size_t>{100, 40, 20, 20, 5};
    data.extra_information = std::vector<std::vector<std::string>>{{"1", "foo"},
                                                                   {"2", "foo"},
                                                                   {"1", "moo"},
                                                                   {"2", "moo"},
                                                                   {"3", "moo"}};
    {
        pack_data copy = data;
        sort_by(copy, 0);

        std::vector<std::string> expected_filenames{"seq1", "seq3", "seq2", "seq4", "seq5"};
        std::vector<size_t> expected_kmer_counts{100, 20, 40, 20, 5};
        std::vector<std::vector<std::string>> expected_extra_information{{"1", "foo"},
                                                                         {"1", "moo"},
                                                                         {"2", "foo"},
                                                                         {"2", "moo"},
                                                                         {"3", "moo"}};

        EXPECT_RANGE_EQ(copy.filenames, expected_filenames);
        EXPECT_RANGE_EQ(copy.kmer_counts, expected_kmer_counts);
        EXPECT_RANGE_EQ(copy.extra_information, expected_extra_information);
    }
}

TEST(aggregate_by_test, small_example)
{
    pack_data data;
    data.filenames = std::vector<std::string>{"seq1", "seq2", "seq3", "seq4", "seq5"};
    data.kmer_counts = std::vector<size_t>{100, 40, 20, 20, 5};
    data.extra_information = std::vector<std::vector<std::string>>{{"1", "foo"},
                                                                   {"2", "foo"},
                                                                   {"1", "moo"},
                                                                   {"2", "moo"},
                                                                   {"3", "moo"}};
    {
        pack_data copy = data;
        aggregate_by(copy, 0);

        std::vector<std::string> expected_filenames{"seq1;seq3", "seq2;seq4", "seq5"};
        std::vector<size_t> expected_kmer_counts{120, 60, 5};
        std::vector<std::vector<std::string>> expected_extra_information{{"1", "foo"}, {"2", "foo"}, {"3", "moo"}};

        EXPECT_RANGE_EQ(copy.filenames, expected_filenames);
        EXPECT_RANGE_EQ(copy.kmer_counts, expected_kmer_counts);
        EXPECT_RANGE_EQ(copy.extra_information, expected_extra_information);
    }

    {
        pack_data copy = data;
        aggregate_by(copy, 1);

        std::vector<std::string> expected_filenames{"seq1;seq2", "seq3;seq4;seq5"};
        std::vector<size_t> expected_kmer_counts{140, 45};
        std::vector<std::vector<std::string>> expected_extra_information{{"1", "foo"}, {"1", "moo"}};

        EXPECT_RANGE_EQ(copy.filenames, expected_filenames);
        EXPECT_RANGE_EQ(copy.kmer_counts, expected_kmer_counts);
        EXPECT_RANGE_EQ(copy.extra_information, expected_extra_information);
    }
}

#ifndef NDEBUG
TEST(aggregate_by_test, column_index_to_sort_too_big)
{
    pack_data data;
    data.filenames = std::vector<std::string>{"seq1", "seq2"};
    data.kmer_counts = std::vector<size_t>{100, 40};
    data.extra_information = std::vector<std::vector<std::string>>{{"1", "foo"}, {"2", "foo"}};

    EXPECT_DEATH((aggregate_by(data, 4)), "");
}
#endif
