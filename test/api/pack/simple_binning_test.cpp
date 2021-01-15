#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sstream>
#include <vector>

#include <chopper/pack/simple_binning.hpp>

TEST(simple_binning_test, small_example)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.kmer_counts = {100, 40, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};

    simple_binning algo{data, "TEST_IBF",  9};
    size_t max_bin = algo.execute();

    std::string expected
    {
        "TEST_IBF_0\tseq4\t1\t20\n"
        "TEST_IBF_1\tseq3\t1\t20\n"
        "TEST_IBF_2\tseq2\t2\t20\n"
        "TEST_IBF_4\tseq1\t5\t20\n"
    };

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0);
}

TEST(simple_binning_test, uniform_distribution)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.kmer_counts = {20, 20, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};

    simple_binning algo{data, "TEST_IBF",  4};
    size_t max_bin = algo.execute();

    std::string expected
    {
        "TEST_IBF_0\tseq4\t1\t20\n"
        "TEST_IBF_1\tseq3\t1\t20\n"
        "TEST_IBF_2\tseq2\t1\t20\n"
        "TEST_IBF_3\tseq1\t1\t20\n"
    };

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.kmer_counts = {100, 40, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};

    EXPECT_THROW((simple_binning{data, "TEST_IBF",  2}), std::runtime_error);
}
