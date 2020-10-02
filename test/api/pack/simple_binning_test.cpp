#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sstream>
#include <vector>

#include <chopper/pack/simple_binning.hpp>

TEST(simple_binning_test, small_example)
{
    std::vector<size_t> const input{100, 40, 20, 20};
    std::vector<std::string> const names{"seq1", "seq2", "seq3", "seq4"};
    std::ostringstream output{};

    simple_binning algo{input, names, "TEST_IBF", output, 9};
    algo.dp_algorithm();

    std::string expected
    {
        "TEST_IBF_0\tseq4\t1\t20\n"
        "TEST_IBF_1\tseq3\t1\t20\n"
        "TEST_IBF_2\tseq2\t2\t20\n"
        "TEST_IBF_3\tseq1\t5\t20\n"
    };

    EXPECT_EQ(output.str(), expected);
}

TEST(simple_binning_test, uniform_distribution)
{
    std::vector<size_t> const input{20, 20, 20, 20};
    std::vector<std::string> const names{"seq1", "seq2", "seq3", "seq4"};
    std::ostringstream output{};

    simple_binning algo{input, names, "TEST_IBF", output, 4};
    algo.dp_algorithm();

    std::string expected
    {
        "TEST_IBF_0\tseq4\t1\t20\n"
        "TEST_IBF_1\tseq3\t1\t20\n"
        "TEST_IBF_2\tseq2\t1\t20\n"
        "TEST_IBF_3\tseq1\t1\t20\n"
    };

    EXPECT_EQ(output.str(), expected);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    std::vector<size_t> const input{100, 40, 20, 20};
    std::vector<std::string> const names{"seq1", "seq2", "seq3", "seq4"};
    std::ostringstream output{};

    simple_binning algo{input, names, "TEST_IBF", output, 2};
    EXPECT_THROW(algo.dp_algorithm(), std::logic_error);
}
