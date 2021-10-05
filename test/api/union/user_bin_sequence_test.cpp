#include <gtest/gtest.h>

#include <chopper/union/user_bin_sequence.hpp>

#include "../api_test.hpp"

// inherits from user_bin_sequence to test private members
struct user_bin_sequence_test : public ::testing::Test, public user_bin_sequence
{
public:

    std::vector<std::string> test_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> test_kmer_counts{500, 600, 700, 800};

    user_bin_sequence_test() :
        user_bin_sequence{test_filenames, test_kmer_counts}
    {}

    using user_bin_sequence::apply_permutation;
};

TEST_F(user_bin_sequence_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_copy_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_move_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_destructible<user_bin_sequence>::value);

    // class has a const pointer
    EXPECT_FALSE(std::is_copy_assignable<user_bin_sequence>::value);
    EXPECT_FALSE(std::is_move_assignable<user_bin_sequence>::value);

    // construction from filenames and kmer_counts
    std::vector<std::string> filenames{"small.fa", "small.fa"};
    std::vector<size_t> kmer_counts{500, 500};
    user_bin_sequence ubs{filenames, kmer_counts};
}

TEST_F(user_bin_sequence_test, apply_permutation)
{
    std::vector<size_t> const permutation{2, 1, 3, 0};

    this->apply_permutation(permutation);

    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small2.fa", "small.fa", "small2.fa", "small.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{700, 600, 800, 500}));
}

TEST_F(user_bin_sequence_test, sort_by_cardinalities)
{
    this->sort_by_cardinalities();

    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small2.fa", "small2.fa", "small.fa", "small.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{800, 700, 600, 500}));
}
