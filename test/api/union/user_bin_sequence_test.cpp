#include <gtest/gtest.h>

#include <unordered_set>

#include <chopper/union/user_bin_sequence.hpp>

#include "../api_test.hpp"

// inherits from user_bin_sequence to test private members
struct user_bin_sequence_test : public ::testing::Test, public user_bin_sequence
{};

TEST_F(user_bin_sequence_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_copy_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_move_constructible<user_bin_sequence>::value);
    EXPECT_TRUE(std::is_destructible<user_bin_sequence>::value);

    EXPECT_FALSE(std::is_copy_assignable<user_bin_sequence>::value); // class has a const pointer
    EXPECT_FALSE(std::is_move_assignable<user_bin_sequence>::value); // class has a const pointer

    // construction from filenames and kmer_counts
    std::vector<std::string> filenames{"small.fa", "small.fa"};
    std::vector<size_t> kmer_counts{500, 500};
    user_bin_sequence ubs{filenames, kmer_counts};
}
