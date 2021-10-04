#include <gtest/gtest.h>

#include <chopper/union/user_bin_sequence.hpp>

#include "../api_test.hpp"

// inherits from user_bin_sequence to test private members
struct user_bin_sequence_test : public ::testing::Test, public user_bin_sequence
{
public:

    std::vector<std::string> global_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> global_kmer_counts{500, 600, 700, 800};

    std::vector<std::string> const & get_filenames()
    {
        return *this->filenames;
    }

    std::vector<size_t> const & get_user_bin_kmer_counts()
    {
        return *this->user_bin_kmer_counts;
    }

    user_bin_sequence_test() :
        user_bin_sequence{global_filenames, global_kmer_counts}
    {}

    using user_bin_sequence::apply_permutation;
};

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

TEST_F(user_bin_sequence_test, apply_permutation)
{
    std::vector<size_t> const permutation{2, 1, 3, 0};

    EXPECT_RANGE_EQ(this->get_filenames(), global_filenames);
    EXPECT_RANGE_EQ(this->get_user_bin_kmer_counts(), global_kmer_counts);

    this->apply_permutation(permutation);

    EXPECT_RANGE_EQ(this->get_filenames(), (std::vector<std::string>{"small2.fa", "small.fa", "small2.fa", "small.fa"}));
    EXPECT_RANGE_EQ(this->get_user_bin_kmer_counts(), (std::vector<size_t>{700, 600, 800, 500}));
}

TEST_F(user_bin_sequence_test, sort_by_cardinalities)
{
    this->sort_by_cardinalities();

    EXPECT_RANGE_EQ(this->get_filenames(), (std::vector<std::string>{"small2.fa", "small2.fa", "small.fa", "small.fa"}));
    EXPECT_RANGE_EQ(this->get_user_bin_kmer_counts(), (std::vector<size_t>{800, 700, 600, 500}));
}
