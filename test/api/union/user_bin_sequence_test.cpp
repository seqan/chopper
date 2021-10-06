#include <gtest/gtest.h>

#include <chopper/union/user_bin_sequence.hpp>

#include "../api_test.hpp"

#include <seqan3/io/sequence_file/input.hpp>
struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};
using sequence_file_type = seqan3::sequence_file_input<input_traits, seqan3::fields<seqan3::field::seq>>;

// inherits from user_bin_sequence to test private members
struct user_bin_sequence_test : public ::testing::Test, public user_bin_sequence
{
public:

    std::vector<std::string> test_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> test_kmer_counts{500, 600, 700, 800};

    std::vector<hyperloglog> const & get_sketches()
    {
        return this->sketches;
    }

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

TEST_F(user_bin_sequence_test, read_hll_files)
{
    size_t const k{16};
    size_t const b{4};
    hyperloglog expected{b};

    std::string const input_file{data("small.fa")};
    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the sketch
    for (auto && [seq] : seq_file)
    {
        // we have to go C-style here for the HyperLogLog Interface
        const char * it = &(*seq.begin());
        char const * const end = it + seq.size() - k + 1;

        for (; it != end; ++it)
            expected.add(it, k);
    }

    this->read_hll_files(data(""));

    EXPECT_EQ(this->get_sketches().size(), 4u);
    for (auto const & sketch : this->get_sketches())
        EXPECT_EQ(sketch.estimate(), expected.estimate());
}

TEST_F(user_bin_sequence_test, read_hll_files_empty_dir)
{
    seqan3::test::tmp_filename const tmp_file{"some_file"};
    auto const parent_dir = tmp_file.get_path().parent_path();
    auto const empty_dir = parent_dir / "empty_dir";
    std::filesystem::create_directory(empty_dir); // create empty dir

    EXPECT_THROW(this->read_hll_files(empty_dir), std::runtime_error);
}

TEST_F(user_bin_sequence_test, read_hll_files_faulty_file)
{
    seqan3::test::tmp_filename const tmp_file{"small.hll"};
    {
        std::ofstream os{tmp_file.get_path()};
        os << "I am not what an hll file looks like\n";
    }

    EXPECT_THROW(this->read_hll_files(tmp_file.get_path().parent_path()), std::runtime_error);
}
