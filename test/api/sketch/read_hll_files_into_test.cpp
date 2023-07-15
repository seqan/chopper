#include <gtest/gtest.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <chopper/sketch/read_hll_files_into.hpp>

#include "../api_test.hpp"

struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};
using sequence_file_type = seqan3::sequence_file_input<input_traits, seqan3::fields<seqan3::field::seq>>;

// inherits from toolbox to test private members
struct read_hll_files_into_test : public ::testing::Test
{
    std::vector<std::string> test_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> test_kmer_counts{500, 600, 700, 800};
    std::vector<size_t> test_positions{0, 1, 2, 3};
    std::vector<hibf::sketch::hyperloglog> test_sketches = [this]()
    {
        std::vector<hibf::sketch::hyperloglog> result;
        chopper::sketch::read_hll_files_into(data(""), test_filenames, result);
        return result;
    }();
};

TEST_F(read_hll_files_into_test, basic)
{
    size_t const k{16};
    size_t const b{5};
    hibf::sketch::hyperloglog expected{b};

    std::vector<hibf::sketch::hyperloglog> target{};

    std::string const input_file{data("small.fa")};
    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the sketch
    for (auto && [seq] : seq_file)
    {
        // we have to go C-style here for the hibf::sketch::hyperloglog Interface
        char const * it = &(*seq.begin());
        char const * const end = it + seq.size() - k + 1;

        for (; it != end; ++it)
            expected.add(it, k);
    }

    chopper::sketch::read_hll_files_into(data(""), test_filenames, target);

    EXPECT_EQ(target.size(), 4u);
    for (auto const & sketch : target)
        EXPECT_EQ(sketch.estimate(), expected.estimate());
}

TEST_F(read_hll_files_into_test, empty_dir)
{
    seqan3::test::tmp_directory empty_dir{};
    ASSERT_TRUE(std::filesystem::exists(empty_dir.path()));
    ASSERT_TRUE(std::filesystem::is_empty(empty_dir.path()));

    std::vector<hibf::sketch::hyperloglog> target{};
    // Will throw in Release, but assert in Debug
#ifdef NDEBUG
    EXPECT_THROW(chopper::sketch::read_hll_files_into(empty_dir.path(), test_filenames, target), std::runtime_error);
#else
    EXPECT_DEATH(chopper::sketch::read_hll_files_into(empty_dir.path(), test_filenames, target), "");
#endif
}

TEST_F(read_hll_files_into_test, file_is_missing)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const tmp_file{tmp_dir.path() / "other.hll"};
    {
        std::ofstream os{tmp_file};
        os << "Doesn't matter I just need to exist\n";
    }

    std::vector<hibf::sketch::hyperloglog> target{};

    EXPECT_THROW(chopper::sketch::read_hll_files_into(tmp_file.parent_path(), test_filenames, target),
                 std::runtime_error);
}

TEST_F(read_hll_files_into_test, read_hll_files_into_faulty_file)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const tmp_file{tmp_dir.path() / "small.hll"};
    {
        std::ofstream os{tmp_file};
        os << "I am not what an hll file looks like\n";
    }

    std::vector<hibf::sketch::hyperloglog> target{};

    EXPECT_THROW(chopper::sketch::read_hll_files_into(tmp_file.parent_path(), test_filenames, target),
                 std::runtime_error);
}
