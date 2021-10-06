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
    using user_bin_sequence::distance_matrix;
    using user_bin_sequence::prio_queue;
    using user_bin_sequence::clustering_node;
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

TEST_F(user_bin_sequence_test, random_shuffle)
{
    prio_queue default_pq{};
    distance_matrix dist{{0, default_pq},{1, default_pq},{2, default_pq},{3, default_pq},{4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> ids{{0,0},{1,1},{2,2},{3,3},{4,4}};

    this->random_shuffle(dist, ids);

    // since randomness is seeded, the output is deterministic
    auto [new_pos_0, new_pos_1, new_pos_2, new_pos_3, new_pos_4] = std::make_tuple(3u, 2u, 1u, 0u, 4u);

    EXPECT_EQ(ids[0], new_pos_0);
    EXPECT_EQ(ids[1], new_pos_1);
    EXPECT_EQ(ids[2], new_pos_2);
    EXPECT_EQ(ids[3], new_pos_3);
    EXPECT_EQ(ids[4], new_pos_4);

    EXPECT_EQ(dist[new_pos_0].id, 0u);
    EXPECT_EQ(dist[new_pos_1].id, 1u);
    EXPECT_EQ(dist[new_pos_2].id, 2u);
    EXPECT_EQ(dist[new_pos_3].id, 3u);
    EXPECT_EQ(dist[new_pos_4].id, 4u);
}

TEST_F(user_bin_sequence_test, prune)
{
    prio_queue default_pq{};
    distance_matrix dist{{0, default_pq},{1, default_pq},{2, default_pq},{3, default_pq},{4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> remaining_ids{{0,0},{1,1},{2,2},{3,3},{4,4}};

    // since remaining_ids contains all_ids, prune shouldn't do anything. All ids are valid.
    this->prune(dist, remaining_ids);

    EXPECT_EQ(remaining_ids[0], 0u);
    EXPECT_EQ(remaining_ids[1], 1u);
    EXPECT_EQ(remaining_ids[2], 2u);
    EXPECT_EQ(remaining_ids[3], 3u);
    EXPECT_EQ(remaining_ids[4], 4u);

    EXPECT_EQ(dist.size(), 5u);
    EXPECT_EQ(dist[0].id, 0u);
    EXPECT_EQ(dist[1].id, 1u);
    EXPECT_EQ(dist[2].id, 2u);
    EXPECT_EQ(dist[3].id, 3u);
    EXPECT_EQ(dist[4].id, 4u);

    remaining_ids.erase(1);
    remaining_ids.erase(3);

    // distance entry 1 and 3 are now invalid, since they do not occur in remaining_ids
    // prune() should therefore remove them from dist.
    this->prune(dist, remaining_ids);

    EXPECT_EQ(remaining_ids[0], 0u);
    EXPECT_EQ(remaining_ids[2], 2u);
    EXPECT_EQ(remaining_ids[4], 1u);

    EXPECT_EQ(dist.size(), 3u);
    EXPECT_EQ(dist[0].id, 0u);
    EXPECT_EQ(dist[1].id, 4u);
    EXPECT_EQ(dist[2].id, 2u);
}

TEST_F(user_bin_sequence_test, rotate)
{
    hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (0,1)     (2,3)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<clustering_node> clustering{{f,f,s}, {f,f,s}, {f,f,s}, {f,f,s}, // the leaves come first
                                            {5,6,s}, {0,1,s}, {2,3,s}};

    // previous_rightmost is already at the very left. Nothing has to be rotated.
    rotate(clustering, 0/*previous_rightmost*/, 0/*interval_start*/, 4/*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(5u, 6u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));

    // now the previous_rightmost is within the tree. Rotation should take place
    rotate(clustering, 2/*previous_rightmost*/, 0/*interval_start*/, 4/*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(6u, 5u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));
}

TEST_F(user_bin_sequence_test, trace)
{
    hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (1,3)     (2,0)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<clustering_node> clustering{{f,f,s}, {f,f,s}, {f,f,s}, {f,f,s}, // the leaves come first
                                            {5,6,s}, {1,3,s}, {2,0,s}};

    std::vector<size_t> permutation{};

    this->trace(clustering, permutation, 2/*previous_rightmost*/, 0/*interval_start*/, 4/*root_id*/);

    EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{1, 3, 0}));
}
