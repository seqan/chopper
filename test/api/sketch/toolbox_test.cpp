#include <gtest/gtest.h>

#include <seqan3/io/sequence_file/input.hpp>

#include <chopper/sketch/toolbox.hpp>

#include "../api_test.hpp"
struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
};
using sequence_file_type = seqan3::sequence_file_input<input_traits, seqan3::fields<seqan3::field::seq>>;

// inherits from toolbox to test private members
struct toolbox_test : public ::testing::Test
{
    std::vector<std::string> test_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> test_kmer_counts{500, 600, 700, 800};
    std::vector<chopper::sketch::hyperloglog> test_sketches = [this]()
    {
        std::vector<chopper::sketch::hyperloglog> result;
        chopper::sketch::toolbox::read_hll_files_into(data(""), test_filenames, result);
        return result;
    }();
};

TEST_F(toolbox_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible<chopper::sketch::toolbox>::value);
    EXPECT_TRUE(std::is_copy_constructible<chopper::sketch::toolbox>::value);
    EXPECT_TRUE(std::is_move_constructible<chopper::sketch::toolbox>::value);
    EXPECT_TRUE(std::is_destructible<chopper::sketch::toolbox>::value);
    EXPECT_TRUE(std::is_copy_assignable<chopper::sketch::toolbox>::value);
    EXPECT_TRUE(std::is_move_assignable<chopper::sketch::toolbox>::value);

    // construction from filenames and kmer_counts and sketches
    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};
}

TEST_F(toolbox_test, apply_permutation)
{
    std::vector<size_t> const permutation{2, 1, 3, 0};

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    ubs.apply_permutation(permutation);

    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small2.fa", "small.fa", "small2.fa", "small.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{700, 600, 800, 500}));
}

TEST_F(toolbox_test, apply_permutation_with_sketches)
{
    std::vector<size_t> const permutation{2, 1, 3, 0};

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    ubs.apply_permutation(permutation);

    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small2.fa", "small.fa", "small2.fa", "small.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{700, 600, 800, 500}));
}

TEST_F(toolbox_test, sort_by_cardinalities)
{
    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    ubs.sort_by_cardinalities();

    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small2.fa", "small2.fa", "small.fa", "small.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{800, 700, 600, 500}));
}

TEST_F(toolbox_test, read_hll_files_into)
{
    size_t const k{16};
    size_t const b{5};
    chopper::sketch::hyperloglog expected{b};

    std::vector<chopper::sketch::hyperloglog> target{};

    std::string const input_file{data("small.fa")};
    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the sketch
    for (auto && [seq] : seq_file)
    {
        // we have to go C-style here for the chopper::sketch::HyperLogLog Interface
        char const * it = &(*seq.begin());
        char const * const end = it + seq.size() - k + 1;

        for (; it != end; ++it)
            expected.add(it, k);
    }

    chopper::sketch::toolbox::read_hll_files_into(data(""), test_filenames, target);

    EXPECT_EQ(target.size(), 4u);
    for (auto const & sketch : target)
        EXPECT_EQ(sketch.estimate(), expected.estimate());
}

TEST_F(toolbox_test, read_hll_files_into_empty_dir)
{
    seqan3::test::tmp_directory empty_dir{};
    ASSERT_TRUE(std::filesystem::exists(empty_dir.path()));
    ASSERT_TRUE(std::filesystem::is_empty(empty_dir.path()));

    std::vector<chopper::sketch::hyperloglog> target{};
    // Will throw in Release, but assert in Debug
#ifdef NDEBUG
    EXPECT_THROW(chopper::sketch::toolbox::read_hll_files_into(empty_dir.path(), test_filenames, target),
                 std::runtime_error);
#else
    EXPECT_DEATH(chopper::sketch::toolbox::read_hll_files_into(empty_dir.path(), test_filenames, target), "");
#endif
}

TEST_F(toolbox_test, read_hll_files_into_file_is_missing)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const tmp_file{tmp_dir.path() / "other.hll"};
    {
        std::ofstream os{tmp_file};
        os << "Doesn't matter I just need to exist\n";
    }

    std::vector<chopper::sketch::hyperloglog> target{};

    EXPECT_THROW(chopper::sketch::toolbox::read_hll_files_into(tmp_file.parent_path(), test_filenames, target),
                 std::runtime_error);
}

TEST_F(toolbox_test, read_hll_files_into_faulty_file)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const tmp_file{tmp_dir.path() / "small.hll"};
    {
        std::ofstream os{tmp_file};
        os << "I am not what an hll file looks like\n";
    }

    std::vector<chopper::sketch::hyperloglog> target{};

    EXPECT_THROW(chopper::sketch::toolbox::read_hll_files_into(tmp_file.parent_path(), test_filenames, target),
                 std::runtime_error);
}

TEST_F(toolbox_test, precompute_union_estimates_for)
{
    std::vector<uint64_t> estimates(4);

    chopper::sketch::toolbox::precompute_union_estimates_for(estimates, test_sketches, test_kmer_counts, 0);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{500, 0, 0, 0}));

    chopper::sketch::toolbox::precompute_union_estimates_for(estimates, test_sketches, test_kmer_counts, 1);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{670, 600, 0, 0}));

    chopper::sketch::toolbox::precompute_union_estimates_for(estimates, test_sketches, test_kmer_counts, 2);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{670, 670, 700, 0}));

    chopper::sketch::toolbox::precompute_union_estimates_for(estimates, test_sketches, test_kmer_counts, 3);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{670, 670, 670, 800}));
}

TEST_F(toolbox_test, random_shuffle)
{
    chopper::sketch::toolbox::prio_queue default_pq{};
    chopper::sketch::toolbox::distance_matrix dist{{0, default_pq},
                                                   {1, default_pq},
                                                   {2, default_pq},
                                                   {3, default_pq},
                                                   {4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> ids{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    ubs.random_shuffle(dist, ids);

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

TEST_F(toolbox_test, prune)
{
    chopper::sketch::toolbox::prio_queue default_pq{};
    chopper::sketch::toolbox::distance_matrix dist{{0, default_pq},
                                                   {1, default_pq},
                                                   {2, default_pq},
                                                   {3, default_pq},
                                                   {4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> remaining_ids{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    // since remaining_ids contains all_ids, prune shouldn't do anything. All ids are valid.
    ubs.prune(dist, remaining_ids);

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
    ubs.prune(dist, remaining_ids);

    EXPECT_EQ(remaining_ids[0], 0u);
    EXPECT_EQ(remaining_ids[2], 2u);
    EXPECT_EQ(remaining_ids[4], 1u);

    EXPECT_EQ(dist.size(), 3u);
    EXPECT_EQ(dist[0].id, 0u);
    EXPECT_EQ(dist[1].id, 4u);
    EXPECT_EQ(dist[2].id, 2u);
}

TEST_F(toolbox_test, rotate)
{
    chopper::sketch::hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (0,1)     (2,3)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<chopper::sketch::toolbox::clustering_node> clustering{{f, f, s},
                                                                      {f, f, s},
                                                                      {f, f, s},
                                                                      {f, f, s}, // the leaves come first
                                                                      {5, 6, s},
                                                                      {0, 1, s},
                                                                      {2, 3, s}};

    // previous_rightmost is already at the very left. Nothing has to be rotated.
    ubs.rotate(clustering, 0 /*previous_rightmost*/, 0 /*interval_start*/, 4 /*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(5u, 6u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));

    // now the previous_rightmost is within the tree. Rotation should take place
    ubs.rotate(clustering, 2 /*previous_rightmost*/, 0 /*interval_start*/, 4 /*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(6u, 5u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));
}

TEST_F(toolbox_test, trace)
{
    chopper::sketch::hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (1,3)     (2,0)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<chopper::sketch::toolbox::clustering_node> clustering{{f, f, s},
                                                                      {f, f, s},
                                                                      {f, f, s},
                                                                      {f, f, s}, // the leaves come first
                                                                      {5, 6, s},
                                                                      {1, 3, s},
                                                                      {2, 0, s}};

    std::vector<size_t> permutation{};

    ubs.trace(clustering, permutation, 2 /*previous_rightmost*/, 0 /*interval_start*/, 4 /*root_id*/);

    EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{1, 3, 0}));
}

TEST_F(toolbox_test, cluster_bins)
{
    chopper::sketch::toolbox ubs{test_filenames, test_kmer_counts, test_sketches};

    { // whole range
        std::vector<size_t> permutation{};
        ubs.cluster_bins(permutation, 0 /*interval start*/, 3 /*interval_end*/, 1 /*number of threads*/);
        // index 3 is not part of current permutation so it can participate in "the next interval"
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{2, 0, 1}));
    }

    { // intervals
        std::vector<size_t> permutation{};
        ubs.cluster_bins(permutation, 0 /*interval start*/, 1 /*interval_end*/, 1 /*number of threads*/);
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{0}));
        ubs.cluster_bins(permutation, 1 /*interval start*/, 3 /*interval_end*/, 1 /*number of threads*/);
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{0, 1, 2}));
    }
}
