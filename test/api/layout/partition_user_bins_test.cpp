// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t
#include <cstdlib> // for exit
#include <numeric> // for iota
#include <string>  // for to_string
#include <vector>  // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/partition_user_bins.hpp>

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace
{

// splitmix64 finaliser: turns consecutive integers into well-distributed hashes. HyperLogLog and the MinHash buckets
// (hash & 15, each needs 40 values) both assume uniformly distributed hashes, which plain consecutive integers are not.
uint64_t scramble(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// User bin `i` consists of `kmer_counts[i]` distinct hashes drawn from content `content_ids[i]`.
// User bins with the same content id share their first min(kmer_counts) hashes, all others are disjoint.
// If `content_ids` is empty, every user bin has its own content.
std::vector<std::vector<size_t>> run_partition_user_bins(std::vector<size_t> const & kmer_counts,
                                                         size_t const tmax,
                                                         std::vector<size_t> content_ids = {})
{
    if (content_ids.empty())
    {
        content_ids.resize(kmer_counts.size());
        std::iota(content_ids.begin(), content_ids.end(), 0u);
    }

    chopper::configuration config{};
    config.hibf_config.tmax = tmax;
    config.hibf_config.number_of_user_bins = kmer_counts.size();
    config.hibf_config.input_fn = [&kmer_counts, &content_ids](size_t const ub, seqan::hibf::insert_iterator it)
    {
        uint64_t const offset = static_cast<uint64_t>(content_ids[ub]) << 32;
        for (uint64_t i = 0; i < kmer_counts[ub]; ++i)
            it = scramble(offset + i);
    };

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

    std::vector<size_t> cardinalities(kmer_counts.size());
    for (size_t i = 0; i < sketches.size(); ++i)
        cardinalities[i] = sketches[i].estimate();

    std::vector<size_t> positions(kmer_counts.size());
    std::iota(positions.begin(), positions.end(), 0u);

    std::vector<std::vector<size_t>> partitions(tmax);
    chopper::layout::partition_user_bins(config, positions, cardinalities, sketches, minHash_sketches, partitions);

    return partitions;
}

} // namespace

// A user bin is split if it is assigned to more than one technical bin.
// A technical bin is merged if it contains more than one user bin.

TEST(partition_user_bins_test, only_split_bins)
{
    // 4 large user bins for 8 technical bins: every user bin exceeds the split threshold.
    auto const partitions = run_partition_user_bins(std::vector<size_t>(4, 10'000), /*tmax*/ 8);

    // each user bin is split into two technical bins
    std::vector<std::vector<size_t>> expected_partitions{{3}, {3}, {1}, {1}, {2}, {2}, {0}, {0}};

    ASSERT_EQ(partitions.size(), expected_partitions.size());
    for (size_t tb = 0; tb < partitions.size(); ++tb)
        EXPECT_EQ(partitions[tb], expected_partitions[tb]) << "technical bin " << tb << " is not correct";
}

TEST(partition_user_bins_test, initial_split_threshold_small)
{
    std::vector<size_t> kmer_counts(2, 200'000);
    kmer_counts.resize(12, 1000);
    std::vector<size_t> content_ids{0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // user bins 0 and 1 are identical

    auto const partitions = run_partition_user_bins(kmer_counts, /*tmax*/ 7, content_ids);

    // ToDo this partitioning is not perfect
    std::vector<std::vector<size_t>> expected_partitions{{8, 4, 6, 10, 9, 3}, {5, 7}, {11, 2}, {0}, {0}, {1}, {1}};

    ASSERT_EQ(partitions.size(), expected_partitions.size());
    for (size_t tb = 0; tb < partitions.size(); ++tb)
        EXPECT_EQ(partitions[tb], expected_partitions[tb]) << "technical bin " << tb << " is not correct";
}

TEST(partition_user_bins_test, only_merged_bins)
{
    std::vector<size_t> kmer_counts(200, 2'000);
    // 20 small user bins for 2 technical bins: no use
    auto const partitions = run_partition_user_bins(kmer_counts, 2);

    // since all user bins have the same cardinality, both partitions should contain equally many
    // SInce not cardinalities but HLL sketches are used, there is some noise
    // Thats why there are not exactly equal
    ASSERT_EQ(partitions.size(), 2);
    EXPECT_EQ(partitions[0].size(), 102);
    EXPECT_EQ(partitions[1].size(), 98);
}

TEST(partition_user_bins_test, another_edge_case)
{
    // one very large user bin and one very small one.
    std::vector<size_t> const kmer_counts{400'000, 2'000};

    auto const partitions = run_partition_user_bins(kmer_counts, 64);

    std::vector<std::vector<size_t>> expected_partitions{{8, 4, 6, 10, 9, 3}, {5, 7}, {11, 2}, {0}, {0}, {1}, {1}};

    ASSERT_EQ(partitions.size(), 64);

    ASSERT_EQ(partitions[0].size(), 1);
    EXPECT_EQ(partitions[0][0], 1);
    for (size_t tb = 1; tb < partitions.size(); ++tb)
    {
        ASSERT_EQ(partitions[tb].size(), 1);
        EXPECT_EQ(partitions[tb][0], 0) << "technical bin " << tb << " is not correct";
    }
}

// KI tests that failed and caught some bugs that are now fixed:

TEST(partition_user_bins_test, one_large_cluster)
{
    // This case simulates that there are more available merged bins than clusters
    // with only one large cluster
    std::vector<size_t> kmer_counts{50'635, 16'021, 1'734, 305'923, 1'844, 7'087};
    std::vector<size_t> content_ids{0, 0, 0, 0, 0, 0};
    auto const partitions = run_partition_user_bins(kmer_counts, 64, content_ids);

    ASSERT_EQ(partitions.size(), 64);
}

TEST(partition_user_bins_test, one_large_one_small_cluster)
{
    // This case simulates that there are more available merged bins than clusters
    // with two clusters, of which one only has size 1
    std::vector<size_t> kmer_counts{4'169, 14'011, 4'311, 2'016, 239'260, 16'111};
    std::vector<size_t> content_ids{0, 0, 0, 0, 1, 0};
    auto const partitions = run_partition_user_bins(kmer_counts, 64, content_ids);

    ASSERT_EQ(partitions.size(), 64);
}

TEST(partition_user_bins_test, several_clusters)
{
    // This case simulates that there are more available merged bins than clusters
    // with several clusters
    std::vector<size_t> kmer_counts{5'829, 1'875, 8'912, 254'423, 6'445, 188'171, 7'621, 7'127, 7'579, 11'915, 344'000};
    std::vector<size_t> content_ids{1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1};
    auto const partitions = run_partition_user_bins(kmer_counts, 64, content_ids);

    ASSERT_EQ(partitions.size(), 64);
}
