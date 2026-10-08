// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <algorithm> // for sort, unique
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <cstdlib>   // for exit
#include <numeric>   // for iota
#include <string>    // for to_string
#include <vector>    // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/lsh_distributed_ibf_layout.hpp>
#include <chopper/workarounds.hpp>

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
// The cardinalities of the user bins in `zero_cardinality_user_bins` are set to 0, regardless of their content.
std::vector<std::vector<size_t>>
run_lsh_distributed_ibf_layout(std::vector<size_t> const & kmer_counts,
                               size_t const tmax,
                               std::vector<size_t> content_ids = {},
                               std::vector<size_t> const & zero_cardinality_user_bins = {})
{
    if (content_ids.empty())
    {
#if CHOPPER_WORKAROUND_GCC_BOGUS_ARRAY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Warray-bounds="
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_ARRAY
        content_ids.resize(kmer_counts.size());
#if CHOPPER_WORKAROUND_GCC_BOGUS_ARRAY
#    pragma GCC diagnostic pop
#endif // CHOPPER_WORKAROUND_GCC_BOGUS_ARRAY
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
    for (size_t const ub : zero_cardinality_user_bins)
        cardinalities[ub] = 0u;

    std::vector<size_t> positions(kmer_counts.size());
    std::iota(positions.begin(), positions.end(), 0u);

    std::vector<std::vector<size_t>> technical_bins(tmax);
    chopper::layout::lsh_distributed_ibf_layout(config,
                                                positions,
                                                cardinalities,
                                                sketches,
                                                minHash_sketches,
                                                technical_bins);

    return technical_bins;
}

// Expects that every user bin in [0, number_of_user_bins) is assigned to at least one technical bin.
void expect_all_user_bins_assigned(std::vector<std::vector<size_t>> const & technical_bins,
                                   size_t const number_of_user_bins)
{
    std::vector<size_t> assigned_user_bins{};
    for (auto const & technical_bin : technical_bins)
        assigned_user_bins.insert(assigned_user_bins.end(), technical_bin.begin(), technical_bin.end());
    std::ranges::sort(assigned_user_bins);
    auto const [first, last] = std::ranges::unique(assigned_user_bins);
    assigned_user_bins.erase(first, last);

    std::vector<size_t> expected_user_bins(number_of_user_bins);
    std::iota(expected_user_bins.begin(), expected_user_bins.end(), 0u);
    EXPECT_EQ(assigned_user_bins, expected_user_bins);
}

} // namespace

// A user bin is split if it is assigned to more than one technical bin.
// A technical bin is merged if it contains more than one user bin.

TEST(lsh_distributed_ibf_layout_test, only_split_bins)
{
    // 4 large user bins for 8 technical bins: every user bin exceeds the split threshold.
    auto const technical_bins = run_lsh_distributed_ibf_layout(std::vector<size_t>(4, 10'000), /*tmax*/ 8);

    // each user bin is split into two technical bins
    std::vector<std::vector<size_t>> expected_technical_bins{{3}, {3}, {1}, {1}, {2}, {2}, {0}, {0}};

    ASSERT_EQ(technical_bins.size(), expected_technical_bins.size());
    for (size_t tb = 0; tb < technical_bins.size(); ++tb)
        EXPECT_EQ(technical_bins[tb], expected_technical_bins[tb]) << "technical bin " << tb << " is not correct";
}

TEST(lsh_distributed_ibf_layout_test, initial_split_threshold_small)
{
    std::vector<size_t> kmer_counts(2, 200'000);
    kmer_counts.resize(12, 1000);
    std::vector<size_t> content_ids{0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // user bins 0 and 1 are identical

    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, /*tmax*/ 7, content_ids);

    // ToDo this distribution is not perfect
    std::vector<std::vector<size_t>> expected_technical_bins{{8, 4, 6, 10, 9, 3}, {5, 7}, {11, 2}, {0}, {0}, {1}, {1}};

    ASSERT_EQ(technical_bins.size(), expected_technical_bins.size());
    for (size_t tb = 0; tb < technical_bins.size(); ++tb)
        EXPECT_EQ(technical_bins[tb], expected_technical_bins[tb]) << "technical bin " << tb << " is not correct";
}

TEST(lsh_distributed_ibf_layout_test, only_merged_bins)
{
    std::vector<size_t> kmer_counts(200, 2'000);
    // 20 small user bins for 2 technical bins: no use
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, 2);

    // since all user bins have the same cardinality, both technical bins should contain equally many
    // SInce not cardinalities but HLL sketches are used, there is some noise
    // Thats why there are not exactly equal
    ASSERT_EQ(technical_bins.size(), 2);
#ifdef _LIBCPP_VERSION
    EXPECT_EQ(technical_bins[0].size(), 99);
    EXPECT_EQ(technical_bins[1].size(), 101);
#else
    EXPECT_EQ(technical_bins[0].size(), 102);
    EXPECT_EQ(technical_bins[1].size(), 98);
#endif
}

TEST(lsh_distributed_ibf_layout_test, another_edge_case)
{
    // one very large user bin and one very small one.
    std::vector<size_t> const kmer_counts{400'000, 2'000};

    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, 64);

    std::vector<std::vector<size_t>> expected_technical_bins{{8, 4, 6, 10, 9, 3}, {5, 7}, {11, 2}, {0}, {0}, {1}, {1}};

    ASSERT_EQ(technical_bins.size(), 64);

    ASSERT_EQ(technical_bins[0].size(), 1);
    EXPECT_EQ(technical_bins[0][0], 1);
    for (size_t tb = 1; tb < technical_bins.size(); ++tb)
    {
        ASSERT_EQ(technical_bins[tb].size(), 1);
        EXPECT_EQ(technical_bins[tb][0], 0) << "technical bin " << tb << " is not correct";
    }
}

// KI tests that failed and caught some bugs that are now fixed:

TEST(lsh_distributed_ibf_layout_test, one_large_cluster)
{
    // This case simulates that there are more available merged bins than clusters
    // with only one large cluster
    std::vector<size_t> kmer_counts{50'635, 16'021, 1'734, 305'923, 1'844, 7'087};
    std::vector<size_t> content_ids{0, 0, 0, 0, 0, 0};
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, 64, content_ids);

    ASSERT_EQ(technical_bins.size(), 64);
}

TEST(lsh_distributed_ibf_layout_test, one_large_one_small_cluster)
{
    // This case simulates that there are more available merged bins than clusters
    // with two clusters, of which one only has size 1
    std::vector<size_t> kmer_counts{4'169, 14'011, 4'311, 2'016, 239'260, 16'111};
    std::vector<size_t> content_ids{0, 0, 0, 0, 1, 0};
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, 64, content_ids);

    ASSERT_EQ(technical_bins.size(), 64);
}

TEST(lsh_distributed_ibf_layout_test, several_clusters)
{
    // This case simulates that there are more available merged bins than clusters
    // with several clusters
    std::vector<size_t> kmer_counts{5'829, 1'875, 8'912, 254'423, 6'445, 188'171, 7'621, 7'127, 7'579, 11'915, 344'000};
    std::vector<size_t> content_ids{1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1};
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, 64, content_ids);

    ASSERT_EQ(technical_bins.size(), 64);
}

TEST(lsh_distributed_ibf_layout_test, cluster_larger_than_tmax_with_small_cardinality)
{
    // User bins 2 to 6 are identical and form one cluster of 5 > tmax user bins. Its cardinality is below
    // 0.05 * sum_of_cardinalities / tmax, so only tmax user bins seed a technical bin. The others must still be
    // assigned. User bins 7 to 9 provide enough clusters, so the large cluster is not broken up beforehand.
    std::vector<size_t> const kmer_counts{300'000, 300'000, 1'000, 1'000, 1'000, 1'000, 1'000, 1'000, 1'000, 1'000};
    std::vector<size_t> const content_ids{0, 1, 2, 2, 2, 2, 2, 3, 4, 5};
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, /*tmax*/ 4, content_ids);

    expect_all_user_bins_assigned(technical_bins, kmer_counts.size());
}

TEST(lsh_distributed_ibf_layout_test, user_bin_with_zero_cardinality)
{
    // User bins 0 to 11 form four clusters of three identical user bins, which leaves eight empty (moved) clusters.
    // User bin 14 has an estimated cardinality of 0. Its non-empty cluster must still be sorted before the empty
    // clusters, which also have the key 0, or it is not assigned.
    std::vector<size_t> const kmer_counts(15, 3'000);
    std::vector<size_t> const content_ids{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 6};
    auto const technical_bins = run_lsh_distributed_ibf_layout(kmer_counts, /*tmax*/ 4, content_ids, {14});

    expect_all_user_bins_assigned(technical_bins, kmer_counts.size());
}
