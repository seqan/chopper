// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <algorithm> // for sort
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <cstdlib>   // for exit
#include <fstream>   // for ifstream
#include <numeric>   // for iota
#include <stdexcept> // for invalid_argument
#include <string>    // for string
#include <vector>    // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/input.hpp>
#include <chopper/layout/phibf/partition_user_bins.hpp>

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

#include "../api_test.hpp"

// Not declared in the header. Tested directly in phibf_regression_test.find_best_partition_intersection_underflow.
namespace chopper::layout::phibf
{
bool find_best_partition(chopper::configuration const & config,
                         size_t const cardinality_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                         std::vector<size_t> & partition_cardinality);
} // namespace chopper::layout::phibf

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

// `number_of_user_bins` user bins in families of 5 that mostly overlap, with sketches and the configuration for
// partition_user_bins. dist: "uniform" (3'000 to 23'000 hashes), "equal" (5'000 hashes), or "skewed" (two user bins
// with 400'000 hashes, the rest 3'000 to 3'600).
struct phibf_data
{
    std::vector<size_t> kmer_counts{};
    chopper::configuration config{};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    std::vector<size_t> cardinalities{};

    phibf_data(int const partitioning_approach,
               size_t const number_of_user_bins,
               size_t const number_of_partitions,
               std::string const & dist = "uniform")
    {
        for (size_t i = 0; i < number_of_user_bins; ++i)
        {
            if (dist == "equal")
                kmer_counts.push_back(5000);
            else if (dist == "skewed")
                kmer_counts.push_back((i < 2) ? 400'000 : 3000 + 100 * (i % 7));
            else
                kmer_counts.push_back(3000 + (scramble(i) % 20000));
        }

        config.number_of_partitions = number_of_partitions;
        config.partitioning_approach = partitioning_approach;
        config.hibf_config.number_of_user_bins = number_of_user_bins;
        config.hibf_config.tmax = 64;
        config.hibf_config.input_fn = [this](size_t const ub, seqan::hibf::insert_iterator it)
        {
            uint64_t const offset = static_cast<uint64_t>(ub / 5) << 32;
            for (uint64_t i = 0; i < kmer_counts[ub]; ++i)
                it = scramble(offset + i + (ub % 5) * 97);
        };
        config.hibf_config.validate_and_set_defaults();

        seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

        for (auto const & sketch : sketches)
            cardinalities.push_back(sketch.estimate());
    }

    std::vector<std::vector<size_t>> partition() const
    {
        std::vector<std::vector<size_t>> partitions(config.number_of_partitions);
        chopper::layout::phibf::partition_user_bins(config, cardinalities, sketches, minHash_sketches, partitions);
        return partitions;
    }

    std::vector<seqan::hibf::layout::layout> layouts()
    {
        seqan3::test::tmp_directory tmp_dir{};
        config.output_filename = tmp_dir.path() / "dummy.layout";
        std::vector<std::vector<std::string>> filenames(kmer_counts.size(), {"ub.fa"});
        chopper::layout::execute(config, filenames, sketches, minHash_sketches);

        std::ifstream layout_stream{config.output_filename};
        return std::get<2>(chopper::layout::read_layouts_file(layout_stream));
    }
};

std::vector<chopper::layout::phibf::partitioning_scheme> const phibf_approaches
{
    chopper::layout::phibf::partitioning_scheme::blocked,
    chopper::layout::phibf::partitioning_scheme::sorted,
    chopper::layout::phibf::partitioning_scheme::folded,
    chopper::layout::phibf::partitioning_scheme::weighted_fold,
    chopper::layout::phibf::partitioning_scheme::similarity,
    chopper::layout::phibf::partitioning_scheme::lsh,
    chopper::layout::phibf::partitioning_scheme::lsh_sim
};

// Expects that every user bin in [0, number_of_user_bins) is assigned to exactly one partition.
bool each_user_bin_assigned_once(std::vector<std::vector<size_t>> const & partitions, size_t const number_of_user_bins)
{
    std::vector<size_t> assigned_user_bins{};
    for (auto const & partition : partitions)
        assigned_user_bins.insert(assigned_user_bins.end(), partition.begin(), partition.end());
    std::ranges::sort(assigned_user_bins);

    std::vector<size_t> expected_user_bins(number_of_user_bins);
    std::iota(expected_user_bins.begin(), expected_user_bins.end(), 0u);
    return assigned_user_bins == expected_user_bins;
}

} // namespace

TEST(phibf_partition_test, each_user_bin_is_assigned_to_exactly_one_partition)
{
    size_t const number_of_user_bins{240};
    size_t const number_of_partitions{4};

    for (chopper::layout::phibf::partitioning_scheme const approach : phibf_approaches)
    {
        phibf_data data{approach, number_of_user_bins, number_of_partitions};
        auto partitions = data.partition();
        ASSERT_EQ(partitions.size(), number_of_partitions);
        EXPECT_TRUE(each_user_bin_assigned_once(partitions, number_of_user_bins));
    }
}

TEST(phibf_partition_test, unknown_partitioning_approach)
{
    phibf_data data{-1 /*unkown approach*/, 16, 4};
    EXPECT_THROW(data.partition(), std::invalid_argument);
}

TEST(phibf_execute_test, writes_one_layout_per_partition)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "phibf.layout"};

    size_t const number_of_user_bins{240};
    size_t const number_of_partitions{3};
    phibf_data data{chopper::layout::phibf::partitioning_scheme::lsh_sim, number_of_user_bins, number_of_partitions};
    data.config.output_filename = layout_file;

    std::vector<std::vector<std::string>> filenames{};
    for (size_t ub = 0; ub < number_of_user_bins; ++ub)
        filenames.push_back({"ub" + std::to_string(ub) + ".fa"});

    chopper::layout::execute(data.config, filenames, data.sketches, data.minHash_sketches);

    std::ifstream layout_stream{layout_file};
    auto [read_filenames, read_config, layouts] = chopper::layout::read_layouts_file(layout_stream);

    EXPECT_EQ(read_filenames, filenames);
    EXPECT_EQ(read_config.number_of_partitions, number_of_partitions);
    EXPECT_EQ(read_config.partitioning_approach, chopper::layout::phibf::partitioning_scheme::lsh_sim);
    ASSERT_EQ(layouts.size(), number_of_partitions);

    // Each layout holds the user bins of one partition, with global user bin indices.
    std::vector<std::vector<size_t>> partitions{};
    for (auto const & layout : layouts)
    {
        EXPECT_FALSE(layout.user_bins.empty());
        partitions.emplace_back();
        for (auto const & user_bin : layout.user_bins)
            partitions.back().push_back(user_bin.idx);
    }
    EXPECT_TRUE(each_user_bin_assigned_once(partitions, number_of_user_bins));
}

TEST(phibf_invalid_input_test, fewer_user_bins_than_partitions)
{
    for (chopper::layout::phibf::partitioning_scheme const approach : phibf_approaches)
    {
        phibf_data const data{approach, 1/*UB*/, 2/*partitions*/, "uniform"};
        EXPECT_THROW(data.partition(), std::invalid_argument) << "approach " << approach;
    }
}

TEST(phibf_regression_test, small_number_of_user_bins_per_part)
{
    // If for small number of user bins the block size could not be optimally chosen (runding errors)
    // The last partition ended up empty. This was fixed.
    phibf_data const data{chopper::layout::phibf::partitioning_scheme::blocked, 100, 16, "uniform"};
    EXPECT_TRUE(each_user_bin_assigned_once(data.partition(), 100));
}

TEST(phibf_regression_test, few_user_bins_with_skewed_sizes)
{
    // If few user bins with skewed sizes are the input,
    // then cardinality_per_part is chosen too high and partitions at the end get no user bins.
    // This has beenf fixed.
    phibf_data const data{chopper::layout::phibf::partitioning_scheme::sorted, 64, 8, "skewed"};
    EXPECT_TRUE(each_user_bin_assigned_once(data.partition(), 64));
}

TEST(phibf_regression_test, zero_cardinality_user_bins)
{
    using chopper::layout::phibf::partitioning_scheme;

    for (chopper::layout::phibf::partitioning_scheme const approach : phibf_approaches)
    {
        phibf_data data{approach, 5, 2, "equal"};
        data.cardinalities = {10, 10, 10, 10, 0};

        auto const partitions = data.partition();
        EXPECT_TRUE(each_user_bin_assigned_once(partitions, 5));
    }
}

TEST(phibf_regression_test, fewer_clusters_than_partitions)
{
    using chopper::layout::phibf::partitioning_scheme;

    for (int const approach : {partitioning_scheme::lsh, partitioning_scheme::lsh_sim})
    {
        phibf_data const data{approach, 9, 8, "equal"};
        EXPECT_TRUE(each_user_bin_assigned_once(data.partition(), 8));
    }
}

TEST(phibf_regression_test, user_bin_does_not_fit)
{
    // which those inputs, the respective partitioning approach failes.
    phibf_data const data_lsh{chopper::layout::phibf::partitioning_scheme::lsh, 8, 5, "equal"};
    EXPECT_TRUE(each_user_bin_assigned_once(data_lsh.partition(), 8));

    phibf_data const data_lsh_sim{chopper::layout::phibf::partitioning_scheme::lsh_sim, 257, 8, "skewed"};
    EXPECT_TRUE(each_user_bin_assigned_once(data_lsh_sim.partition(), 257));

    phibf_data const data_similarity{chopper::layout::phibf::partitioning_scheme::similarity, 4, 3, "equal"};
    EXPECT_TRUE(each_user_bin_assigned_once(data_similarity.partition(), 4));
}

TEST(phibf_regression_test, weighted_fold_crash)
{
    // About two user bins per partition: some partitions stay empty, which is reported.
    phibf_data const data3{chopper::layout::phibf::partitioning_scheme::weighted_fold, 33, 16, "uniform"};
    EXPECT_TRUE(each_user_bin_assigned_once(data3.partition(), 33));
}

TEST(phibf_regression_test, weighted_fold_score_edge_case)
{
    // 8 user bins, 2 partitions: 125 k-mers and 4 user bins per partition.
    phibf_data data{chopper::layout::phibf::partitioning_scheme::weighted_fold, 8, 2, "equal"};
    data.cardinalities = {100, 90, 10, 10, 10, 10, 10, 10};

    // Partition 0 takes the large user bins 100 and 90 (score |1 - 190/125| + |1 - 2/4| = 1.02).
    // Removing 90 and adding small user bins while the score improves: 110/2 (0.62), 120/3 (0.29), 130/4 (0.04);
    // a fourth one would be worse (140/5: 0.37). 0.04 < 1.02, so this is accepted.
    // Removing 100 as well cannot improve on that (at best 30/4: 0.76). The rest goes to partition 1.
    auto const partitions = data.partition();
    ASSERT_EQ(partitions.size(), 2u);

    auto partition_cardinalities = [&data](std::vector<size_t> const & partition)
    {
        std::vector<size_t> result;
        for (size_t const user_bin : partition)
            result.push_back(data.cardinalities[user_bin]);
        std::ranges::sort(result, std::ranges::greater{});
        return result;
    };

    EXPECT_EQ(partition_cardinalities(partitions[0]), (std::vector<size_t>{100, 10, 10, 10}));
    EXPECT_EQ(partition_cardinalities(partitions[1]), (std::vector<size_t>{90, 10, 10, 10}));
    EXPECT_TRUE(each_user_bin_assigned_once(partitions, 8));
}

TEST(phibf_execute_test, tmaxset_or_not_set)
{
    // The number of technical bins of the top-level IBF of a layout.
    auto top_level_technical_bins = [](seqan::hibf::layout::layout const & layout)
    {
        size_t result{};
        for (auto const & user_bin : layout.user_bins)
        {
            if (user_bin.previous_TB_indices.empty()) // single or split bin on the top level
                result = std::max(result, user_bin.storage_TB_id + user_bin.number_of_technical_bins);
            else // merged bin on the top level
                result = std::max(result, user_bin.previous_TB_indices[0] + 1);
        }
        return result;
    };

    size_t const number_of_partitions{3};

    // tmax not set: about 80 user bins per partition, next_multiple_of_64(ceil(sqrt(80))) = 64.
    {
        phibf_data data{chopper::layout::phibf::partitioning_scheme::sorted, 240, number_of_partitions};
        auto const layouts = data.layouts();
        ASSERT_EQ(layouts.size(), number_of_partitions);
        for (auto const & layout : layouts)
            EXPECT_EQ(top_level_technical_bins(layout), 64u);
    }

    // tmax set: every partition uses it.
    {
        phibf_data data{chopper::layout::phibf::partitioning_scheme::sorted, 240, number_of_partitions};
        data.config.hibf_config.tmax = 128;
        data.config.tmax_is_set = true;
        auto const layouts = data.layouts();
        ASSERT_EQ(layouts.size(), number_of_partitions);
        for (auto const & layout : layouts)
            EXPECT_EQ(top_level_technical_bins(layout), 128u);
    }
}

namespace
{

// A HyperLogLog sketch of the hashes scramble(offset), ..., scramble(offset + count - 1).
seqan::hibf::sketch::hyperloglog make_sketch(uint64_t const offset, size_t const count, uint8_t const sketch_bits)
{
    seqan::hibf::sketch::hyperloglog sketch{sketch_bits};
    for (uint64_t i = 0; i < count; ++i)
        sketch.add(scramble(offset + i));
    return sketch;
}

} // namespace

// lsh and lsh_sim initialise the partitions with the first clusters, which post_process_clusters orders by their
// number of user bins, not by cardinality. A cluster that exceeds the target cardinality per partition spills over
// into the next partition, but `++p` was only guarded by an assertion: a heavy cluster in the last seeded position
// moved p past the last partition (assertion in Debug, out-of-bounds write in Release).
TEST(phibf_regression_test, seeding_spills_past_last_partition)
{
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    using chopper::layout::phibf::partitioning_scheme;

    // User bins 0-3 (cardinality 1) and 4-6 (cardinality 60) have identical MinHash tables within each group and form
    // the clusters {0, 1, 2, 3} and {4, 5, 6}. The first has more user bins and seeds partition 0, the second seeds
    // partition 1. The target cardinality per partition is ceil(184 / 2) = 92, which {4, 5, 6} exceeds after two of
    // its user bins, so the third one spilled over into partition 2 of 2.
    size_t const number_of_user_bins{7};
    std::vector<size_t> const cardinalities{1, 1, 1, 1, 60, 60, 60};

    chopper::configuration config{};
    config.number_of_partitions = 2;
    config.hibf_config.number_of_user_bins = number_of_user_bins;
    config.hibf_config.tmax = 64;

    auto minhashes_of_group = [](uint64_t const group)
    {
        using seqan::hibf::sketch::minhashes;
        minhashes result{};
        result.table.assign(minhashes::num_sketches, std::vector<uint64_t>(minhashes::sketch_size));
        uint64_t value{group << 32};
        for (auto & sketch : result.table)
            for (auto & hash : sketch)
                hash = scramble(value++);
        return result;
    };

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    for (size_t ub = 0; ub < number_of_user_bins; ++ub)
    {
        sketches.push_back(make_sketch(uint64_t{ub} << 32, cardinalities[ub], config.hibf_config.sketch_bits));
        minHash_sketches.push_back(minhashes_of_group(ub < 4 ? 1 : 2));
    }

    for (int const approach : {partitioning_scheme::lsh, partitioning_scheme::lsh_sim})
    {
        config.partitioning_approach = approach;

        // Debug builds abort, so this runs in a child process. It exits with 0 if every user bin is assigned once.
        EXPECT_EXIT(
            {
                std::vector<std::vector<size_t>> partitions(config.number_of_partitions);
                chopper::layout::phibf::partition_user_bins(config,
                                                            cardinalities,
                                                            sketches,
                                                            minHash_sketches,
                                                            partitions);
                std::exit(each_user_bin_assigned_once(partitions, number_of_user_bins) ? 0 : 1);
            },
            ::testing::ExitedWithCode(0),
            "")
            << "approach " << approach;
    }
}

// find_best_partition estimated the intersection of a cluster and a partition as
// est(cluster) - (est(union) - est(partition)) in size_t. For disjoint sets, HyperLogLog noise makes the difference
// larger than est(cluster) about half of the time. The subtraction then wrapped around, and the partition got a huge
// similarity, so user bins went to a partition they share nothing with.
TEST(phibf_regression_test, find_best_partition_intersection_underflow)
{
    chopper::configuration config{};
    config.number_of_partitions = 2;
    uint8_t const sketch_bits{config.hibf_config.sketch_bits};

    // Partition 0 holds 200'000 k-mers that no user bin shares. Partition 1 holds 100'000 k-mers, and each user bin is
    // a disjoint subset of 1'000 of them. Each user bin must go to partition 1: there, the union is the partition
    // itself, so the similarity is est(user bin) / est(partition 1). Partition 0 can score at most
    // est(user bin) / est(partition 0), which is smaller.
    uint64_t const offset_1{uint64_t{1} << 40};
    seqan::hibf::sketch::hyperloglog const partition_0 = make_sketch(0, 200'000, sketch_bits);
    seqan::hibf::sketch::hyperloglog const partition_1 = make_sketch(offset_1, 100'000, sketch_bits);

    size_t const number_of_user_bins{50};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<size_t> cardinalities{};
    for (size_t ub = 0; ub < number_of_user_bins; ++ub)
    {
        sketches.push_back(make_sketch(offset_1 + ub * 1'000, 1'000, sketch_bits));
        cardinalities.push_back(sketches.back().estimate());
    }

    size_t wrong_partition{0};
    for (size_t ub = 0; ub < number_of_user_bins; ++ub)
    {
        std::vector<std::vector<size_t>> partitions(config.number_of_partitions);
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches{partition_0, partition_1};
        std::vector<size_t> partition_cardinality{0, 0}; // Both partitions have room for the user bin.

        bool const found = chopper::layout::phibf::find_best_partition(config,
                                                                       /*cardinality_per_part*/ 1'000'000,
                                                                       {ub},
                                                                       cardinalities,
                                                                       sketches,
                                                                       partitions,
                                                                       partition_sketches,
                                                                       partition_cardinality);
        EXPECT_TRUE(found) << "user bin " << ub;
        wrong_partition += (partitions[1] != std::vector<size_t>{ub});
    }
    EXPECT_EQ(wrong_partition, 0u) << "of " << number_of_user_bins << " user bins";
}

// The check that every user bin is assigned exactly once was only compiled in Debug builds, although
// partition_user_bins documents a std::logic_error. lsh bounds its last loop by config.hibf_config.number_of_user_bins.
// If that is not set, the remaining clusters are never assigned, and Release builds returned the incomplete partitions.
TEST(phibf_regression_test, unassigned_user_bins_in_release)
{
#ifndef NDEBUG
    GTEST_SKIP() << "Debug builds have the check, and an assertion on number_of_user_bins fires before it.";
#else
    size_t const number_of_user_bins{50};
    phibf_data data{chopper::layout::phibf::partitioning_scheme::lsh, number_of_user_bins, 2, "uniform"};
    data.config.hibf_config.number_of_user_bins = 0;

    try
    {
        // Either all user bins are assigned ...
        EXPECT_TRUE(each_user_bin_assigned_once(data.partition(), number_of_user_bins));
    }
    catch (std::logic_error const &)
    {
        // ... or partition_user_bins reports that they are not.
    }
#endif
}
