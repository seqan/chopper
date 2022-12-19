#include <gtest/gtest.h>

#include <cmath>

#include "../api_test.hpp"
#include <chopper/configuration.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/detail_apply_prefix.hpp>

TEST(count_kmers_test, small_example)
{
    chopper::configuration config;
    config.k = 15;
    config.threads = 1;
    config.disable_sketch_output = true;

    chopper::data_store store{.filenames = {data("small.fa")}};

    chopper::count::count_kmers(config, store);

    ASSERT_EQ(store.sketches.size(), 1);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);
}

TEST(count_kmers_test, with_file_output)
{
    seqan3::test::tmp_filename sketch_dir{"small"};

    chopper::configuration config;
    config.k = 15;
    config.threads = 1;
    config.sketch_directory = sketch_dir.get_path().string();

    chopper::data_store store{.filenames = {data("small.fa")}};

    chopper::count::count_kmers(config, store);

    ASSERT_EQ(store.sketches.size(), 1);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);

    // check files
    ASSERT_TRUE(std::filesystem::exists(config.sketch_directory));
    ASSERT_TRUE(std::filesystem::exists(config.sketch_directory / "small.hll"));

    std::ifstream hll_file(config.sketch_directory / "small.hll", std::ios::binary);
    chopper::sketch::hyperloglog sketch;
    sketch.restore(hll_file); // the sketch bits will be automatically read from the files
    EXPECT_EQ(std::lround(sketch.estimate()), 571);
}

TEST(count_kmers_test, small_example_parallel_2_threads)
{
    chopper::configuration config;
    config.k = 15;
    config.threads = 2;
    config.disable_sketch_output = true;

    chopper::data_store store{.filenames = {data("small.fa"), data("small.fa"), data("small.fa"), data("small.fa")}};

    chopper::count::count_kmers(config, store);

    ASSERT_EQ(store.sketches.size(), 4);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);
    EXPECT_EQ(std::lround(store.sketches[1].estimate()), 571);
    EXPECT_EQ(std::lround(store.sketches[2].estimate()), 571);
    EXPECT_EQ(std::lround(store.sketches[3].estimate()), 571);
}

TEST(count_kmers_test, read_in_precomputed_binary_files)
{
    chopper::configuration config;
    config.k = 15;
    config.threads = 1;
    config.precomputed_files = true;
    config.disable_sketch_output = true;

    chopper::data_store store{.filenames = {data("small.minimizer")}};

    chopper::count::count_kmers(config, store);

    ASSERT_EQ(store.sketches.size(), 1);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);
}
