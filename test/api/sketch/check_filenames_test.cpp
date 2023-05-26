#include <gtest/gtest.h>

#include <chopper/sketch/check_filenames.hpp>

#include "../api_test.hpp"

TEST(check_filenames_test, sequence_filenames)
{
    std::vector<std::string> filenames{data("seq1.fa").string(), data("seq2.fa").string(), data("seq3.fa").string()};

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::sketch::check_filenames(filenames, config));

    EXPECT_FALSE(config.precomputed_files);
}

TEST(check_filenames_test, minimiser_filenames)
{
    std::vector<std::string> filenames{data("small.minimiser").string(),
                                       data("small.minimiser").string(),
                                       data("small.minimiser").string()};

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::sketch::check_filenames(filenames, config));

    EXPECT_TRUE(config.precomputed_files);
}

TEST(check_filenames_test, mixed_filenames_sequence_files)
{
    std::vector<std::string> filenames{data("seq2.fa").string(), data("small.minimiser").string()};

    chopper::configuration config;

    EXPECT_THROW(chopper::sketch::check_filenames(filenames, config), std::invalid_argument);
}

TEST(check_filenames_test, mixed_filenames_minimiser_files)
{
    std::vector<std::string> filenames{data("small.minimiser").string(), data("seq2.fa").string()};

    chopper::configuration config;

    EXPECT_THROW(chopper::sketch::check_filenames(filenames, config), std::invalid_argument);
}
