// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <chopper/sketch/check_filenames.hpp>

#include "../api_test.hpp"

TEST(check_filenames_test, sequence_filenames)
{
    std::vector<std::string> filenames{data("seq1.fa").string(), data("seq2.fa").string(), data("seq3.fa").string()};

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::sketch::check_filenames(filenames, config));

    EXPECT_FALSE(config.precomputed_files);
}

TEST(check_filenames_test, overload)
{
    std::vector<std::vector<std::string>> filenames{{data("seq1.fa").string()},
                                                    {data("seq2.fa").string()},
                                                    {data("seq3.fa").string()}};

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
