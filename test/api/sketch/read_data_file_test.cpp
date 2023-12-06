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

#include <seqan3/test/tmp_directory.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include "../api_test.hpp"

TEST(read_data_file_test, file_open_error)
{
    chopper::configuration config{};
    std::vector<std::vector<std::string>> filenames{};
    config.data_file = data("non_existing.file");
    EXPECT_THROW(chopper::sketch::read_data_file(config, filenames), std::runtime_error);
}

TEST(read_data_file_test, small_example)
{
    chopper::configuration config;
    std::vector<std::vector<std::string>> filenames{};
    config.data_file = data("seqinfo.tsv");

    chopper::sketch::read_data_file(config, filenames);

    std::vector<std::vector<std::string>> expected_filenames{{"file1"}, {"file2"}, {"file3"}, {"file4"}, {"file5"}};
    EXPECT_RANGE_EQ(filenames, expected_filenames);
}

TEST(read_data_file_test, multi_filenames)
{
    chopper::configuration config;
    std::vector<std::vector<std::string>> filenames{};

    seqan3::test::tmp_directory tmp_dir{};
    config.data_file = tmp_dir.path() / "multi_files.txt";

    {
        std::ofstream of{config.data_file};
        of << "file1a file1b\nfile2\nfile3a file3b file3c\n";
    }

    chopper::sketch::read_data_file(config, filenames);

    std::vector<std::vector<std::string>> expected_filenames{{"file1a", "file1b"},
                                                             {"file2"},
                                                             {"file3a", "file3b", "file3c"}};
    EXPECT_RANGE_EQ(filenames, expected_filenames);
}
