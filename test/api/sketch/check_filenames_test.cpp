#include <gtest/gtest.h>

#include <chopper/sketch/check_filenames.hpp>

TEST(check_filenames_test, sequence_filenames)
{
    std::vector<std::string> filenames{
        "/path/to/file1.fa",
        "/path/to/file2.fasta",
        "/path/to/file3.fq",
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::sketch::check_filenames(filenames, config));

    EXPECT_FALSE(config.precomputed_files);
}

TEST(check_filenames_test, minimiser_filenames)
{
    std::vector<std::string> filenames{
        "/path/to/file1.minimiser",
        "/path/to/file2.minimiser",
        "/path/to/file3.minimiser",
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::sketch::check_filenames(filenames, config));

    EXPECT_TRUE(config.precomputed_files);
}

TEST(check_filenames_test, mixed_filenames_sequence_files)
{
    std::vector<std::string> filenames{
        "/path/to/file2.fa",
        "/path/to/file1.minimiser",
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::sketch::check_filenames(filenames, config), std::invalid_argument);
}

TEST(check_filenames_test, mixed_filenames_minimiser_files)
{
    std::vector<std::string> filenames{
        "/path/to/file1.minimiser",
        "f.fa",
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::sketch::check_filenames(filenames, config), std::invalid_argument);
}
