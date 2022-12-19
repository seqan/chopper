#include <gtest/gtest.h>

#include <chopper/count/check_filenames.hpp>

TEST(check_filenames_test, sequence_filenames)
{
    std::vector<std::string> filenames{
        "/path/to/file1.fa",
        "/path/to/file2.fasta",
        "/path/to/file3.fq",
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::count::check_filenames(filenames, config));

    EXPECT_FALSE(config.precomputed_files);
}

TEST(check_filenames_test, minimizer_filenames)
{
    std::vector<std::string> filenames{
        "/path/to/file1.minimizer",
        "/path/to/file2.minimizer",
        "/path/to/file3.minimizer",
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::count::check_filenames(filenames, config));

    EXPECT_TRUE(config.precomputed_files);
}

TEST(check_filenames_test, mixed_filenames_sequence_files)
{
    std::vector<std::string> filenames{
        "/path/to/file2.fa",
        "/path/to/file1.minimizer",
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::count::check_filenames(filenames, config), std::invalid_argument);
}

TEST(check_filenames_test, mixed_filenames_minimizer_files)
{
    std::vector<std::string> filenames{
        "/path/to/file1.minimizer",
        "f.fa",
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::count::check_filenames(filenames, config), std::invalid_argument);
}
