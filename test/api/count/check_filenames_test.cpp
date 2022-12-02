#include <gtest/gtest.h>

#include <chopper/count/check_filenames.hpp>

TEST(check_filenames_test, sequence_filenames)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters
    {
        {"key1", {"/path/to/file1.fa"}},
        {"key2", {"/path/to/file2.fasta", "/path/to/file3.fq"}},
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::count::check_filenames(filename_clusters, config));

    EXPECT_FALSE(config.precomputed_files);
}

TEST(check_filenames_test, minimizer_filenames)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters
    {
        {"key1", {"/path/to/file1.minimizer"}},
        {"key2", {"/path/to/file2.minimizer", "/path/to/file3.minimizer"}},
    };

    chopper::configuration config;

    EXPECT_NO_THROW(chopper::count::check_filenames(filename_clusters, config));

    EXPECT_TRUE(config.precomputed_files);
}

TEST(check_filenames_test, mixed_filenames_sequence_files)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters
    {
        {"key2", {"/path/to/file2.fa"}},
        {"key1", {"/path/to/file1.minimizer"}},
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::count::check_filenames(filename_clusters, config), std::invalid_argument);
}

TEST(check_filenames_test, mixed_filenames_minimizer_files)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters
    {
        {"key1", {"/path/to/file1.minimizer"}},
        {"key2", {"f.fa"}},
    };

    chopper::configuration config;

    EXPECT_THROW(chopper::count::check_filenames(filename_clusters, config), std::invalid_argument);
}
