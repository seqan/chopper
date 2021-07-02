#include <gtest/gtest.h>

#include <chopper/count/count_config.hpp>
#include <chopper/count/count_kmers.hpp>

#include "../api_test.hpp"

TEST(count_kmers_test, small_example_serial)
{
    seqan3::test::tmp_filename output_filename{"kmer_counts.txt"};

    count_config config;
    config.k = 15;
    config.w = 25;
    config.num_threads = 1;
    config.output_filename = output_filename.get_path();

    std::string input_file = data("small.fa");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}},
        {"TAX2", {input_file, input_file}}
    };

    std::string expected
    {
        input_file + "\t95\tTAX1\n" +
        input_file + ";" + input_file + "\t95\tTAX2\n"
    };

    count_kmers(filename_clusters, config);

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}

TEST(count_kmers_test, small_example_hll)
{
    seqan3::test::tmp_filename output_filename{"kmer_counts.txt"};
    seqan3::test::tmp_filename hll_dir{"hll"};

    count_config config;
    config.k = 15;
    config.w = 25;
    config.num_threads = 1;
    config.exclusively_hlls = true;
    config.hll_dir = hll_dir.get_path();
    config.output_filename = output_filename.get_path();

    std::string input_file = data("small.fa");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}}
    };

    std::string expected
    {
        input_file + "\t96\tTAX1\n"
    };

    count_kmers(filename_clusters, config);

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}

TEST(count_kmers_test, small_example_parallel_2_threads)
{
    seqan3::test::tmp_filename output_filename{"kmer_counts.txt"};

    count_config config;
    config.k = 15;
    config.w = 25;
    config.num_threads = 2;
    config.output_filename = output_filename.get_path();

    std::string input_file = data("small.fa");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}},
        {"TAX2", {input_file, input_file}}
    };

    count_kmers(filename_clusters, config);

    std::vector<std::string> expected_components
    {
        input_file + "\t95\tTAX1",
        input_file + ";" + input_file + "\t95\tTAX2"
    };

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : output_file_str | std::views::split('\n') | seqan3::views::to<std::vector<std::string>>)
    {
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end());
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);
}
