#include <gtest/gtest.h>

#include <chopper/configuration.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/detail_apply_prefix.hpp>

#include "../api_test.hpp"

TEST(count_kmers_test, small_example)
{
    seqan3::test::tmp_filename output_prefix{"small"};

    chopper::configuration config;
    config.k = 15;
    config.threads = 1;
    config.output_prefix = output_prefix.get_path().string();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    std::string input_file = data("small.fa");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}}
    };

    std::string expected
    {
        input_file + "\t571\tTAX1\n"
    };

    chopper::count::count_kmers(filename_clusters, config);

    ASSERT_TRUE(std::filesystem::exists(config.count_filename));
    std::ifstream output_file{config.count_filename};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}

TEST(count_kmers_test, small_example_parallel_2_threads)
{
    seqan3::test::tmp_filename output_prefix{"parallel"};

    chopper::configuration config;
    config.k = 15;
    config.threads = 2;
    config.output_prefix = output_prefix.get_path().string();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    std::string input_file = data("small.fa");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}},
        {"TAX2", {input_file/* , input_file -- not supported by hll sketches yet*/}}
    };

    chopper::count::count_kmers(filename_clusters, config);

    std::vector<std::string> expected_components
    {
        input_file + "\t571\tTAX1",
        input_file + /* ";" + input_file + */ "\t571\tTAX2"
    };

    ASSERT_TRUE(std::filesystem::exists(config.count_filename));
    std::ifstream output_file{config.count_filename};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : output_file_str | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "missing: " << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count);
}

TEST(count_kmers_test, read_in_precomputed_binary_files)
{
    seqan3::test::tmp_filename output_prefix{"small"};

    chopper::configuration config;
    config.k = 15;
    config.threads = 1;
    config.output_prefix = output_prefix.get_path().string();
    config.precomputed_files = true;
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    std::string input_file = data("small.minimizer");

    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters
    {
        {"TAX1", {input_file}}
    };

    std::string expected
    {
        input_file + "\t571\tTAX1\n"
    };

    chopper::count::count_kmers(filename_clusters, config);

    ASSERT_TRUE(std::filesystem::exists(config.count_filename));
    std::ifstream output_file{config.count_filename};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}
