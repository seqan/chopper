#include <gtest/gtest.h>

#include <vector>

#include <chopper/sketch/execute.hpp>

#include "../api_test.hpp"

TEST(execute_test, small_example_parallel_2_threads)
{
    std::string input_filename = data("small.fa");

    chopper::configuration config{};
    config.threads = 2;
    config.k = 15;
    config.disable_sketch_output = true;

    std::vector<std::string> filenames{input_filename, input_filename};
    std::vector<chopper::sketch::hyperloglog> sketches{};

    EXPECT_NO_THROW(chopper::sketch::execute(config, filenames, sketches));

    ASSERT_EQ(sketches.size(), 2);
    EXPECT_EQ(std::lround(sketches[0].estimate()), 571);
    EXPECT_EQ(std::lround(sketches[1].estimate()), 571);
}

TEST(execute_test, some_test)
{
    std::string input_filename = data("small.fa");

    chopper::configuration config{};
    config.threads = 1;
    config.k = 25;
    config.disable_sketch_output = true;

    std::vector<std::string> filenames{input_filename, input_filename};
    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::execute(config, filenames, sketches);

    ASSERT_EQ(sketches.size(), 2);
    EXPECT_EQ(std::lround(sketches[0].estimate()), 591);
    EXPECT_EQ(std::lround(sketches[1].estimate()), 591);
}
