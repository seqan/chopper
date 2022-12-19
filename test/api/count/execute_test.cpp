#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/utility/range/to.hpp>

#include "../api_test.hpp"
#include <chopper/count/execute.hpp>
#include <chopper/detail_apply_prefix.hpp>

TEST(execute_test, small_example_parallel_2_threads)
{
    std::string input_filename = data("small.fa");
    seqan3::test::tmp_filename data_filename{"data.tsv"};
    seqan3::test::tmp_filename output_prefix{"small_example"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
            /* << input_filename << '\t' << "TAX2\n" */;
    }

    chopper::configuration config{};
    config.threads = 2;
    config.k = 15;
    config.column_index_to_cluster = 2;
    config.disable_sketch_output = true;
    config.data_file = data_filename.get_path();
    config.output_prefix = output_prefix.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::data_store store{};

    EXPECT_NO_THROW(chopper::count::execute(config, store));

    EXPECT_RANGE_EQ(store.filenames, (std::vector<std::string>{input_filename, input_filename}));
    ASSERT_EQ(store.sketches.size(), 2);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);
    EXPECT_EQ(std::lround(store.sketches[1].estimate()), 571);
}

TEST(execute_test, some_test)
{
    std::string input_filename = data("small.fa");
    seqan3::test::tmp_filename data_filename{"data.tsv"};
    seqan3::test::tmp_filename output_prefix{"small_example"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
            /* << input_filename << '\t' << "TAX2\n" */;
    }

    chopper::configuration config{};
    config.threads = 1;
    config.k = 25;
    config.column_index_to_cluster = 2;
    config.disable_sketch_output = true;
    config.data_file = data_filename.get_path();
    config.output_prefix = output_prefix.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::data_store store{};

    chopper::count::execute(config, store);

    EXPECT_RANGE_EQ(store.filenames, (std::vector<std::string>{input_filename, input_filename}));
    ASSERT_EQ(store.sketches.size(), 2);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 591);
    EXPECT_EQ(std::lround(store.sketches[1].estimate()), 591);
}
