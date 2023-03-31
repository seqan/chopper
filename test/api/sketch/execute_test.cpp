#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/utility/range/to.hpp>

#include <chopper/sketch/execute.hpp>
#include <chopper/detail_apply_prefix.hpp>

#include "../api_test.hpp"

TEST(execute_test, small_example_parallel_2_threads)
{
    std::string input_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path data_filename{tmp_dir.path() / "data.tsv"};

    // generate data filename
    {
        std::ofstream fout{data_filename};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
            /* << input_filename << '\t' << "TAX2\n" */;
    }

    chopper::configuration config{};
    config.threads = 2;
    config.k = 15;
    config.column_index_to_cluster = 2;
    config.disable_sketch_output = true;
    config.data_file = data_filename;

    chopper::data_store store{};

    EXPECT_NO_THROW(chopper::sketch::execute(config, store));

    EXPECT_RANGE_EQ(store.filenames, (std::vector<std::string>{input_filename, input_filename}));
    ASSERT_EQ(store.sketches.size(), 2);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 571);
    EXPECT_EQ(std::lround(store.sketches[1].estimate()), 571);
}

TEST(execute_test, some_test)
{
    std::string input_filename = data("small.fa");
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path data_filename{tmp_dir.path() / "data.tsv"};

    // generate data filename
    {
        std::ofstream fout{data_filename};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
            /* << input_filename << '\t' << "TAX2\n" */;
    }

    chopper::configuration config{};
    config.threads = 1;
    config.k = 25;
    config.column_index_to_cluster = 2;
    config.disable_sketch_output = true;
    config.data_file = data_filename;

    chopper::data_store store{};

    chopper::sketch::execute(config, store);

    EXPECT_RANGE_EQ(store.filenames, (std::vector<std::string>{input_filename, input_filename}));
    ASSERT_EQ(store.sketches.size(), 2);
    EXPECT_EQ(std::lround(store.sketches[0].estimate()), 591);
    EXPECT_EQ(std::lround(store.sketches[1].estimate()), 591);
}
