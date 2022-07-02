#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/utility/range/to.hpp>

#include <chopper/count/execute.hpp>

#include "../api_test.hpp"

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

    const char * argv[] = {"./chopper-count",
                           "--kmer-size", "15",
                           "--threads", "2",
                           "--column-index", "2",
                           "--disable-sketch-output",
                           "--input-file", data_filename.get_path().c_str(),
                           "--output-prefix", output_prefix.get_path().c_str()};
    int argc = 12;
    sharg::parser count_parser{"chopper-count", argc, argv, sharg::update_notifications::off};

    std::vector<std::string> expected_components
    {
        input_filename + "\t571\tTAX1",
        input_filename + /* ";" + input_filename */ + "\t571\tTAX2"
    };

    EXPECT_NO_THROW(chopper::count::execute(count_parser));

    std::ifstream output_file{output_prefix.get_path().string() + ".count"};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : output_file_str | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "missing:" << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count) << "File: " <<  output_file_str;
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

    const char * argv[] = {"./chopper-count",
                           "--kmer-size", "25",
                           "--threads", "1",
                           "--column-index", "2",
                           "--disable-sketch-output",
                           "--input-file", data_filename.get_path().c_str(),
                           "--output-prefix", output_prefix.get_path().c_str()};
    int argc = 12;
    sharg::parser count_parser{"chopper-count", argc, argv, sharg::update_notifications::off};

    std::vector<std::string> expected_components
    {
        input_filename + "\t590\tTAX1",
        input_filename + /* ";" + input_filename */ + "\t590\tTAX2"
    };

    chopper::count::execute(count_parser);

    std::ifstream output_file{output_prefix.get_path().string() + ".count"};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    size_t line_count{};
    for (auto && line : output_file_str | std::views::split('\n') | seqan3::ranges::to<std::vector<std::string>>())
    {
#if defined(__GNUC__) && (__GNUC__ == 12)
        if (line.empty())
            continue;
#endif
        EXPECT_TRUE(std::ranges::find(expected_components, line) != expected_components.end()) << "missing:" << line;
        ++line_count;
    }

    EXPECT_EQ(expected_components.size(), line_count) << "File: " <<  output_file_str;
}
