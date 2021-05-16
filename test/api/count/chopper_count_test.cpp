#include <gtest/gtest.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/count/count_config.hpp>
#include <chopper/count/chopper_count.hpp>

TEST(chopper_count_test, small_example_parallel_2_threads)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
             << input_filename << '\t' << "TAX2\n";
    }

    const char * argv[] = {"./chopper-count", "-k", "15", "-w", "25", "-t", "2", "-c", "2",
                           "-f", data_filename.get_path().c_str()};
    int argc = 11;
    seqan3::argument_parser count_parser{"chopper-count", argc, argv, seqan3::update_notifications::off};

    std::string expected
    {
        input_filename + "\t95\tTAX1\n" +
        input_filename + ";" + input_filename + "\t95\tTAX2\n"
    };

    testing::internal::CaptureStdout();
    chopper_count(count_parser);
    std::string std_cout = testing::internal::GetCapturedStdout();

    EXPECT_EQ(expected, std_cout);
}

TEST(chopper_count_test, disable_minimizers)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
             << input_filename << '\t' << "TAX2\n";
    }

    const char * argv[] = {"./chopper-count",
                           "-k", "25",
                           "-t", "1",
                           "-c", "2",
                           "--disable-minimizers",
                           "-f", data_filename.get_path().c_str()};
    int argc = 10;
    seqan3::argument_parser count_parser{"chopper-count", argc, argv, seqan3::update_notifications::off};

    std::string expected
    {
        input_filename + "\t585\tTAX1\n" +
        input_filename + ";" + input_filename + "\t585\tTAX2\n"
    };

    testing::internal::CaptureStdout();
    chopper_count(count_parser);
    std::string std_cout = testing::internal::GetCapturedStdout();

    EXPECT_EQ(expected, std_cout);
}
