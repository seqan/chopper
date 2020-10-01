#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("chopper");
    std::string expected
    {
        "chopper\n"
        "=======\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_out_file)
{
    seqan3::test::tmp_filename output_filename{"small_traverse.out"};

    cli_test_result result = execute_app("chopper", "split",
                                         "-k", "15", "-w", "25", "-b", "3",
                                         "-s", data("small.fa"),
                                         "-o", output_filename.get_path().c_str());

    // compare results
    std::ifstream expected_file{DATADIR"small_traverse.out"};
    std::ifstream output_file{output_filename.get_path()};

    std::string expected_line;
    std::string output_line;

    while (std::getline(expected_file, expected_line) && std::getline(output_file, output_line))
        EXPECT_EQ(expected_line, output_line);

    EXPECT_FALSE(std::getline(expected_file, expected_line)); // both files are exhausted
    EXPECT_FALSE(std::getline(output_file, output_line)); // both files are exhausted

    EXPECT_EQ(result.exit_code, 0);
}
