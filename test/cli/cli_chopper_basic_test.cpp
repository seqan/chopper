#include <ranges>     // range comparisons
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
        "chopper - Compute an HIBF layout\n"
        "================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, chopper_cmd_error_unknown_option)
{
    cli_test_result result = execute_app("chopper", "--unkown-option");
    std::string expected
    {
        "[CHOPPER ERROR] Option --input-file is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_tmax_missing)
{
    cli_test_result result = execute_app("chopper", "--input-file", "foo.txt");
    std::string expected
    {
        "[CHOPPER ERROR] Option --tmax is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_empty_file)
{
    seqan3::test::tmp_filename empty_file{"empty.count"};

    {
        std::ofstream ofs{empty_file.get_path().string()}; // opens file, s.t. it exists but is empty
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax", "64", /* required option */
                                         "--input-file", empty_file.get_path().c_str());

    std::string expected
    {
        "terminate called after throwing an instance of 'seqan3::argument_parser_error'\n"
        "  what():  [CHOPPER ERROR] The file " + empty_file.get_path().string() +  " appears to be empty.\n"
        "Aborted (core dumped)\n"
    };
    EXPECT_EQ(result.exit_code, 34304);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
