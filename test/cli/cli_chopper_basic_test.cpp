#include <fstream>
#include <ranges> // range comparisons
#include <string> // strings
#include <vector> // vectors

#include <seqan3/test/tmp_directory.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("chopper");
    std::string expected{"chopper - Compute an HIBF layout\n"
                         "================================\n"
                         "    Try -h or --help for more information.\n"};
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, chopper_cmd_error_unknown_option)
{
    cli_test_result result = execute_app("chopper", "--unkown-option");
    std::string expected{"[CHOPPER ERROR] Option --input-file is required but not set.\n"};
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_empty_file)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path empty_file{tmp_dir.path() / "empty.count"};

    {
        std::ofstream ofs{empty_file.string()}; // opens file, s.t. it exists but is empty
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax",
                                         "64", /* required option */
                                         "--input-file",
                                         empty_file.c_str());

    std::string expected{"[CHOPPER ERROR] The file " + empty_file.string() + " appears to be empty.\n"};
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
