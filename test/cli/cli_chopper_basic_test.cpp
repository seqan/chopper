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
        "chopper\n"
        "=======\n"
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
        "[CHOPPER LAYOUT ERROR] Option --tmax is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_empty_file)
{
    seqan3::test::tmp_filename empty_file_prefix{"empty"};

    {
        std::ofstream ofs{empty_file_prefix.get_path().string() + ".count"}; // opens file, s.t. it exists but is empty
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax", "64", /* required option */
                                         "--input-file", empty_file_prefix.get_path().c_str());

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] The file " + empty_file_prefix.get_path().string() +  ".count appears to be empty.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_no_extra_information)
{
    seqan3::test::tmp_filename prefix{"no_extra_information"};

    {
        std::ofstream ofs{prefix.get_path().string() + ".count"};
        ofs << "seq1\t500\n"
            << "seq2\t600\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax", "64",
                                         "--input-file", prefix.get_path(),
                                         "--aggregate-by-column", "3");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] Aggregate Error: You want to aggregate by something but your "
        "file does not contain any extra information columns.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_column_index_out_of_bounds)
{
    seqan3::test::tmp_filename prefix{"no_extra_information"};

    {
        std::ofstream ofs{prefix.get_path().string() + ".count"};
        ofs << "seq1\t500\tinformation1\n"
            << "seq2\t600\tinformation1\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax", "64",
                                         "--input-file", prefix.get_path(),
                                         "--aggregate-by-column", "4");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] Aggregate Error: You want to aggregate by a column index that is "
        "larger than the number of extra information columns.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_cmd_error_no_hll_dir)
{
    seqan3::test::tmp_filename prefix{"foo"};

    {
        std::ofstream ofs{prefix.get_path().string() + ".count"};
        ofs << "seq1\t500\n"
            << "seq2\t600\n";
    }

    cli_test_result result = execute_app("chopper",
                                         "--tmax", "64",
                                         "--input-file", prefix.get_path(),
                                         "--estimate-union");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] The directory " + prefix.get_path().string() + "_sketches must be present and not "
        "empty in order to enable --estimate-union or --rearrange-user-bins (created with chopper count).\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
