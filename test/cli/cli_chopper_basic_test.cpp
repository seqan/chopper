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

TEST_F(cli_test, chopper_cmd_error)
{
    cli_test_result result = execute_app("chopper", "nonexistingsubmodule");
    std::string expected
    {
        "[CHOPPER ERROR] You either forgot or misspelled the subcommand! "
        "Please specify which sub-program you want to use: one of [count,layout]. "
        "Use -h/--help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_count_cmd_error_unknown_option)
{
    cli_test_result result = execute_app("chopper", "count", "--unkown-option");
    std::string expected
    {
        "[CHOPPER COUNT ERROR] Option -f/--data_file is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_layout_cmd_error_unknown_option)
{
    cli_test_result result = execute_app("chopper", "layout", "--unkown-option");
    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] Option -f/--filenames is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_layout_cmd_error_empty_file)
{
    seqan3::test::tmp_filename empty_file{"empty.file"};

    {
        std::ofstream ofs{empty_file.get_path()}; // opens file, s.t. it exists but is empty
    }

    cli_test_result result = execute_app("chopper", "layout", "-f", empty_file.get_path());

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] The file " + empty_file.get_path().string() +  " appears to be empty.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_layout_cmd_error_no_extra_information)
{
    seqan3::test::tmp_filename empty_file{"no_extra_information.data"};

    {
        std::ofstream ofs{empty_file.get_path()};
        ofs << "seq1\t500\n"
            << "seq2\t600\n";
    }

    cli_test_result result = execute_app("chopper", "layout", "-f", empty_file.get_path(), "--aggregate-by", "3");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] Aggregate Error: You want to aggregate by something but your "
        "file does not contain any extra information columns.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_layout_cmd_error_column_index_out_of_bounds)
{
    seqan3::test::tmp_filename empty_file{"no_extra_information.data"};

    {
        std::ofstream ofs{empty_file.get_path()};
        ofs << "seq1\t500\tinformation1\n"
            << "seq2\t600\tinformation1\n";
    }

    cli_test_result result = execute_app("chopper", "layout", "-f", empty_file.get_path(), "--aggregate-by", "4");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] Aggregate Error: You want to aggregate by a column index that is "
        "larger than the number of extra information columns.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, chopper_layout_cmd_error_no_hll_dir)
{
    seqan3::test::tmp_filename empty_file{"my.data"};

    {
        std::ofstream ofs{empty_file.get_path()};
        ofs << "seq1\t500\n"
            << "seq2\t600\n";
    }

    cli_test_result result = execute_app("chopper", "layout", "-f", empty_file.get_path(), "--estimate-union");

    std::string expected
    {
        "[CHOPPER LAYOUT ERROR] An hll dir needs to be provided when enabling -u or -r.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
