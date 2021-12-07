#include <gtest/gtest.h>

#include <chopper/detail_apply_prefix.hpp>

#include "../api_test.hpp"

TEST(detail_apply_prefix_test, simple_prefix)
{
    std::string const prefix{"foo"};
    std::filesystem::path filename;
    std::filesystem::path dir;

    chopper::detail::apply_prefix(prefix, filename, dir);

    EXPECT_EQ(filename, std::filesystem::path{"foo.count"});
    EXPECT_EQ(dir, std::filesystem::path{"foo_sketches"});
}

TEST(detail_apply_prefix_test, prefix_with_trailing_slash)
{
    std::string const prefix{"path/"};
    std::filesystem::path filename;
    std::filesystem::path dir;

    chopper::detail::apply_prefix(prefix, filename, dir);

    EXPECT_EQ(filename, std::filesystem::path{"path/chopper.count"});
    EXPECT_EQ(dir, std::filesystem::path{"path/chopper_sketches"});
}

TEST(detail_apply_prefix_test, prefix_with_file_extension)
{
    std::string const prefix{"foo.txt"};
    std::filesystem::path filename;
    std::filesystem::path dir;

    chopper::detail::apply_prefix(prefix, filename, dir);

    EXPECT_EQ(filename, std::filesystem::path{"foo.count"});
    EXPECT_EQ(dir, std::filesystem::path{"foo_sketches"});
}

TEST(detail_apply_prefix_test, prefix_with_absolut_path)
{
    std::string const prefix{"/path/to/foo.txt"};
    std::filesystem::path filename;
    std::filesystem::path dir;

    chopper::detail::apply_prefix(prefix, filename, dir);

    EXPECT_EQ(filename, std::filesystem::path{"/path/to/foo.count"});
    EXPECT_EQ(dir, std::filesystem::path{"/path/to/foo_sketches"});
}
