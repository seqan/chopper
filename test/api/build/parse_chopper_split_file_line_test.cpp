#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <chopper/build/parse_chopper_split_file_line.hpp>

TEST(parse_chopper_split_file_line_test, small_example)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\t5\n"};

    auto && [filename, id, reg] = parse_chopper_split_file_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.hidx, 10);
    EXPECT_EQ(reg.lidx, 5);
}

TEST(parse_chopper_split_file_line_test, low_level_idx_not_given)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\t-\n"};

    auto && [filename, id, reg] = parse_chopper_split_file_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.hidx, 10);
    EXPECT_EQ(reg.lidx, -1);
}
