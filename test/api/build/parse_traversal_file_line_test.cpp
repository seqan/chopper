#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <chopper/build/parse_traversal_file_line.hpp>

TEST(parse_traversal_file_line_test, small_example)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\n"};

    auto && [filename, id, begin, end, idx] = parse_traversal_file_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_EQ(begin, 216);
    EXPECT_EQ(end, 400);
    EXPECT_EQ(idx, 10);
}
