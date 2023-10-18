#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <sstream> // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>  // for allocator, string
#include <vector>  // for vector

#include <chopper/layout/input.hpp>
#include <chopper/layout/output.hpp>

TEST(output, user_bins)
{
    std::vector<std::string> const filenames{"file1.fa", "file2.fa", "path/to/file3.fa", "file4.fastq"};

    std::stringstream ss{};
    chopper::layout::write_user_bins_to(filenames, ss);

    std::string const expected{"@CHOPPER_USER_BINS\n"
                               "@0 file1.fa\n"
                               "@1 file2.fa\n"
                               "@2 path/to/file3.fa\n"
                               "@3 file4.fastq\n"
                               "@CHOPPER_USER_BINS_END\n"};

    EXPECT_EQ(ss.str(), expected);
}

TEST(input, user_bins)
{
    std::stringstream ss{"@CHOPPER_USER_BINS\n"
                         "@0 file1.fa\n"
                         "@1 file2.fa\n"
                         "@2 path/to/file3.fa\n"
                         "@3 file4.fastq\n"
                         "@CHOPPER_USER_BINS_END\n"};

    std::vector<std::vector<std::string>> filenames = chopper::layout::read_filenames_from(ss);
    std::vector<std::vector<std::string>> const expected{{"file1.fa"},
                                                         {"file2.fa"},
                                                         {"path/to/file3.fa"},
                                                         {"file4.fastq"}};

    EXPECT_EQ(filenames, expected);
}
