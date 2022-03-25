#include <gtest/gtest.h>

#include <chopper/helper.hpp>

TEST(byte_size_to_formatted_str, correct_output)
{
    // `ULL` means unsigned long long; ensures that shifting works since `55` will usually be a 32bit-int
    // `<< 10` is the same as `* 1024`, `<< 20` the same as `* 1024 * 1024`, and so on
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL), std::string{"55Bytes"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL << 10), std::string{"55KiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL << 20), std::string{"55MiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL << 30), std::string{"55GiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL << 40), std::string{"55TiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(55ULL << 50), std::string{"55PiB"});
}
