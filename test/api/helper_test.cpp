#include <gtest/gtest.h>

#include <chopper/helper.hpp>

#include "api_test.hpp"

TEST(byte_size_to_formatted_str, correct_output)
{
    EXPECT_EQ(chopper::byte_size_to_formatted_str(size_t{55}), std::string{"55 Bytes"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(size_t{55} * 1024), std::string{"55 KiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(size_t{55} * 1024 * 1024), std::string{"55 MiB"});
    EXPECT_EQ(chopper::byte_size_to_formatted_str(size_t{55} * 1024 * 1024 * 1024), std::string{"55 GiB"});
}
