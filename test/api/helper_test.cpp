#include <gtest/gtest.h>

#include <chopper/helper.hpp>

TEST(byte_size_to_formatted_str, storage_unit)
{
    // `ULL` means unsigned long long; ensures that shifting works since `55` will usually be a 32bit-int
    // `<< 10` is the same as `* 1024`, `<< 20` the same as `* 1024 * 1024`, and so on
    EXPECT_EQ("8Bytes", chopper::byte_size_to_formatted_str(8ULL));
    EXPECT_EQ("8.0KiB", chopper::byte_size_to_formatted_str(8ULL << 10));
    EXPECT_EQ("8.0MiB", chopper::byte_size_to_formatted_str(8ULL << 20));
    EXPECT_EQ("8.0GiB", chopper::byte_size_to_formatted_str(8ULL << 30));
    EXPECT_EQ("8.0TiB", chopper::byte_size_to_formatted_str(8ULL << 40));
    EXPECT_EQ("8.0PiB", chopper::byte_size_to_formatted_str(8ULL << 50));
    EXPECT_EQ("8.0EiB", chopper::byte_size_to_formatted_str(8ULL << 60));
}

TEST(byte_size_to_formatted_str, rounding)
{
    EXPECT_EQ("5.8GiB", chopper::byte_size_to_formatted_str(6'174'015'488ULL));
    EXPECT_EQ("5.7GiB", chopper::byte_size_to_formatted_str(6'174'015'487ULL));
    // This is 1 bytes short of 1MiB. Rounding should change the unit from KiB to MiB.
    EXPECT_EQ("1.0MiB", chopper::byte_size_to_formatted_str(1'048'575ULL));
}

TEST(byte_size_to_formatted_str, edge_cases)
{
    EXPECT_EQ("0Bytes", chopper::byte_size_to_formatted_str(0ULL));
    EXPECT_EQ("16.0EiB", chopper::byte_size_to_formatted_str(std::numeric_limits<size_t>::max()));
}
