#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <chopper/detail_parse_chopper_pack_header_line.hpp>

TEST(parse_chopper_pack_header_line_test, max_hibf_bin_line)
{
    std::string const line{"#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2"};

    build_data data;
    parse_chopper_pack_header_line(line, data);

    EXPECT_EQ(data.hibf_max_bin, 2);
}

TEST(parse_chopper_pack_header_line_test, max_libf_bin_line)
{
    std::string const line1{"#MERGED_BIN_2 max_bin_id:16"};
    std::string const line2{"#MERGED_BIN_24 max_bin_id:160"};

    build_data data;
    parse_chopper_pack_header_line(line1, data);
    parse_chopper_pack_header_line(line2, data);

    EXPECT_EQ(data.merged_max_bin_map.size(), 2);
    EXPECT_EQ(data.merged_max_bin_map.at(2), 16);
    EXPECT_EQ(data.merged_max_bin_map.at(24), 160);
}
