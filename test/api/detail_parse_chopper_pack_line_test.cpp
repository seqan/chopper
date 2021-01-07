#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/detail_parse_chopper_pack_line.hpp>

TEST(chopper_pack_record_test, euality_operator)
{
    chopper_pack_record r1{"SPLIT_BIN_0", 0, -1, {"file7"}, 2, 500};
    chopper_pack_record r2{"SPLIT_BIN_0", 0, -1, {"file7"}, 2, 500};

    chopper_pack_record r3{"SPLIT_BIN_0", 0, -1, {"FOO"}, 2, 500};

    EXPECT_EQ(r1, r2);
    EXPECT_NE(r1, r3);
    EXPECT_NE(r2, r3);
}

TEST(parse_chopper_pack_line_test, high_level_data_file)
{
    std::vector<std::string> const lines
    {
        "SPLIT_BIN_0\tseq7\t2\t500\n",
        "SPLIT_BIN_1\tseq6\t1\t500\n",
        "MERGED_BIN_2_0\tseq0\t16\t32\n",
        "MERGED_BIN_2_16\tseq2\t12\t42\n",
        "MERGED_BIN_2_28\tseq3.1;seq3.2;seq3.3\t12\t42\n",
        "MERGED_BIN_20_0\tseq4\t12\t42\n",
        "MERGED_BIN_20_12\tseq5\t12\t42\n",
        "SPLIT_BIN_3\tseq1.1;seq1.2\t2\t1000\n"
    };

    std::vector<std::vector<std::string>> expected_filenames
    {
        {"seq7"}, {"seq6"}, {"seq0"}, {"seq2"}, {"seq3.1", "seq3.2", "seq3.3"}, {"seq4"}, {"seq5"}, {"seq1.1", "seq1.2"}
    };

    std::vector<int> expected_bins{2, 1, 16, 12, 12, 12, 12, 2};

    std::vector<std::string> expected_bin_names
    {
        "SPLIT_BIN_0",
        "SPLIT_BIN_1",
        "MERGED_BIN_2_0",
        "MERGED_BIN_2_16",
        "MERGED_BIN_2_28",
        "MERGED_BIN_20_0",
        "MERGED_BIN_20_12",
        "SPLIT_BIN_3"
    };

    std::vector<int> expected_max_size{500, 500, 32, 42, 42, 42, 42, 1000};

    std::vector<int64_t> expected_hidxs{0, 1, 2, 2, 2, 20, 20, 3};

    std::vector<int64_t> expected_lidxs{-1, -1, 0, 16, 28, 0, 12, -1};

    for (size_t i = 0; i < expected_filenames.size(); ++i)
    {
        chopper_pack_record && record = parse_chopper_pack_line(lines[i]);

        EXPECT_EQ(record.bin_name, expected_bin_names[i]);
        EXPECT_RANGE_EQ(record.filenames, expected_filenames[i]);
        EXPECT_EQ(record.bins, expected_bins[i]) << " failed at " << expected_bin_names[i] << std::endl;
        EXPECT_EQ(record.max_size, expected_max_size[i]) << " failed at " << expected_bin_names[i] << std::endl;
        EXPECT_EQ(record.hidx, expected_hidxs[i]) << " failed at " << expected_bin_names[i] << std::endl;
        EXPECT_EQ(record.lidx, expected_lidxs[i]) << " failed at " << expected_bin_names[i] << std::endl;
    }
}
