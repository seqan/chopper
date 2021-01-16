#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/detail_parse_chopper_pack_line.hpp>

TEST(chopper_pack_record_test, euality_operator)
{
    chopper_pack_record r1{{"file7"}, {0, 1}, {1, 2}, {500, 20}};
    chopper_pack_record r2{{"file7"}, {0, 1}, {1, 2}, {500, 20}};

    chopper_pack_record r3{{"FOO"}, {0, 1}, {1, 2}, {500, 20}};

    EXPECT_EQ(r1, r2);
    EXPECT_NE(r1, r3);
    EXPECT_NE(r2, r3);
}

TEST(parse_chopper_pack_line_test, high_level_data_file)
{
    std::vector<std::string> const lines
    {
        "seq7\t0\t2\t500\n",
        "seq6\t1\t1\t500\n",
        "seq0\t2;0\t1;16\t400;32\n",
        "seq2\t2;16\t1;12\t400;42\n",
        "seq3.1;seq3.2;seq3.3\t2;28\t1;12\t400;42\n",
        "seq4\t20;0\t1;12\t200;42\n",
        "seq5\t20;12\t1;12\t200;42\n",
        "seq1.1;seq1.2\t3\t2\t1000\n"
    };

    std::vector<std::vector<std::string>> expected_filenames
    {
        {"seq7"}, {"seq6"}, {"seq0"}, {"seq2"}, {"seq3.1", "seq3.2", "seq3.3"}, {"seq4"}, {"seq5"}, {"seq1.1", "seq1.2"}
    };

    std::vector<std::vector<size_t>> expected_bin_indices
    {
        {0}, {1}, {2, 0}, {2, 16}, {2, 28}, {20, 0}, {20, 12}, {3}
    };

    std::vector<std::vector<size_t>> expected_number_of_bins
    {
        {2}, {1}, {1, 16}, {1, 12}, {1, 12}, {1, 12}, {1, 12}, {2}
    };

    std::vector<std::vector<size_t>> expected_estimated_sizes
    {
        {500}, {500}, {400, 32}, {400, 42}, {400, 42}, {200, 42}, {200, 42}, {1000}
    };

    for (size_t i = 0; i < expected_filenames.size(); ++i)
    {
        chopper_pack_record && record = parse_chopper_pack_line(lines[i]);

        EXPECT_RANGE_EQ(record.filenames, expected_filenames[i]);
        EXPECT_RANGE_EQ(record.bin_indices, expected_bin_indices[i]);
        EXPECT_RANGE_EQ(record.number_of_bins, expected_number_of_bins[i]);
        EXPECT_RANGE_EQ(record.estimated_sizes, expected_estimated_sizes[i]);
    }
}
