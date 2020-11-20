#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/detail_parse_binning_line.hpp>

TEST(filename_batches_range_test, high_level_data_file)
{
    std::vector<std::string> const lines
    {
        "SPLIT_BIN_0\tseq7\t2\t500\n",
        "SPLIT_BIN_1\tseq6\t1\t500\n",
        "MERGED_BIN_2_0\tseq0\t16\t32\n",
        "MERGED_BIN_2_1\tseq2\t12\t42\n",
        "MERGED_BIN_2_2\tseq3.1;seq3.2;seq3.3\t12\t42\n",
        "MERGED_BIN_20_3\tseq4\t12\t42\n",
        "MERGED_BIN_20_4\tseq5\t12\t42\n",
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
        "MERGED_BIN_2_1",
        "MERGED_BIN_2_2",
        "MERGED_BIN_20_3",
        "MERGED_BIN_20_4",
        "SPLIT_BIN_3"
    };

    std::vector<int> expected_max_size{500, 500, 32, 42, 42, 42, 42, 1000};

    for (size_t i = 0; i < expected_filenames.size(); ++i)
    {
        data_file_record && record = parse_binning_line(lines[i]);

        EXPECT_EQ(record.bin_name, expected_bin_names[i]);
        EXPECT_RANGE_EQ(record.filenames, expected_filenames[i]);
        EXPECT_EQ(record.bins, expected_bins[i]) << " failed at " << expected_bin_names[i] << std::endl;
        EXPECT_EQ(record.max_size, expected_max_size[i]) << " failed at " << expected_bin_names[i] << std::endl;
    }
}
