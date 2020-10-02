#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/split/split_config.hpp>
#include <chopper/split/filename_batches_range.hpp>

TEST(filename_batches_range_test, high_level_data_file)
{
    split_config config;
    seqan3::test::tmp_filename binning_filename{"binning.out"};

    {
        std::ofstream fout{binning_filename.get_path()};
        fout << "BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
             << "SPLIT_BIN_0\tseq7\t2\t500\n"
             << "SPLIT_BIN_1\tseq6\t1\t500\n"
             << "COLORFUL_MERGED_BIN_2_0\tseq0\t16\t32\n"
             << "COLORFUL_MERGED_BIN_2_1\tseq2\t12\t42\n"
             << "COLORFUL_MERGED_BIN_2_2\tseq3.1;seq3.2;seq3.3\t12\t42\n"
             << "COLORFUL_MERGED_BIN_2_3\tseq4\t12\t42\n"
             << "COLORFUL_MERGED_BIN_2_4\tseq5\t12\t42\n"
             << "SPLIT_BIN_3\tseq1.1;seq1.2\t2\t1000\n";
    }

    config.data_filename = binning_filename.get_path().string();

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::data_file);

    std::vector<std::vector<std::string>> expected_seqfiles_range
    {
        {"seq7"}, {"seq0"}, {"seq2"}, {"seq3.1", "seq3.2", "seq3.3"}, {"seq4"}, {"seq5"}, {"seq1.1", "seq1.2"}
    };

    std::vector<int> expected_bins_range{2, 16, 12, 12, 12, 12, 2};

    auto it = r.begin();
    for (size_t i = 0; i < expected_seqfiles_range.size(); ++i, ++it)
    {
        EXPECT_RANGE_EQ((*it).seqfiles, expected_seqfiles_range[i]);
        EXPECT_EQ((*it).bins, expected_bins_range[i]);
    }
    EXPECT_TRUE(it == r.end());
}

TEST(filename_batches_range_test, no_data_file)
{
    split_config config;
    config.seqfiles = {"seq3.1", "seq3.2", "seq3.3"};

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::seqfiles_given);

    auto it = r.begin();

    EXPECT_RANGE_EQ((*it).seqfiles, config.seqfiles);
    EXPECT_EQ((*it).bins, config.bins);

    ++it;

    EXPECT_TRUE(it == r.end());
}
