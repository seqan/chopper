#include <gtest/gtest.h>

#include <fstream>

#include <chopper/split/filename_batches_range.hpp>
#include <chopper/split/split_config.hpp>

#include "../api_test.hpp"

TEST(filename_batches_range_test, high_level_data_file)
{
    seqan3::test::tmp_filename binning_filename{"binning.out"};

    split_config config;
    config.out_path = "/dummy/"; // since no files are written
    config.data_filename = binning_filename.get_path().string();

    {
        std::ofstream fout{binning_filename.get_path()};
        fout << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << "seq7\t0\t2\t500\n"
             << "seq6\t2\t1\t500\n"
             << "seq0\t3;0\t1;16\t100;32\n"
             << "seq2\t3;16\t1;12\t100;42\n"
             << "seq3.1;seq3.2;seq3.3\t3;28\t1;12\t100;42\n"
             << "seq4\t4;0\t1;12\t100;42\n"
             << "seq5\t4;12\t1;12\t100;42\n"
             << "seq1.1;seq1.2\t5\t2\t1000\n";
    }

    filename_batches_range r{config};
    EXPECT_TRUE(r.current_file_type == filename_batches_range::file_type::data_file);

    std::vector<std::vector<std::string>> expected_seqfiles_range
    {
        {"seq7"}, {"seq6"}, {"seq0"}, {"seq2"}, {"seq3.1", "seq3.2", "seq3.3"}, {"seq4"}, {"seq5"}, {"seq1.1", "seq1.2"}
    };

    std::vector<int> expected_bins_range{2, 1, 16, 12, 12, 12, 12, 2};

    std::vector<bool> expected_merged{false, false, true, true, true, true, true, false};

    std::vector<std::vector<size_t>> expected_bin_indices
    {
        {0}, {2}, {3, 0}, {3, 16}, {3, 28}, {4, 0}, {4, 12}, {5}
    };

    auto it = r.begin();
    for (size_t i = 0; i < expected_seqfiles_range.size(); ++i, ++it)
    {
        EXPECT_RANGE_EQ((*it).seqfiles, expected_seqfiles_range[i]);
        EXPECT_EQ((*it).bins, expected_bins_range[i]) << " failed at " << i << std::endl;
        EXPECT_EQ((*it).merged_bin, expected_merged[i]) << " failed at " << i << std::endl;
        EXPECT_RANGE_EQ((*it).bin_indices, expected_bin_indices[i]);
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
