#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/build/read_chopper_pack_file.hpp>

TEST(read_chopper_pack_file_test, small_example)
{
    seqan3::test::tmp_filename chopper_split_filename{"test.split"};

    { // write example file
        std::ofstream fout{chopper_split_filename.get_path()};
        fout << "#MERGED_BIN_3 max_bin_id:16\n"
             << "#MERGED_BIN_6 max_bin_id:0\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_3\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             << "SPLIT_BIN_0\tfile7\t2\t500\n"
             << "SPLIT_BIN_2\tfile6\t1\t500\n"
             << "MERGED_BIN_3_0\tfile0\t16\t32\n"
             << "MERGED_BIN_3_16\tfile2\t12\t42\n"
             << "MERGED_BIN_3_28\tfile3.1;file3.2;file3.3\t12\t42\n"
             << "SPLIT_BIN_4\tfile1.1;file1.2\t2\t1000\n"
             << "MERGED_BIN_6_0\tfile4\t12\t42\n"
             << "MERGED_BIN_6_12\tfile5\t12\t42\n";
    }

    auto const && [data, records_per_hibf_bin] = read_chopper_pack_file(chopper_split_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 3);
    EXPECT_EQ(data.merged_max_bin_map.size(), 2);

    EXPECT_EQ(data.merged_max_bin_map.at(3), 16);
    EXPECT_EQ(data.merged_max_bin_map.at(6), 0);

    EXPECT_EQ(data.hibf_num_technical_bins, 7);
    EXPECT_EQ(data.hibf_max_bin, 3);

    // check records
    std::vector<std::vector<chopper_pack_record>> expected_records_per_hibf_bin
    {
        /*BIN 0*/ {{"SPLIT_BIN_0", 0, -1, {"file7"}, 2, 500}},
        /*BIN 1*/ {}, // belongs to SPLIT_BIN_0 since it is split into 2 bins
        /*BIN 2*/ {{"SPLIT_BIN_2", 2, -1, {"file6"}, 1, 500}},
        /*BIN 3*/ {{"MERGED_BIN_3_0", 3, 0, {"file0"}, 16, 32},
                   {"MERGED_BIN_3_16", 3, 16, {"file2"}, 12, 42},
                   {"MERGED_BIN_3_28", 3, 28, {"file3.1", "file3.2", "file3.3"}, 12, 42}},
        /*BIN 4*/ {{"SPLIT_BIN_4", 4, -1, {"file1.1", "file1.2"}, 2, 1000}},
        /*BIN 5*/ {}, // belongs to SPLIT_BIN_4 since it is split into 2 bins
        /*BIN 6*/ {{"MERGED_BIN_6_0", 6, 0, {"file4"}, 12, 42},
                   {"MERGED_BIN_6_12", 6, 12, {"file5"}, 12, 42}}
    };

    for (size_t i = 0; i < expected_records_per_hibf_bin.size(); ++i)
        EXPECT_RANGE_EQ(records_per_hibf_bin[i], expected_records_per_hibf_bin[i]);
}

TEST(read_chopper_pack_file_test, split_bin_with_multiple_bins_at_the_end)
{
    seqan3::test::tmp_filename chopper_split_filename{"test.split"};

    { // write example file
        std::ofstream fout{chopper_split_filename.get_path()};
        fout << "#MERGED_BIN_3 max_bin_id:16\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_3\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             << "SPLIT_BIN_0\tfile7\t2\t500\n"
             << "SPLIT_BIN_2\tfile6\t1\t500\n"
             << "MERGED_BIN_3_0\tfile0\t16\t32\n"
             << "MERGED_BIN_3_16\tfile2\t12\t42\n"
             << "MERGED_BIN_3_28\tfile3.1;file3.2;file3.3\t12\t42\n"
             << "SPLIT_BIN_4\tfile1.1;file1.2\t2\t1000\n";
    }

    auto const && [data, records_per_hibf_bin] = read_chopper_pack_file(chopper_split_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 3);
    EXPECT_EQ(data.merged_max_bin_map.size(), 1);

    EXPECT_EQ(data.merged_max_bin_map.at(3), 16);

    EXPECT_EQ(data.hibf_num_technical_bins, 6);

    // check records
    std::vector<std::vector<chopper_pack_record>> expected_records_per_hibf_bin
    {
        /*BIN 0*/ {{"SPLIT_BIN_0", 0, -1, {"file7"}, 2, 500}},
        /*BIN 1*/ {}, // belongs to SPLIT_BIN_0 since it is split into 2 bins
        /*BIN 2*/ {{"SPLIT_BIN_2", 2, -1, {"file6"}, 1, 500}},
        /*BIN 3*/ {{"MERGED_BIN_3_0", 3, 0, {"file0"}, 16, 32},
                   {"MERGED_BIN_3_16", 3, 16, {"file2"}, 12, 42},
                   {"MERGED_BIN_3_28", 3, 28, {"file3.1", "file3.2", "file3.3"}, 12, 42}},
        /*BIN 4*/ {{"SPLIT_BIN_4", 4, -1, {"file1.1", "file1.2"}, 2, 1000}},
        /*BIN 5*/ {}, // belongs to SPLIT_BIN_4 since it is split into 2 bins
    };

    for (size_t i = 0; i < expected_records_per_hibf_bin.size(); ++i)
        EXPECT_RANGE_EQ(records_per_hibf_bin[i], expected_records_per_hibf_bin[i]);
}
