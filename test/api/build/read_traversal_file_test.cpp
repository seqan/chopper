#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/build/read_traversal_file.hpp>

TEST(parse_traversal_header_line_test, max_hibf_bin_line)
{
    std::string const line{"#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2"};

    build_data data;
    parse_traversal_header_line(line, data);

    EXPECT_EQ(data.hibf_max_bin, 2);
}

TEST(parse_traversal_header_line_test, max_libf_bin_line)
{
    std::string const line1{"#MERGED_BIN_2 max_bin_id:16"};
    std::string const line2{"#MERGED_BIN_24 max_bin_id:160"};

    build_data data;
    parse_traversal_header_line(line1, data);
    parse_traversal_header_line(line2, data);

    EXPECT_EQ(data.merged_max_bin_map.size(), 2);
    EXPECT_EQ(data.merged_max_bin_map.at(2), 16);
    EXPECT_EQ(data.merged_max_bin_map.at(24), 160);
}

TEST(parse_traversal_line_test, small_example)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\t5\n"};

    auto && [filename, id, reg] = parse_traversal_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.hidx, 10);
    EXPECT_EQ(reg.lidx, 5);
}

TEST(parse_traversal_line_test, low_level_idx_not_given)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\t-\n"};

    auto && [filename, id, reg] = parse_traversal_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.hidx, 10);
    EXPECT_EQ(reg.lidx, -1);
}

TEST(read_traversal_file_test, small_example)
{
    std::string seq_filename = "test.fa";
    std::string seq_filename2 = "test1.fa";
    std::string seq_filename3 = "test2.fa";
    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    { // write example file
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_2 max_bin_id:1\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*SPLIT_BIN_0*/
             << seq_filename + "\tseq1\t0\t209\t0\t-\n"
             << seq_filename + "\tseq2\t0\t289\t0\t-\n"
             << seq_filename + "\tseq3\t0\t209\t0\t-\n"
             << seq_filename + "\tseq1\t209\t400\t1\t-\n"
             << seq_filename + "\tseq2\t289\t480\t1\t-\n"
             << seq_filename + "\tseq3\t209\t481\t1\t-\n"
             /*MERGED_BIN_2_0*/
             << seq_filename2 + "\tseq10\t0\t209\t2\t0\n"
             << seq_filename2 + "\tseq20\t0\t289\t2\t0\n"
             << seq_filename2 + "\tseq30\t0\t209\t2\t0\n"
             << seq_filename2 + "\tseq10\t209\t400\t2\t1\n"
             << seq_filename2 + "\tseq20\t289\t480\t2\t1\n"
             << seq_filename2 + "\tseq30\t209\t481\t2\t1\n"
             /*MERGED_BIN_2_1*/
             << seq_filename2 + "\tseq10\t0\t400\t2\t2\n"
             << seq_filename2 + "\tseq20\t0\t480\t2\t2\n"
             << seq_filename2 + "\tseq30\t0\t481\t2\t2\n"
             /*SPLIT_BIN_3*/
             << seq_filename3 + "\tseq1\t0\t163\t3\t-\n"
             << seq_filename3 + "\tseq2\t0\t186\t3\t-\n"
             << seq_filename3 + "\tseq3\t0\t163\t3\t-\n"
             << seq_filename3 + "\tseq1\t163\t247\t4\t-\n"
             << seq_filename3 + "\tseq2\t186\t327\t4\t-\n"
             << seq_filename3 + "\tseq3\t163\t284\t4\t-\n"
             << seq_filename3 + "\tseq1\t247\t400\t5\t-\n"
             << seq_filename3 + "\tseq2\t327\t480\t5\t-\n"
             << seq_filename3 + "\tseq3\t284\t481\t5\t-\n";
    }

    auto const && [data, batches] = read_traversal_file(traversal_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 2);
    EXPECT_EQ(data.merged_max_bin_map.size(), 1);
    EXPECT_EQ(data.merged_max_bin_map.at(2), 1);

    // check region_map
    std::vector<region> expected_file1_seq1_regions
    {
        /*SPLIT_BIN_0*/    {0, -1, 0, 209}, {1, -1, 209, 400},
    };
    std::vector<region> expected_file1_seq2_regions
    {
        /*SPLIT_BIN_0*/    {0, -1, 0, 289}, {1, -1, 289, 480},
    };
    std::vector<region> expected_file1_seq3_regions
    {
        /*SPLIT_BIN_0*/    {0, -1, 0, 209}, {1, -1, 209, 481},
    };
    std::vector<region> expected_file2_seq10_regions
    {
        /*MERGED_BIN_2_0*/ {2, 0, 0, 209}, {2, 1, 209, 400},
        /*MERGED_BIN_2_1*/ {2, 2, 0, 400},
    };
    std::vector<region> expected_file2_seq20_regions
    {
        /*MERGED_BIN_2_0*/ {2, 0, 0, 289}, {2, 1, 289, 480},
        /*MERGED_BIN_2_1*/ {2, 2, 0, 480},
    };
    std::vector<region> expected_file2_seq30_regions
    {
        /*MERGED_BIN_2_0*/ {2, 0, 0, 209}, {2, 1, 209, 481},
        /*MERGED_BIN_2_1*/ {2, 2, 0, 481},
    };
    std::vector<region> expected_file3_seq1_regions
    {
        /*SPLIT_BIN_3*/    {3, -1, 0, 163}, {4, -1, 163, 247}, {5, -1, 247, 400}
    };
    std::vector<region> expected_file3_seq2_regions
    {
        /*SPLIT_BIN_3*/    {3, -1, 0, 186}, {4, -1, 186, 327}, {5, -1, 327, 480}
    };
    std::vector<region> expected_file3_seq3_regions
    {
        /*SPLIT_BIN_3*/    {3, -1, 0, 163}, {4, -1, 163, 284}, {5, -1, 284, 481}
    };

    EXPECT_RANGE_EQ(data.region_map.at(seq_filename + "seq1"), expected_file1_seq1_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename + "seq2"), expected_file1_seq2_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename + "seq3"), expected_file1_seq3_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename2 + "seq10"), expected_file2_seq10_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename2 + "seq20"), expected_file2_seq20_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename2 + "seq30"), expected_file2_seq30_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename3 + "seq1"), expected_file3_seq1_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename3 + "seq2"), expected_file3_seq2_regions);
    EXPECT_RANGE_EQ(data.region_map.at(seq_filename3 + "seq3"), expected_file3_seq3_regions);

    // check if HIBF bins were calculated c
    EXPECT_EQ(data.hibf_num_technical_bins, 6);

    // check batches
    std::vector<batch> expected_batches
    {
        {{seq_filename2}, {2}, 3},
        {{seq_filename}, {0, 1}, 0},
        {{seq_filename3}, {3, 4, 5}, 0}
    };

    EXPECT_RANGE_EQ(batches, expected_batches);

    // check max high level IBF record
    EXPECT_EQ(*data.hibf_max_batch_record, expected_batches[0]); // MERGED_BIN_2
}

TEST(read_traversal_file_test, merged_file_in_the_end)
{
    std::string seq_filename = "test.fa";
    std::string seq_filename2 = "test2.fa";
    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    { // write example file
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_2 max_bin_id:1\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*SPLIT_BIN_0*/
             << seq_filename + "\tseq1\t0\t209\t0\t-\n"
             << seq_filename + "\tseq2\t0\t289\t0\t-\n"
             << seq_filename + "\tseq3\t0\t209\t0\t-\n"
             << seq_filename + "\tseq1\t209\t400\t1\t-\n"
             << seq_filename + "\tseq2\t289\t480\t1\t-\n"
             << seq_filename + "\tseq3\t209\t481\t1\t-\n"
             /*MERGED_BIN_2_0*/
             << seq_filename2 + "\tseq10\t0\t209\t2\t0\n"
             << seq_filename2 + "\tseq20\t0\t289\t2\t0\n"
             << seq_filename2 + "\tseq30\t0\t209\t2\t0\n"
             << seq_filename2 + "\tseq10\t209\t400\t2\t1\n"
             << seq_filename2 + "\tseq20\t289\t480\t2\t1\n"
             << seq_filename2 + "\tseq30\t209\t481\t2\t1\n"
             /*MERGED_BIN_2_1*/
             << seq_filename2 + "\tseq10\t0\t400\t2\t2\n"
             << seq_filename2 + "\tseq20\t0\t480\t2\t2\n"
             << seq_filename2 + "\tseq30\t0\t481\t2\t2\n";
    }

    auto const && [data, batches] = read_traversal_file(traversal_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 2);
    EXPECT_EQ(data.merged_max_bin_map.size(), 1);
    EXPECT_EQ(data.merged_max_bin_map.at(2), 1);

    // check if HIBF bins were calculated c
    EXPECT_EQ(data.hibf_num_technical_bins, 3);

    // check batches
    std::vector<batch> expected_batches
    {
        {{seq_filename2}, {2}, 3},
        {{seq_filename}, {0, 1}, 0}
    };

    EXPECT_RANGE_EQ(batches, expected_batches);
}

TEST(read_traversal_file_test, merged_file_with_distinct_files)
{
    std::string seq_filename = "test.fa";
    std::string seq_filename2 = "test2.fa";
    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    { // write example file
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_0 max_bin_id:1\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_0\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*MERGED_BIN_0_0*/
             << seq_filename2 + "\tseq10\t0\t209\t0\t0\n"
             << seq_filename2 + "\tseq20\t0\t289\t0\t0\n"
             << seq_filename2 + "\tseq30\t0\t209\t0\t0\n"
             << seq_filename2 + "\tseq10\t209\t400\t0\t1\n"
             << seq_filename2 + "\tseq20\t289\t480\t0\t1\n"
             << seq_filename2 + "\tseq30\t209\t481\t0\t1\n"
             /*MERGED_BIN_0_1*/
             << seq_filename + "\tseq1\t0\t400\t0\t2\n"
             << seq_filename + "\tseq2\t0\t480\t0\t2\n"
             << seq_filename + "\tseq3\t0\t481\t0\t2\n";
    }

    auto const && [data, batches] = read_traversal_file(traversal_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 0);
    EXPECT_EQ(data.merged_max_bin_map.size(), 1);
    EXPECT_EQ(data.merged_max_bin_map.at(0), 1);

    // check if HIBF bins were calculated c
    EXPECT_EQ(data.hibf_num_technical_bins, 1);

    // check batches
    std::vector<batch> expected_batches
    {
        {{seq_filename, seq_filename2}, {0}, 3}
    };

    EXPECT_RANGE_EQ(batches, expected_batches);
}

TEST(read_traversal_file_test, split_file_with_distinct_files)
{
    std::string seq_filename = "test.fa";
    std::string seq_filename2 = "test2.fa";
    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    { // write example file
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_3 max_bin_id:1\n"
             << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_3\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*SPLIT_BIN_0*/
             << seq_filename + "\tseq1\t0\t209\t0\t-\n"
             << seq_filename + "\tseq2\t0\t289\t0\t-\n"
             << seq_filename + "\tseq3\t0\t209\t0\t-\n"
             << seq_filename2 + "\tseq10\t0\t400\t2\t-\n"
             << seq_filename2 + "\tseq20\t0\t480\t2\t-\n"
             << seq_filename2 + "\tseq30\t0\t481\t2\t-\n"
             << seq_filename + "\tseq1\t209\t400\t1\t-\n"
             << seq_filename + "\tseq2\t289\t480\t1\t-\n"
             << seq_filename + "\tseq3\t209\t481\t1\t-\n"
             /*MERGED_BIN_3_0*/
             << seq_filename + "\tseq1\t0\t209\t3\t0\n"
             << seq_filename + "\tseq2\t0\t289\t3\t0\n"
             << seq_filename + "\tseq3\t0\t209\t3\t0\n"
             << seq_filename + "\tseq1\t209\t400\t3\t1\n"
             << seq_filename + "\tseq2\t289\t480\t3\t1\n"
             << seq_filename + "\tseq3\t209\t481\t3\t1\n";
    }

    auto const && [data, batches] = read_traversal_file(traversal_filename.get_path().string());

    // check data from header information
    EXPECT_EQ(data.hibf_max_bin, 3);
    EXPECT_EQ(data.merged_max_bin_map.size(), 1);
    EXPECT_EQ(data.merged_max_bin_map.at(3), 1);

    // check if HIBF bins were calculated c
    EXPECT_EQ(data.hibf_num_technical_bins, 4);

    // check batches
    std::vector<batch> expected_batches
    {
        {{seq_filename}, {3}, 2},
        {{seq_filename}, {0, 1}, 0},
        {{seq_filename2}, {2}, 0}
    };

    EXPECT_RANGE_EQ(batches, expected_batches);
}

TEST(read_traversal_file_test, split_bin_is_highest_bin)
{
    std::string seq_filename = "test.fa";
    seqan3::test::tmp_filename traversal_filename{"test.traverse"};

    { // write example file
        std::ofstream fout{traversal_filename.get_path()};
        fout << "#MERGED_BIN_3 max_bin_id:1\n"
             << "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_0\n"
             << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
             /*SPLIT_BIN_0*/
             << seq_filename + "\tseq1\t0\t209\t0\t-\n"
             << seq_filename + "\tseq1\t209\t400\t1\t-\n"
             << seq_filename + "\tseq10\t0\t400\t2\t-\n"
             /*MERGED_BIN_3_0*/
             << seq_filename + "\tseq1\t0\t209\t3\t0\n"
             << seq_filename + "\tseq2\t0\t289\t3\t0\n"
             << seq_filename + "\tseq3\t0\t209\t3\t0\n"
             << seq_filename + "\tseq1\t209\t400\t3\t1\n"
             << seq_filename + "\tseq2\t289\t480\t3\t1\n"
             << seq_filename + "\tseq3\t209\t481\t3\t1\n";
    }

    auto const && [data, batches] = read_traversal_file(traversal_filename.get_path().string());

    EXPECT_EQ(data.hibf_max_bin, 0);

    std::vector<batch> expected_batches
    {
        {{seq_filename}, {3}, 2},
        {{seq_filename}, {0, 1, 2}, 0}
    };

    EXPECT_RANGE_EQ(batches, expected_batches);

    ASSERT_NE(data.hibf_max_batch_record, nullptr);
    EXPECT_EQ(data.hibf_max_batch_record, &batches[1]);
    EXPECT_EQ(*data.hibf_max_batch_record, (batch{{seq_filename}, {0, 1, 2}, 0}));
}
