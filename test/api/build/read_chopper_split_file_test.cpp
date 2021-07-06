#include <gtest/gtest.h>

#include <chopper/build/read_chopper_split_file.hpp>

#include "../api_test.hpp"

TEST(parse_chopper_split_line_test, single_bin)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10\n"};

    auto && [filename, id, bin_indices, reg] = parse_chopper_split_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_RANGE_EQ(bin_indices, (std::vector<size_t>{10}));
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.bin_index, 10);
}

TEST(parse_chopper_split_line_test, multiple_bins)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10;5;3\n"};

    auto && [filename, id, bin_indices, reg] = parse_chopper_split_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_RANGE_EQ(bin_indices, (std::vector<size_t>{10, 5, 3}));
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.bin_index, 3);
}

TEST(parse_chopper_split_line_test, without_new_line)
{
    std::string const example_line{"file1\tseq1\t216\t400\t10;5;3"};

    auto && [filename, id, bin_indices, reg] = parse_chopper_split_line(example_line);

    EXPECT_EQ(filename, "file1");
    EXPECT_EQ(id, "seq1");
    EXPECT_RANGE_EQ(bin_indices, (std::vector<size_t>{10, 5, 3}));
    EXPECT_EQ(reg.begin, 216);
    EXPECT_EQ(reg.end, 400);
    EXPECT_EQ(reg.bin_index, 3);
}

TEST(read_chopper_split_file_test, small_example)
{
    seqan3::test::tmp_filename chopper_split_filename{"test.split"};

    { // write example file
        std::ofstream fout{chopper_split_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:2\n"
             << "#MERGED_BIN_2 max_bin_id:1\n"
             << "#FILES\tSEQ_ID\tBEGIN\tEND\tBIN_INDICES\n"
             << "test1.fa\tseq1\t0\t209\t0\n"
             << "test1.fa\tseq2\t0\t289\t0\n"
             << "test1.fa\tseq3\t0\t209\t0\n"
             << "test1.fa\tseq1\t209\t400\t1\n"
             << "test1.fa\tseq2\t289\t480\t1\n"
             << "test1.fa\tseq3\t209\t481\t1\n"

             << "test2.fa\tseq10\t0\t209\t2;0\n"
             << "test2.fa\tseq20\t0\t289\t2;0\n"
             << "test2.fa\tseq30\t0\t209\t2;0\n"
             << "test2.fa\tseq10\t209\t400\t2;1\n"
             << "test2.fa\tseq20\t289\t480\t2;1\n"
             << "test2.fa\tseq30\t209\t481\t2;1\n"
             << "test2.fa\tseq10\t0\t400\t2;2\n"
             << "test2.fa\tseq20\t0\t480\t2;2\n"
             << "test2.fa\tseq30\t0\t481\t2;2\n"

             << "test3.fa\tseq1\t0\t163\t3\n"
             << "test3.fa\tseq2\t0\t186\t3\n"
             << "test3.fa\tseq3\t0\t163\t3\n"
             << "test3.fa\tseq1\t163\t247\t4\n"
             << "test3.fa\tseq2\t186\t327\t4\n"
             << "test3.fa\tseq3\t163\t284\t4\n"
             << "test3.fa\tseq1\t247\t400\t5\n"
             << "test3.fa\tseq2\t327\t480\t5\n"
             << "test3.fa\tseq3\t284\t481\t5\n";
    }

    build_data<chopper_split_record> data{};

    read_chopper_split_file(data, chopper_split_filename.get_path().string());

    /*graph: naming is [node_id: ibf name based on header information]
     *
     *  [0: HIGH_LEVEL_IBF]--> [1: MERGED_BIN_2]
     */

    EXPECT_EQ(lemon::countNodes(data.ibf_graph), 2);
    EXPECT_EQ(lemon::countArcs(data.ibf_graph), 1);

    /* node_data:
     * size_t parent_bin_index{};
     * size_t max_bin_index{};
     * size_t number_of_technical_bins{};
     * lemon::ListDigraph::Node favourite_child{lemon::INVALID};
     * std::vector<chopper_pack_record> remaining_records{};
     */

    /* chopper split record:
     * std::vector<std::string> filenames{};
     * std::vector<size_t> bin_indices{};
     * std::vector<size_t> number_of_bins{};
     * robin_hood::unordered_map<std::string, std::vector<region>> region_map{};
     */

    std::vector<node_data<chopper_split_record>> expected_node_data
    {
        /*0*/ {0, 2, 6, data.ibf_graph.nodeFromId(1),
                    {{{"test1.fa"}, {0,1}, {1}, {{"test1.faseq1", {{0,0,209},{1,209,400}}},
                                                 {"test1.faseq2", {{0,0,289},{1,289,480}}},
                                                 {"test1.faseq3", {{0,0,209},{1,209,481}}}}},
                     {{"test3.fa"}, {3,4,5}, {1}, {{"test3.faseq1", {{3,0,163},{4,163,247},{5,247,400}}},
                                                   {"test3.faseq2", {{3,0,186},{4,186,327},{5,327,480}}},
                                                   {"test3.faseq3", {{3,0,163},{4,163,284},{5,284,481}}}}}}},
        /*1*/ {2, 1, 3, lemon::INVALID,
                    {{{"test2.fa"}, {1,0,2}, {1}, {{"test2.faseq10", {{0,0,209},{1,209,400},{2,0,400}}},
                                                   {"test2.faseq20", {{0,0,289},{1,289,480},{2,0,480}}},
                                                   {"test2.faseq30", {{0,0,209},{1,209,481},{2,0,481}}}}}}}
    };


    ASSERT_EQ(expected_node_data.size(), lemon::countNodes(data.ibf_graph)); // sanity check for the next loop
    for (size_t i = 0; i < expected_node_data.size(); ++i)
    {
        lemon::ListDigraph::Node nd = data.ibf_graph.nodeFromId(i);
        auto & nd_data = data.node_map[nd];
        EXPECT_EQ(nd_data, expected_node_data[i]);
    }
}

// TEST(read_chopper_split_file_test, merged_file_in_the_end)
// {
//     std::string seq_filename = "test.fa";
//     std::string seq_filename2 = "test2.fa";
//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};

//     { // write example file
//         std::ofstream fout{chopper_split_filename.get_path()};
//         fout << "#MERGED_BIN_2 max_bin_id:1\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_2\n"
//              << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//              /*SPLIT_BIN_0*/
//              << seq_filename + "\tseq1\t0\t209\t0\t-\n"
//              << seq_filename + "\tseq2\t0\t289\t0\t-\n"
//              << seq_filename + "\tseq3\t0\t209\t0\t-\n"
//              << seq_filename + "\tseq1\t209\t400\t1\t-\n"
//              << seq_filename + "\tseq2\t289\t480\t1\t-\n"
//              << seq_filename + "\tseq3\t209\t481\t1\t-\n"
//              /*MERGED_BIN_2_0*/
//              << seq_filename2 + "\tseq10\t0\t209\t2\t0\n"
//              << seq_filename2 + "\tseq20\t0\t289\t2\t0\n"
//              << seq_filename2 + "\tseq30\t0\t209\t2\t0\n"
//              << seq_filename2 + "\tseq10\t209\t400\t2\t1\n"
//              << seq_filename2 + "\tseq20\t289\t480\t2\t1\n"
//              << seq_filename2 + "\tseq30\t209\t481\t2\t1\n"
//              /*MERGED_BIN_2_1*/
//              << seq_filename2 + "\tseq10\t0\t400\t2\t2\n"
//              << seq_filename2 + "\tseq20\t0\t480\t2\t2\n"
//              << seq_filename2 + "\tseq30\t0\t481\t2\t2\n";
//     }

//     auto const && [data, batches] = read_chopper_split_file(chopper_split_filename.get_path().string());

//     // check data from header information
//     EXPECT_EQ(data.hibf.hibf_max_bin, 2);
//     EXPECT_EQ(data.merged_max_bin_map.size(), 1);
//     EXPECT_EQ(data.merged_max_bin_map.at(2), 1);

//     // check if HIBF bins were calculated c
//     EXPECT_EQ(data.hibf.hibf_num_technical_bins, 3);

//     // check batches
//     std::vector<batch> expected_batches
//     {
//         {{seq_filename2}, {2}, 3},
//         {{seq_filename}, {0, 1}, 0}
//     };

//     EXPECT_RANGE_EQ(batches, expected_batches);
// }

// TEST(read_chopper_split_file_test, merged_file_with_distinct_files)
// {
//     std::string seq_filename = "test.fa";
//     std::string seq_filename2 = "test2.fa";
//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};

//     { // write example file
//         std::ofstream fout{chopper_split_filename.get_path()};
//         fout << "#MERGED_BIN_0 max_bin_id:1\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_0\n"
//              << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//              /*MERGED_BIN_0_0*/
//              << seq_filename2 + "\tseq10\t0\t209\t0\t0\n"
//              << seq_filename2 + "\tseq20\t0\t289\t0\t0\n"
//              << seq_filename2 + "\tseq30\t0\t209\t0\t0\n"
//              << seq_filename2 + "\tseq10\t209\t400\t0\t1\n"
//              << seq_filename2 + "\tseq20\t289\t480\t0\t1\n"
//              << seq_filename2 + "\tseq30\t209\t481\t0\t1\n"
//              /*MERGED_BIN_0_1*/
//              << seq_filename + "\tseq1\t0\t400\t0\t2\n"
//              << seq_filename + "\tseq2\t0\t480\t0\t2\n"
//              << seq_filename + "\tseq3\t0\t481\t0\t2\n";
//     }

//     auto const && [data, batches] = read_chopper_split_file(chopper_split_filename.get_path().string());

//     // check data from header information
//     EXPECT_EQ(data.hibf.hibf_max_bin, 0);
//     EXPECT_EQ(data.merged_max_bin_map.size(), 1);
//     EXPECT_EQ(data.merged_max_bin_map.at(0), 1);

//     // check if HIBF bins were calculated c
//     EXPECT_EQ(data.hibf.hibf_num_technical_bins, 1);

//     // check batches
//     std::vector<batch> expected_batches
//     {
//         {{seq_filename, seq_filename2}, {0}, 3}
//     };

//     EXPECT_RANGE_EQ(batches, expected_batches);
// }

// TEST(read_chopper_split_file_test, split_file_with_distinct_files)
// {
//     std::string seq_filename = "test.fa";
//     std::string seq_filename2 = "test2.fa";
//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};

//     { // write example file
//         std::ofstream fout{chopper_split_filename.get_path()};
//         fout << "#MERGED_BIN_3 max_bin_id:1\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:MERGED_BIN_3\n"
//              << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//              /*SPLIT_BIN_0*/
//              << seq_filename + "\tseq1\t0\t209\t0\t-\n"
//              << seq_filename + "\tseq2\t0\t289\t0\t-\n"
//              << seq_filename + "\tseq3\t0\t209\t0\t-\n"
//              << seq_filename2 + "\tseq10\t0\t400\t2\t-\n"
//              << seq_filename2 + "\tseq20\t0\t480\t2\t-\n"
//              << seq_filename2 + "\tseq30\t0\t481\t2\t-\n"
//              << seq_filename + "\tseq1\t209\t400\t1\t-\n"
//              << seq_filename + "\tseq2\t289\t480\t1\t-\n"
//              << seq_filename + "\tseq3\t209\t481\t1\t-\n"
//              /*MERGED_BIN_3_0*/
//              << seq_filename + "\tseq1\t0\t209\t3\t0\n"
//              << seq_filename + "\tseq2\t0\t289\t3\t0\n"
//              << seq_filename + "\tseq3\t0\t209\t3\t0\n"
//              << seq_filename + "\tseq1\t209\t400\t3\t1\n"
//              << seq_filename + "\tseq2\t289\t480\t3\t1\n"
//              << seq_filename + "\tseq3\t209\t481\t3\t1\n";
//     }

//     auto const && [data, batches] = read_chopper_split_file(chopper_split_filename.get_path().string());

//     // check data from header information
//     EXPECT_EQ(data.hibf.hibf_max_bin, 3);
//     EXPECT_EQ(data.merged_max_bin_map.size(), 1);
//     EXPECT_EQ(data.merged_max_bin_map.at(3), 1);

//     // check if HIBF bins were calculated c
//     EXPECT_EQ(data.hibf.hibf_num_technical_bins, 4);

//     // check batches
//     std::vector<batch> expected_batches
//     {
//         {{seq_filename}, {3}, 2},
//         {{seq_filename}, {0, 1}, 0},
//         {{seq_filename2}, {2}, 0}
//     };

//     EXPECT_RANGE_EQ(batches, expected_batches);
// }

// TEST(read_chopper_split_file_test, split_bin_is_highest_bin)
// {
//     std::string seq_filename = "test.fa";
//     seqan3::test::tmp_filename chopper_split_filename{"test.split"};

//     { // write example file
//         std::ofstream fout{chopper_split_filename.get_path()};
//         fout << "#MERGED_BIN_3 max_bin_id:1\n"
//              << "#HIGH_LEVEL_IBF max_bin_id:SPLIT_BIN_0\n"
//              << "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n"
//              /*SPLIT_BIN_0*/
//              << seq_filename + "\tseq1\t0\t209\t0\t-\n"
//              << seq_filename + "\tseq1\t209\t400\t1\t-\n"
//              << seq_filename + "\tseq10\t0\t400\t2\t-\n"
//              /*MERGED_BIN_3_0*/
//              << seq_filename + "\tseq1\t0\t209\t3\t0\n"
//              << seq_filename + "\tseq2\t0\t289\t3\t0\n"
//              << seq_filename + "\tseq3\t0\t209\t3\t0\n"
//              << seq_filename + "\tseq1\t209\t400\t3\t1\n"
//              << seq_filename + "\tseq2\t289\t480\t3\t1\n"
//              << seq_filename + "\tseq3\t209\t481\t3\t1\n";
//     }

//     auto const && [data, batches] = read_chopper_split_file(chopper_split_filename.get_path().string());

//     EXPECT_EQ(data.hibf.hibf_max_bin, 0);

//     std::vector<batch> expected_batches
//     {
//         {{seq_filename}, {3}, 2},
//         {{seq_filename}, {0, 1, 2}, 0}
//     };

//     EXPECT_RANGE_EQ(batches, expected_batches);

//     ASSERT_NE(data.hibf.hibf_max_batch_record, nullptr);
//     EXPECT_EQ(data.hibf.hibf_max_batch_record, &batches[1]);
//     EXPECT_EQ(*data.hibf.hibf_max_batch_record, (batch{{seq_filename}, {0, 1, 2}, 0}));
// }
