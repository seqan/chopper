#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/build/read_chopper_pack_file.hpp>

TEST(read_chopper_pack_file_test, small_example)
{
    seqan3::test::tmp_filename chopper_split_filename{"test.split"};

    { // write example file
        std::ofstream fout{chopper_split_filename.get_path()};
        fout << "#HIGH_LEVEL_IBF max_bin_id:1\n"
             << "#MERGED_BIN_0;0;0 max_bin_id:30\n"
             << "#MERGED_BIN_0;0 max_bin_id:4\n"
             << "#MERGED_BIN_0;1 max_bin_id:34\n"
             << "#MERGED_BIN_0 max_bin_id:1\n"
             << "#MERGED_BIN_1 max_bin_id:26\n"
             << "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
             << "user_bin_15\t0;0;0;0\t1;1;1;30\t1650;350;80;1\n"
             << "user_bin_16\t0;0;0;30\t1;1;1;11\t1650;350;80;2\n"
             << "user_bin_17\t0;0;0;41\t1;1;1;11\t1650;350;80;2\n"
             << "user_bin_18\t0;0;0;52\t1;1;1;6\t1650;350;80;2\n"
             << "user_bin_19\t0;0;0;58\t1;1;1;6\t1650;350;80;2\n"
             << "user_bin_14\t0;0;1\t1;1;1\t1650;350;40\n"
             << "user_bin_13\t0;0;2\t1;1;1\t1650;350;50\n"
             << "user_bin_12\t0;0;3\t1;1;1\t1650;350;80\n"
             << "user_bin_11\t0;0;4\t1;1;1\t1650;350;100\n"
             << "user_bin_8\t0;1;0\t1;1;34\t1650;400;6\n"
             << "user_bin_9\t0;1;34\t1;1;15\t1650;400;7\n"
             << "user_bin_10\t0;1;49\t1;1;15\t1650;400;7\n"
             << "user_bin_7\t0;2\t1;1\t1650;200\n"
             << "user_bin_6\t0;3\t1;1\t1650;300\n"
             << "user_bin_5\t0;4\t1;1\t1650;400\n"
             << "user_bin_1\t2\t1\t1200\n"
             << "user_bin_0\t3\t2\t1500\n"
             << "user_bin_2\t1;0\t1;26\t2000;31\n"
             << "user_bin_3\t1;26\t1;19\t2000;32\n"
             << "user_bin_4\t1;45\t1;19\t2000;32\n";
    }

    build_data<chopper_pack_record> data{};

    read_chopper_pack_file(data, chopper_split_filename.get_path().string());

    /*graph: naming is [node_id: ibf name based on header information]
     *
     *  [0: HIGH_LEVEL_IBF]--> [1: MERGED_BIN_0]--> [3: MERGED_BIN_0;0]--> [5: MERGED_BIN_0;0;0]
     *                   |                    |
     *                   |                     ---> [4: MERGED_BIN_0;1]
     *                   |
     *                    ---> [2: MERGED_BIN_1]
     */

    EXPECT_EQ(lemon::countNodes(data.ibf_graph), 6);
    EXPECT_EQ(lemon::countArcs(data.ibf_graph), 5);

    /* node_data:
     * size_t parent_bin_index{};
     * size_t max_bin_index{};
     * size_t number_of_technical_bins{};
     * lemon::ListDigraph::Node favourite_child{lemon::INVALID};
     * std::vector<chopper_pack_record> remaining_records{};
     */

    /* chopper pack record:
     * std::vector<std::string> filenames{};
     * std::vector<size_t> bin_indices{};
     * std::vector<size_t> number_of_bins{};
     * std::vector<size_t> estimated_sizes{};
     */

    std::vector<node_data<chopper_pack_record>> expected_node_data
    {
        /*0*/ {0, 1, 5, data.ibf_graph.nodeFromId(2), {{{"user_bin_1"}, {2}, {1}, {1200}},
                                                       {{"user_bin_0"}, {3}, {2}, {1500}}}},
        /*1*/ {0, 1, 5, data.ibf_graph.nodeFromId(4), {{{"user_bin_7"}, {0,2}, {1,1}, {1650,200}},
                                                       {{"user_bin_6"}, {0,3}, {1,1}, {1650,300}},
                                                       {{"user_bin_5"}, {0,4}, {1,1}, {1650,400}}}},
        /*2*/ {1, 26, 64, lemon::INVALID, {{{"user_bin_3"}, {1,26}, {1,19}, {2000,32}},
                                           {{"user_bin_2"}, {1,0}, {1,26}, {2000,31}},
                                           {{"user_bin_4"}, {1,45}, {1,19}, {2000,32}}}},
        /*3*/ {0, 4, 5, lemon::INVALID, {{{"user_bin_11"}, {0,0,4}, {1,1,1}, {1650,350,100}},
                                         {{"user_bin_14"}, {0,0,1}, {1,1,1}, {1650,350,40}},
                                         {{"user_bin_13"}, {0,0,2}, {1,1,1}, {1650,350,50}},
                                         {{"user_bin_12"}, {0,0,3}, {1,1,1}, {1650,350,80}}}},
        /*4*/ {1, 34, 64, lemon::INVALID, {{{"user_bin_9"}, {0,1,34}, {1,1,15}, {1650,400,7}},
                                           {{"user_bin_8"}, {0,1,0}, {1,1,34}, {1650,400,6}},
                                           {{"user_bin_10"}, {0,1,49}, {1,1,15}, {1650,400,7}}}},
        /*5*/ {0, 30, 64, lemon::INVALID, {{{"user_bin_16"}, {0,0,0,30}, {1,1,1,11}, {1650,350,80,2}},
                                           {{"user_bin_15"}, {0,0,0,0}, {1,1,1,30}, {1650,350,80,1}},
                                           {{"user_bin_17"}, {0,0,0,41}, {1,1,1,11}, {1650,350,80,2}},
                                           {{"user_bin_18"}, {0,0,0,52}, {1,1,1,6}, {1650,350,80,2}},
                                           {{"user_bin_19"}, {0,0,0,58}, {1,1,1,6}, {1650,350,80,2}}}}
    };

    ASSERT_EQ(expected_node_data.size(), lemon::countNodes(data.ibf_graph)); // sanity check for the next loop
    for (size_t i = 0; i < expected_node_data.size(); ++i)
    {
        lemon::ListDigraph::Node nd = data.ibf_graph.nodeFromId(i);
        auto & nd_data = data.node_map[nd];
        EXPECT_EQ(nd_data, expected_node_data[i]);
    }
}
