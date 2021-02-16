#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/build/read_chopper_pack_file.hpp>

void write_graph(lemon::ListDigraph const & g,
                 lemon::ListDigraph::NodeMap<node_data<chopper_pack_record>> const & node_map,
                 std::filesystem::path const & graph_file_name)
{
    std::ofstream sout{graph_file_name};
    seqan3::debug_stream_type fout{sout};
    std::cerr << "[LOG] Writing graph to " << graph_file_name << std::endl;

    std::string header = R"(graph G {
/* Graph Attributes */
graph [rankdir = LR];

/* Node Attributes */
node [shape = rectangle, fillcolor = white, style = filled, fontname = "Times-Italic"];

/* Edge Attributes */
edge [fontname = "Times-Italic", arrowsize = 0.75, fontsize = 16];

/* Nodes */
)";

    fout << header;

    for (lemon::ListDigraph::NodeIt u(g); u != lemon::INVALID; ++u)
        fout << g.id(u) << " [label = \"" << g.id(u) << ":"<< node_map[u] << "\", group = " << 0 << "];" << std::endl;

    fout << std::endl << "/*directed edges*/" << std::endl;

    for (lemon::ListDigraph::ArcIt a(g); a != lemon::INVALID; ++a)
        fout << g.id(g.source(a)) << " -- " << g.id(g.target(a)) << " [];" << std::endl;

    fout << std::endl << "}" << std::endl;
}

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
             << "user_bin_0\t3\t2\t1500;\n"
             << "user_bin_2\t1;0\t1;26\t2000;31\n"
             << "user_bin_3\t1;26\t1;19\t2000;32\n"
             << "user_bin_4\t1;45\t1;19\t2000;32\n";
    }

    build_data<chopper_pack_record> data{};
    read_chopper_pack_file(data, chopper_split_filename.get_path().string());

    write_graph(data.ibf_graph, data.node_map, "/tmp/ibf_graph.dot");

    // check data from header information

}
