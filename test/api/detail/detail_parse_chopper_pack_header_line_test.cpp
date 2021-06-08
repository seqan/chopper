#include <gtest/gtest.h>

#include <seqan3/std/filesystem>
#include <sstream>
#include <fstream>

#include <lemon/graph_to_eps.h>

#include <chopper/detail_parse_chopper_pack_header_line.hpp>

#include "../api_test.hpp"

template <typename record_type>
void write_graph(lemon::ListDigraph const & g,
                 lemon::ListDigraph::NodeMap<node_data<record_type>> const & node_map,
                 std::filesystem::path const & graph_file_name)
{
    std::ofstream sout{graph_file_name};
    seqan3::debug_stream_type fout{sout};
    // std::cerr << "[LOG] Writing graph to " << graph_file_name << std::endl;

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

TEST(parse_bin_indices_test, one_number)
{
    std::string str{"123"};
    auto v = parse_bin_indices(str);

    EXPECT_EQ(v.front(), 123);
}

TEST(parse_bin_indices_test, with_semicolons)
{
    std::string str{"123;3;54"};
    auto v = parse_bin_indices(str);

    EXPECT_RANGE_EQ(v, (std::vector<size_t>{123, 3, 54}));
}

TEST(parse_chopper_pack_header_test, foo)
{
    std::string header
    {
        "#HIGH_LEVEL_IBF max_bin_id:0\n"
        "#MERGED_BIN_0 max_bin_id:1\n"
        "#MERGED_BIN_0;0;0 max_bin_id:30\n"
        "#MERGED_BIN_0;0 max_bin_id:4\n"
        "#MERGED_BIN_1 max_bin_id:26\n"
        "#MERGED_BIN_0;1 max_bin_id:34\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
    };

    std::stringstream header_stream{header};

    lemon::ListDigraph g;
    lemon::ListDigraph::NodeMap<node_data<chopper_pack_record>> node_map{g};

    parse_chopper_pack_header(g, node_map, header_stream);

    seqan3::test::tmp_filename dot_file{"header_graph.dot"};
    write_graph(g, node_map, dot_file.get_path().c_str());
}
