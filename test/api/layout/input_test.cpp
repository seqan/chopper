#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <sstream> // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>  // for allocator, string
#include <vector>  // for vector

#include <chopper/layout/input.hpp>

std::string const layout_header{"@CHOPPER_USER_BINS\n"
                                "@0 seq0.fa\n"
                                "@1 seq1.fa\n"
                                "@2 seq2.fa\n"
                                "@3 seq3.fa\n"
                                "@4 seq4.fa\n"
                                "@5 seq5.fa\n"
                                "@6 seq6.fa\n"
                                "@7 seq7.fa\n"
                                "@CHOPPER_USER_BINS_END\n"
                                "@CHOPPER_CONFIG\n"
                                "@{\n"
                                "@    \"chopper_config\": {\n"
                                "@        \"version\": 2,\n"
                                "@        \"data_file\": {\n"
                                "@            \"value0\": \"input.txt\"\n"
                                "@        },\n"
                                "@        \"debug\": false,\n"
                                "@        \"sketch_directory\": {\n"
                                "@            \"value0\": \"\"\n"
                                "@        },\n"
                                "@        \"k\": 15,\n"
                                "@        \"window_size\": 15,\n"
                                "@        \"disable_sketch_output\": true,\n"
                                "@        \"precomputed_files\": false,\n"
                                "@        \"maximum_index_size\": 0,\n"
                                "@        \"number_of_partitions\": 0,\n"
                                "@        \"output_filename\": {\n"
                                "@            \"value0\": \"foo.layout\"\n"
                                "@        },\n"
                                "@        \"determine_best_tmax\": false,\n"
                                "@        \"force_all_binnings\": false\n"
                                "@    }\n"
                                "@}\n"
                                "@CHOPPER_CONFIG_END\n"
                                "@HIBF_CONFIG\n"
                                "@{\n"
                                "@    \"hibf_config\": {\n"
                                "@        \"version\": 1,\n"
                                "@        \"number_of_user_bins\": 3,\n"
                                "@        \"number_of_hash_functions\": 2,\n"
                                "@        \"maximum_fpr\": 0.05,\n"
                                "@        \"relaxed_fpr\": 0.3,\n"
                                "@        \"threads\": 2,\n"
                                "@        \"sketch_bits\": 12,\n"
                                "@        \"tmax\": 64,\n"
                                "@        \"alpha\": 1.2,\n"
                                "@        \"max_rearrangement_ratio\": 0.5,\n"
                                "@        \"disable_estimate_union\": false,\n"
                                "@        \"disable_rearrangement\": false\n"
                                "@    }\n"
                                "@}\n"
                                "@HIBF_CONFIG_END\n"};

TEST(layout_test, read_single_layout)
{
    // layout consists of three partitions, written one after the other
    std::stringstream ss{layout_header + R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file"};

    auto [filenames, chopper_config, layouts] = chopper::layout::read_layouts_file(ss);

    auto const & layout = layouts[0];

    EXPECT_EQ(layout.top_level_max_bin_id, 111);
    EXPECT_EQ(layout.max_bins[0], (seqan::hibf::layout::layout::max_bin{{0}, 0}));
    EXPECT_EQ(layout.max_bins[1], (seqan::hibf::layout::layout::max_bin{{2}, 2}));
    EXPECT_EQ(layout.max_bins[2], (seqan::hibf::layout::layout::max_bin{{1, 2, 3, 4}, 22}));
    EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 1, 7}));
    EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 22, 4}));
    EXPECT_EQ(layout.user_bins[2], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5}));
}

TEST(layout_test, read_from_partitioned_layout)
{
    // layout consists of three partitions, written one after the other
    std::stringstream ss{layout_header + R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file"};

    auto [filenames, chopper_config, hibf_layouts] = chopper::layout::read_layouts_file(ss);

    for (size_t i = 0; i < 3; ++i)
    {
        auto layout = hibf_layouts[i];

        EXPECT_EQ(layout.top_level_max_bin_id, 111);
        EXPECT_EQ(layout.max_bins[0], (seqan::hibf::layout::layout::max_bin{{0}, 0}));
        EXPECT_EQ(layout.max_bins[1], (seqan::hibf::layout::layout::max_bin{{2}, 2}));
        EXPECT_EQ(layout.max_bins[2], (seqan::hibf::layout::layout::max_bin{{1, 2, 3, 4}, 22}));
        EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 1, 7}));
        EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 22, 4}));
        EXPECT_EQ(layout.user_bins[2],
                  (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5}));
    }
}
