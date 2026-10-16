#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <sstream> // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>  // for allocator, string
#include <vector>  // for vector

#include <chopper/layout/input.hpp>

// A layout file consists of one section per (partitioned HIBF) layout. Each section has its own user bin filenames,
// chopper config, and layout. User bin indices within a section are local to that section.
std::string layout_section(std::vector<std::string> const & filenames,
                           size_t const number_of_partitions,
                           std::string const & layout)
{
    std::string section{"@CHOPPER_USER_BINS\n"};
    for (size_t i = 0; i < filenames.size(); ++i)
        section += "@" + std::to_string(i) + " " + filenames[i] + "\n";

    section += "@CHOPPER_USER_BINS_END\n"
               "@CHOPPER_CONFIG\n"
               "@{\n"
               "@    \"chopper_config\": {\n"
               "@        \"version\": 4,\n"
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
               "@        \"output_filename\": {\n"
               "@            \"value0\": \"foo.layout\"\n"
               "@        },\n"
               "@        \"determine_best_tmax\": false,\n"
               "@        \"force_all_binnings\": false,\n"
               "@        \"fast_layout\": false,\n"
               "@        \"number_of_partitions\": "
             + std::to_string(number_of_partitions)
             + ",\n"
               "@        \"partitioning_approach\": 1\n"
               "@    }\n"
               "@}\n"
               "@CHOPPER_CONFIG_END\n"
               "@HIBF_CONFIG\n"
               "@{\n"
               "@    \"hibf_config\": {\n"
               "@        \"version\": 3,\n"
               "@        \"number_of_user_bins\": "
             + std::to_string(filenames.size())
             + ",\n"
               "@        \"number_of_hash_functions\": 2,\n"
               "@        \"maximum_fpr\": 0.05,\n"
               "@        \"relaxed_fpr\": 0.3,\n"
               "@        \"threads\": 2,\n"
               "@        \"sketch_bits\": 12,\n"
               "@        \"tmax\": 64,\n"
               "@        \"empty_bin_fraction\": 0.0,\n"
               "@        \"track_occupancy\": false,\n"
               "@        \"alpha\": 1.2,\n"
               "@        \"max_rearrangement_ratio\": 0.5,\n"
               "@        \"disable_estimate_union\": false,\n"
               "@        \"disable_rearrangement\": false\n"
               "@    }\n"
               "@}\n"
               "@HIBF_CONFIG_END\n"
             + layout;
    return section;
}

std::string const layout_with_three_user_bins{R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
2	0	1
0	1;0	1;22
1	1;2;3;4;22	1;1;1;1;21
)layout_file"};

std::string const layout_with_two_user_bins{R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:0
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
1	0	40
0	40	24
)layout_file"};

void check_layout_with_three_user_bins(seqan::hibf::layout::layout const & layout)
{
    EXPECT_EQ(layout.top_level_max_bin_id, 111);
    ASSERT_EQ(layout.max_bins.size(), 3u);
    EXPECT_EQ(layout.max_bins[0], (seqan::hibf::layout::layout::max_bin{{0}, 0}));
    EXPECT_EQ(layout.max_bins[1], (seqan::hibf::layout::layout::max_bin{{2}, 2}));
    EXPECT_EQ(layout.max_bins[2], (seqan::hibf::layout::layout::max_bin{{1, 2, 3, 4}, 22}));
    ASSERT_EQ(layout.user_bins.size(), 3u);
    EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 1, 2}));
    EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 22, 0}));
    EXPECT_EQ(layout.user_bins[2], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1, 2, 3, 4}, 22, 21, 1}));
}

void check_layout_with_two_user_bins(seqan::hibf::layout::layout const & layout)
{
    EXPECT_EQ(layout.top_level_max_bin_id, 0);
    EXPECT_TRUE(layout.max_bins.empty());
    ASSERT_EQ(layout.user_bins.size(), 2u);
    EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 40, 1}));
    EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 40, 24, 0}));
}

TEST(layout_test, read_single_layout)
{
    std::stringstream ss{layout_section({"seq0.fa", "seq1.fa", "seq2.fa"}, 1u, layout_with_three_user_bins)};

    auto [filenames, chopper_config, layout] = chopper::layout::read_layout_file(ss);

    EXPECT_EQ(filenames, (std::vector<std::vector<std::string>>{{"seq0.fa"}, {"seq1.fa"}, {"seq2.fa"}}));
    EXPECT_EQ(chopper_config.k, 15u);
    EXPECT_EQ(chopper_config.number_of_partitions, 1u);
    EXPECT_EQ(chopper_config.hibf_config.number_of_user_bins, 3u);
    check_layout_with_three_user_bins(layout);

    // The whole file has been consumed.
    EXPECT_FALSE(ss.good());
}

TEST(layout_test, read_from_partitioned_layout)
{
    // The layout consists of three partitions, each written as its own section, one after the other.
    std::stringstream ss{layout_section({"seq0.fa", "seq1.fa", "seq2.fa"}, 3u, layout_with_three_user_bins)
                         + layout_section({"seq3.fa", "seq4.fa"}, 3u, layout_with_two_user_bins)
                         + layout_section({"seq5.fa", "seq6.fa", "seq7.fa"}, 3u, layout_with_three_user_bins)};

    std::vector<std::vector<std::vector<std::string>>> filenames{};
    std::vector<chopper::configuration> configs{};
    std::vector<seqan::hibf::layout::layout> layouts{};

    while (ss.good())
    {
        auto [section_filenames, section_config, section_layout] = chopper::layout::read_layout_file(ss);
        filenames.push_back(std::move(section_filenames));
        configs.push_back(std::move(section_config));
        layouts.push_back(std::move(section_layout));
    }

    ASSERT_EQ(layouts.size(), 3u);

    // Each section has its own filenames.
    EXPECT_EQ(filenames[0], (std::vector<std::vector<std::string>>{{"seq0.fa"}, {"seq1.fa"}, {"seq2.fa"}}));
    EXPECT_EQ(filenames[1], (std::vector<std::vector<std::string>>{{"seq3.fa"}, {"seq4.fa"}}));
    EXPECT_EQ(filenames[2], (std::vector<std::vector<std::string>>{{"seq5.fa"}, {"seq6.fa"}, {"seq7.fa"}}));

    // Each section has its own config, whose number of user bins is the one of its partition.
    for (auto const & config : configs)
    {
        EXPECT_EQ(config.k, 15u);
        EXPECT_EQ(config.number_of_partitions, 3u);
        EXPECT_EQ(config.partitioning_approach, chopper::partitioning_scheme::sorted);
    }
    EXPECT_EQ(configs[0].hibf_config.number_of_user_bins, 3u);
    EXPECT_EQ(configs[1].hibf_config.number_of_user_bins, 2u);
    EXPECT_EQ(configs[2].hibf_config.number_of_user_bins, 3u);

    // Each layout keeps the user bin indices local to its section.
    check_layout_with_three_user_bins(layouts[0]);
    check_layout_with_two_user_bins(layouts[1]);
    check_layout_with_three_user_bins(layouts[2]);
}
