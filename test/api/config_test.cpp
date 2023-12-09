#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <chopper/configuration.hpp> // for config
#include <chopper/layout/input.hpp>

chopper::configuration generate_config()
{
    chopper::configuration config{};

    config.data_file = "/path/to/data.file";
    config.sketch_directory = "/path/to/sketch/directory";
    config.k = 20;
    config.window_size = 24;
    config.disable_sketch_output = true;
    config.precomputed_files = true;
    config.output_filename = "file.layout";
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.output_verbose_statistics = true;

    config.hibf_config.number_of_user_bins = 123456789;
    config.hibf_config.number_of_hash_functions = 4;
    config.hibf_config.maximum_fpr = 0.0001;
    config.hibf_config.relaxed_fpr = 0.2;
    config.hibf_config.threads = 31;
    config.hibf_config.sketch_bits = 8;
    config.hibf_config.tmax = 128;
    config.hibf_config.alpha = 1.0;
    config.hibf_config.max_rearrangement_ratio = 0.333;
    config.hibf_config.disable_estimate_union = true;
    config.hibf_config.disable_rearrangement = false;

    return config;
}

namespace chopper
{

bool operator==(chopper::configuration const & lhs, chopper::configuration const & rhs)
{
    return lhs.data_file == rhs.data_file &&                                                       //
           lhs.debug == rhs.debug &&                                                               //
           lhs.sketch_directory == rhs.sketch_directory &&                                         //
           lhs.k == rhs.k &&                                                                       //
           lhs.window_size == rhs.window_size &&                                                   //
           lhs.disable_sketch_output == rhs.disable_sketch_output &&                               //
           lhs.precomputed_files == rhs.precomputed_files &&                                       //
           lhs.output_filename == rhs.output_filename &&                                           //
           lhs.determine_best_tmax == rhs.determine_best_tmax &&                                   //
           lhs.force_all_binnings == rhs.force_all_binnings &&                                     //
           lhs.hibf_config.number_of_user_bins == rhs.hibf_config.number_of_user_bins &&           //
           lhs.hibf_config.number_of_hash_functions == rhs.hibf_config.number_of_hash_functions && //
           lhs.hibf_config.maximum_fpr == rhs.hibf_config.maximum_fpr &&                           //
           lhs.hibf_config.relaxed_fpr == rhs.hibf_config.relaxed_fpr &&                           //
           lhs.hibf_config.threads == rhs.hibf_config.threads &&                                   //
           lhs.hibf_config.sketch_bits == rhs.hibf_config.sketch_bits &&                           //
           lhs.hibf_config.tmax == rhs.hibf_config.tmax &&                                         //
           lhs.hibf_config.alpha == rhs.hibf_config.alpha &&                                       //
           lhs.hibf_config.max_rearrangement_ratio == rhs.hibf_config.max_rearrangement_ratio &&   //
           lhs.hibf_config.disable_estimate_union == rhs.hibf_config.disable_estimate_union &&     //
           lhs.hibf_config.disable_rearrangement == rhs.hibf_config.disable_rearrangement;
}

} // namespace chopper

static constexpr std::string_view config_string_view{"@CHOPPER_CONFIG\n"
                                                     "@{\n"
                                                     "@    \"chopper_config\": {\n"
                                                     "@        \"version\": 2,\n"
                                                     "@        \"data_file\": {\n"
                                                     "@            \"value0\": \"/path/to/data.file\"\n"
                                                     "@        },\n"
                                                     "@        \"debug\": false,\n"
                                                     "@        \"sketch_directory\": {\n"
                                                     "@            \"value0\": \"/path/to/sketch/directory\"\n"
                                                     "@        },\n"
                                                     "@        \"k\": 20,\n"
                                                     "@        \"window_size\": 24,\n"
                                                     "@        \"disable_sketch_output\": true,\n"
                                                     "@        \"precomputed_files\": true,\n"
                                                     "@        \"maximum_index_size\": 0,\n"
                                                     "@        \"number_of_partitions\": 0,\n"
                                                     "@        \"output_filename\": {\n"
                                                     "@            \"value0\": \"file.layout\"\n"
                                                     "@        },\n"
                                                     "@        \"determine_best_tmax\": true,\n"
                                                     "@        \"force_all_binnings\": true\n"
                                                     "@    }\n"
                                                     "@}\n"
                                                     "@CHOPPER_CONFIG_END\n"
                                                     "@HIBF_CONFIG\n"
                                                     "@{\n"
                                                     "@    \"hibf_config\": {\n"
                                                     "@        \"version\": 1,\n"
                                                     "@        \"number_of_user_bins\": 123456789,\n"
                                                     "@        \"number_of_hash_functions\": 4,\n"
                                                     "@        \"maximum_fpr\": 0.0001,\n"
                                                     "@        \"relaxed_fpr\": 0.2,\n"
                                                     "@        \"threads\": 31,\n"
                                                     "@        \"sketch_bits\": 8,\n"
                                                     "@        \"tmax\": 128,\n"
                                                     "@        \"alpha\": 1.0,\n"
                                                     "@        \"max_rearrangement_ratio\": 0.333,\n"
                                                     "@        \"disable_estimate_union\": true,\n"
                                                     "@        \"disable_rearrangement\": false\n"
                                                     "@    }\n"
                                                     "@}\n"
                                                     "@HIBF_CONFIG_END\n"};

TEST(config_test, write_to)
{
    chopper::configuration const config{generate_config()};

    std::stringstream ss{};
    config.write_to(ss);

    EXPECT_EQ(ss.str(), config_string_view);
}

TEST(config_test, read_from)
{
    std::string config_string{config_string_view};
    std::stringstream ss{config_string};

    chopper::configuration config;
    config.read_from(ss);

    EXPECT_EQ(config, generate_config());
}

TEST(config_test, read_from_with_more_meta)
{
    std::string config_string{"@blah some other stuff\n"
                              "@blah some other stuff\n"
                              "@blah some other stuff\n"
                              "@blah some other stuff\n"
                              "@blah some other stuff\n"};
    config_string += config_string_view;
    std::stringstream ss{config_string};

    chopper::configuration config;
    config.read_from(ss);

    EXPECT_EQ(config, generate_config());
}

// Easier to do in the config_test because of existing helper functions
TEST(input, read_layouts_file)
{
    std::string config_string{"@CHOPPER_USER_BINS\n"
                              "@0 file1.fa\n"
                              "@1 file2.fa\n"
                              "@2 path/to/file3.fa\n"
                              "@3 file4.fastq\n"
                              "@CHOPPER_USER_BINS_END\n"};
    config_string += config_string_view;

    {
        seqan::hibf::layout::layout hibf_layout{};
        hibf_layout.top_level_max_bin_id = 42u;
        std::stringstream config_out{};
        hibf_layout.write_to(config_out);
        config_string += config_out.str();
    }

    std::stringstream ss{config_string};

    auto [filenames, config, layouts] = chopper::layout::read_layouts_file(ss);

    auto const & layout = layouts[0];

    std::vector<std::vector<std::string>> const expected_filenames{{"file1.fa"},
                                                                   {"file2.fa"},
                                                                   {"path/to/file3.fa"},
                                                                   {"file4.fastq"}};

    EXPECT_EQ(filenames, expected_filenames);
    EXPECT_EQ(config, generate_config());
    EXPECT_EQ(layout.top_level_max_bin_id, 42u);
    EXPECT_TRUE(layout.max_bins.empty());
    EXPECT_TRUE(layout.user_bins.empty());
}
