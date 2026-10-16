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
    config.number_of_partitions = 4;
    config.partitioning_approach = chopper::partitioning_scheme::lsh_sim;
    config.output_filename = "file.layout";
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.output_verbose_statistics = true;
    config.fast_layout = true;

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
           lhs.number_of_partitions == rhs.number_of_partitions &&                                 //
           lhs.partitioning_approach == rhs.partitioning_approach &&                               //
           lhs.output_filename == rhs.output_filename &&                                           //
           lhs.determine_best_tmax == rhs.determine_best_tmax &&                                   //
           lhs.force_all_binnings == rhs.force_all_binnings &&                                     //
           lhs.fast_layout == rhs.fast_layout &&                                                   //
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
                                                     "@        \"version\": 4,\n"
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
                                                     "@        \"output_filename\": {\n"
                                                     "@            \"value0\": \"file.layout\"\n"
                                                     "@        },\n"
                                                     "@        \"determine_best_tmax\": true,\n"
                                                     "@        \"force_all_binnings\": true,\n"
                                                     "@        \"fast_layout\": true,\n"
                                                     "@        \"number_of_partitions\": 4,\n"
                                                     "@        \"partitioning_approach\": 6\n"
                                                     "@    }\n"
                                                     "@}\n"
                                                     "@CHOPPER_CONFIG_END\n"
                                                     "@HIBF_CONFIG\n"
                                                     "@{\n"
                                                     "@    \"hibf_config\": {\n"
                                                     "@        \"version\": 3,\n"
                                                     "@        \"number_of_user_bins\": 123456789,\n"
                                                     "@        \"number_of_hash_functions\": 4,\n"
                                                     "@        \"maximum_fpr\": 0.0001,\n"
                                                     "@        \"relaxed_fpr\": 0.2,\n"
                                                     "@        \"threads\": 31,\n"
                                                     "@        \"sketch_bits\": 8,\n"
                                                     "@        \"tmax\": 128,\n"
                                                     "@        \"empty_bin_fraction\": 0.0,\n"
                                                     "@        \"track_occupancy\": false,\n"
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

// Layout files written before version 3 do not contain fast_layout. They must stay readable.
TEST(config_test, read_from_version_2)
{
    std::string config_string{config_string_view};

    auto replace = [&config_string](std::string_view const from, std::string_view const to)
    {
        size_t const pos = config_string.find(from);
        ASSERT_NE(pos, std::string::npos) << from;
        config_string.replace(pos, from.size(), to);
    };

    replace("\"version\": 4,\n@        \"data_file\"", "\"version\": 2,\n@        \"data_file\"");
    replace("\"force_all_binnings\": true,\n"
            "@        \"fast_layout\": true,\n"
            "@        \"number_of_partitions\": 4,\n"
            "@        \"partitioning_approach\": 6\n",
            "\"force_all_binnings\": true\n");

    std::stringstream ss{config_string};

    chopper::configuration config;
    config.read_from(ss);

    chopper::configuration expected{generate_config()};
    expected.fast_layout = false;
    expected.number_of_partitions = 1;
    expected.partitioning_approach = chopper::partitioning_scheme::blocked;

    EXPECT_EQ(config, expected);
}

// Configs written before version 4 do not contain the partitioned HIBF fields.
TEST(config_test, read_from_version_3)
{
    std::string config_string{config_string_view};

    auto erase_line = [&config_string](std::string_view const line)
    {
        size_t const pos = config_string.find(line);
        ASSERT_NE(pos, std::string::npos) << line;
        config_string.erase(pos, line.size());
    };
    erase_line(",\n@        \"number_of_partitions\": 4");
    erase_line(",\n@        \"partitioning_approach\": 6");

    std::string_view const version_4{"@        \"version\": 4,\n"};
    size_t const pos = config_string.find(version_4);
    ASSERT_NE(pos, std::string::npos);
    config_string.replace(pos, version_4.size(), "@        \"version\": 3,\n");

    std::stringstream ss{config_string};

    chopper::configuration config;
    config.read_from(ss);

    // The fields are not read and keep their default values.
    chopper::configuration expected{generate_config()};
    expected.number_of_partitions = 1;
    expected.partitioning_approach = chopper::partitioning_scheme::blocked;

    EXPECT_EQ(config, expected);
}
