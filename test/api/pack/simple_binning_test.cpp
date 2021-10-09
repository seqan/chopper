#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <chopper/pack/hibf_model.hpp>
#include <chopper/pack/pack_config.hpp>
#include <chopper/pack/simple_binning.hpp>

TEST(simple_binning_test, small_example)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &output_buffer;
    data.kmer_counts = {100, 40, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};
    data.fp_correction = std::vector<double>(65, 1.0);

    hibf_model hibf(pack_config{}, data.fp_correction);
    simple_binning algo{data,hibf.get_top_level_ibf(), 9};
    size_t max_bin = algo.execute();

    std::string expected
    {
        "seq4\t;0\t;1\n"
        "seq3\t;1\t;1\n"
        "seq2\t;2\t;2\n"
        "seq1\t;4\t;5\n"
    };

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, uniform_distribution)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &output_buffer;
    data.kmer_counts = {20, 20, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};
    data.fp_correction = std::vector<double>(65, 1.0);

    hibf_model hibf(pack_config{}, data.fp_correction);
    simple_binning algo{data, hibf.get_top_level_ibf(),  4u};
    size_t max_bin = algo.execute();

    std::string expected
    {
        "seq4\t;0\t;1\n"
        "seq3\t;1\t;1\n"
        "seq2\t;2\t;1\n"
        "seq1\t;3\t;1\n"
    };

    EXPECT_EQ(output_buffer.str(), expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    std::stringstream output_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &output_buffer;
    data.kmer_counts = {100, 40, 20, 20};
    data.filenames = {"seq1", "seq2", "seq3", "seq4"};
    data.fp_correction = std::vector<double>(65, 1.0);

    hibf_model hibf(pack_config{}, data.fp_correction);
    EXPECT_THROW((simple_binning{data, hibf.get_top_level_ibf(), 2}), std::runtime_error);
}
