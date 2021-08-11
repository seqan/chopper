#include <gtest/gtest.h>

#include <fstream>

#include <chopper/hierarchical_interleaved_bloom_filter.hpp>

#include "../api_test.hpp"

using hibf_user_bins = typename hibf::hierarchical_interleaved_bloom_filter<>::user_bins;

TEST(hibf_user_bins_test, access_vector)
{
    hibf_user_bins user_bins{};

    user_bins.set_ibf_count(2);
    user_bins.set_user_bin_count(2);

    user_bins.filename_of_user_bin(0) = "foo;bar";
    user_bins.filename_of_user_bin(1) = "bar";

    user_bins.bin_indices_of_ibf(0) = std::vector<int64_t>{0, 0, 1};
    user_bins.bin_indices_of_ibf(1) = std::vector<int64_t>{-1, 0, 1};

    EXPECT_RANGE_EQ(user_bins[0], (std::vector<std::string>{"foo;bar", "foo;bar", "bar"}));
    EXPECT_RANGE_EQ(user_bins[1], (std::vector<std::string>{"", "foo;bar", "bar"}));
}

TEST(hibf_user_bins_test, cerealize)
{
    hibf_user_bins user_bins{};

    user_bins.set_ibf_count(2);
    user_bins.set_user_bin_count(2);

    user_bins.filename_of_user_bin(0) = "foo;bar";
    user_bins.filename_of_user_bin(1) = "bar";

    user_bins.bin_indices_of_ibf(0) = std::vector<int64_t>{0, 0, 1};
    user_bins.bin_indices_of_ibf(1) = std::vector<int64_t>{-1, 0, 1};

    seqan3::test::tmp_filename filename{"user_bins.out"};

    { // write
        std::ofstream fout(filename.get_path(), std::ios::binary);

        cereal::BinaryOutputArchive archive(fout);
        archive(user_bins);
    }

    hibf_user_bins read_user_bins{};

    { // read
        std::ifstream is{filename.get_path(), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(read_user_bins);
    }

    EXPECT_RANGE_EQ(read_user_bins[0], (std::vector<std::string>{"foo;bar", "foo;bar", "bar"}));
    EXPECT_RANGE_EQ(read_user_bins[1], (std::vector<std::string>{"", "foo;bar", "bar"}));
}
