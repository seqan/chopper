#include <gtest/gtest.h>

#include <fstream>

#include <chopper/detail_hibf_user_bins.hpp>

#include "../api_test.hpp"


TEST(hibf_user_bins_test, access_vector)
{
    hibf_user_bins user_bins{};

    user_bins.resize_bins(2);
    user_bins.resize_filename(2);

    user_bins.filename_at(0) = "foo;bar";
    user_bins.filename_at(1) = "bar";

    user_bins.bin_at(0) = std::vector<int64_t>{0, 0, 1};
    user_bins.bin_at(1) = std::vector<int64_t>{-1, 0, 1};

    EXPECT_RANGE_EQ(user_bins[0], (std::vector<std::string>{"foo;bar", "foo;bar", "bar"}));
    EXPECT_RANGE_EQ(user_bins[1], (std::vector<std::string>{"", "foo;bar", "bar"}));
}

TEST(hibf_user_bins_test, cerealize)
{
    hibf_user_bins user_bins{};

    user_bins.resize_bins(2);
    user_bins.resize_filename(2);

    user_bins.filename_at(0) = "foo;bar";
    user_bins.filename_at(1) = "bar";

    user_bins.bin_at(0) = std::vector<int64_t>{0, 0, 1};
    user_bins.bin_at(1) = std::vector<int64_t>{-1, 0, 1};

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
