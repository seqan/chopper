#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/detail_hibf_user_bins.hpp>

TEST(hibf_user_bins_test, small_test)
{
    hibf_user_bins user_bins{};

    user_bins.add_user_bin("foo");
    user_bins.add_user_bin("bar");

    user_bins.add_user_bin_positions({0, 0, 1});

    EXPECT_EQ((user_bins[{0, 0}]), std::string{"foo"});
    EXPECT_EQ((user_bins[{0, 1}]), std::string{"foo"});
    EXPECT_EQ((user_bins[{0, 2}]), std::string{"bar"});
}

TEST(hibf_user_bins_test, with_concatenation)
{
    hibf_user_bins user_bins{};

    user_bins.add_user_bin(std::vector<std::string>{"foo", "bar"});
    user_bins.add_user_bin("bar");

    user_bins.add_user_bin_positions({0, 0, 1});

    EXPECT_EQ((user_bins[{0, 0}]), std::string{"foo;bar"});
    EXPECT_EQ((user_bins[{0, 1}]), std::string{"foo;bar"});
    EXPECT_EQ((user_bins[{0, 2}]), std::string{"bar"});
}

TEST(hibf_user_bins_test, access_vector)
{
    hibf_user_bins user_bins{};

    user_bins.add_user_bin(std::vector<std::string>{"foo", "bar"});
    user_bins.add_user_bin("bar");

    user_bins.add_user_bin_positions({0, 0, 1});
    user_bins.add_user_bin_positions({-1, 0, 1});

    EXPECT_RANGE_EQ(user_bins[0], (std::vector<std::string>{"foo;bar", "foo;bar", "bar"}));
    EXPECT_RANGE_EQ(user_bins[1], (std::vector<std::string>{"", "foo;bar", "bar"}));
}

TEST(hibf_user_bins_test, cerealize)
{
    hibf_user_bins user_bins{};

    user_bins.add_user_bin(std::vector<std::string>{"foo", "bar"});
    user_bins.add_user_bin("bar");

    user_bins.add_user_bin_positions({0, 0, 1});
    user_bins.add_user_bin_positions({-1, 0, 1});

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
