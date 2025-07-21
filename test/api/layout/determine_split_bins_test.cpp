#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <numeric>  // for allocator, string
#include <vector>   // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_split_bins.hpp>

TEST(simple_split, first)
{
    chopper::configuration config{};

    std::vector<size_t> cardinalities{};
    for (size_t i{0}; i < 20; ++i)
        cardinalities.push_back(2400);
    cardinalities.push_back(200);
    for (size_t i{0}; i < 84; ++i)
        cardinalities.push_back(50);

    std::vector<size_t> positions(cardinalities.size());
    std::iota(positions.begin(), positions.end(), 0u);

    size_t num_technical_bins{63};
    size_t num_user_bins{21};

    std::vector<std::vector<size_t>> partitions(64);

    auto const [num_splits, max_size] =
        chopper::layout::determine_split_bins(config, positions, cardinalities, num_technical_bins, num_user_bins, partitions);

    EXPECT_EQ(num_splits, 61);
    EXPECT_EQ(max_size, 1452);

    size_t ub_idx{20};
    size_t tb_idx{partitions.size() - 1};

    // TBs:  0, 1, 2, 3, ...
    // UBs: 20,19,19,19,18,18,18,17,17,17,16,16,16,15,15,15,14,14,14,13,13,13,12,12,12,11,11,11,10,10,10,9,9,9,8,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,2,2,2,1,1,1,0,0,0,

    EXPECT_EQ(partitions[tb_idx][0], ub_idx); // single bin with card = 200
    --tb_idx;
    --ub_idx;

    for (; tb_idx > partitions.size() - num_splits - 1; --tb_idx)
    {
        for (size_t three = 0; three < 3; ++three)
        {
            EXPECT_EQ(partitions[tb_idx][0], ub_idx);
            --tb_idx;
        }
        ++tb_idx; // one too much in loop before
        --ub_idx;
    }
}