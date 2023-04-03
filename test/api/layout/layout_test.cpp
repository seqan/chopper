#include <gtest/gtest.h>

#include <sstream>

#include <chopper/layout/layout.hpp>

TEST(layout_test, printing_max_bins)
{
    std::stringstream ss{};

    chopper::layout::layout layout;

    layout.max_bins.emplace_back(std::vector<size_t>{}, 0);
    layout.max_bins.emplace_back(std::vector<size_t>{2}, 2);
    layout.max_bins.emplace_back(std::vector<size_t>{1,2,3,4}, 22);

    for (auto const & mb : layout.max_bins)
        ss << mb << "\n";

    std::string expected = R"mb(#MERGED_BIN_ max_bin_id:0
#MERGED_BIN_2 max_bin_id:2
#MERGED_BIN_1;2;3;4 max_bin_id:22
)mb";

    EXPECT_EQ(ss.str(), expected);
}

TEST(layout_test, printing_user_bins)
{
    std::stringstream ss{};

    chopper::layout::layout layout;

    layout.user_bins.emplace_back("seq7", std::vector<size_t>{}, 1, 0);
    layout.user_bins.emplace_back("seq4", std::vector<size_t>{1}, 22, 0);
    layout.user_bins.emplace_back("seq5", std::vector<size_t>{1,2,3,4}, 21, 22);

    for (auto const & ub : layout.user_bins)
        ss << ub << "\n";

    std::string expected = R"ub(seq7	0	1
seq4	1;0	1;22
seq5	1;2;3;4;22	1;1;1;1;21
)ub";

    EXPECT_EQ(ss.str(), expected);
}
