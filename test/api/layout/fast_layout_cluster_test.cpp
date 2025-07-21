#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <chopper/layout/fast_layout_cluster.hpp>

TEST(Cluster_test, ctor_from_id)
{
    size_t user_bin_idx{5};
    chopper::layout::Cluster const cluster{user_bin_idx};

    EXPECT_EQ(cluster.id(), user_bin_idx);
    EXPECT_FALSE(cluster.empty());
    EXPECT_EQ(cluster.size(), 1u);
    ASSERT_EQ(cluster.contained_user_bins().size(), 1u);
    EXPECT_EQ(cluster.contained_user_bins()[0], user_bin_idx);
    EXPECT_TRUE(cluster.is_valid(user_bin_idx));
}

TEST(Cluster_test, move_to)
{
    size_t user_bin_idx1{5};
    size_t user_bin_idx2{7};
    chopper::layout::Cluster cluster1{user_bin_idx1};
    chopper::layout::Cluster cluster2{user_bin_idx2};

    EXPECT_TRUE(cluster1.is_valid(user_bin_idx1));
    EXPECT_TRUE(cluster2.is_valid(user_bin_idx2));

    cluster2.move_to(cluster1);

    // cluster1 now contains user bins 5 and 7
    EXPECT_EQ(cluster1.size(), 2u);
    ASSERT_EQ(cluster1.contained_user_bins().size(), 2u);
    EXPECT_EQ(cluster1.contained_user_bins()[0], user_bin_idx1);
    EXPECT_EQ(cluster1.contained_user_bins()[1], user_bin_idx2);

    // cluster 2 is empty
    EXPECT_TRUE(cluster2.has_been_moved());
    EXPECT_TRUE(cluster2.empty());
    EXPECT_EQ(cluster2.size(), 0u);
    EXPECT_EQ(cluster2.contained_user_bins().size(), 0u);
    EXPECT_EQ(cluster2.moved_to_cluster_id(), cluster1.id());

    // both should still be valid
    EXPECT_TRUE(cluster1.is_valid(user_bin_idx1));
    EXPECT_TRUE(cluster2.is_valid(user_bin_idx2));
}

TEST(LSH_find_representative_cluster_test, cluster_one_move)
{
    std::vector<chopper::layout::Cluster> clusters{chopper::layout::Cluster{0}, chopper::layout::Cluster{1}};
    clusters[1].move_to(clusters[0]);

    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(clusters, clusters[1].id()), clusters[0].id());
}

TEST(LSH_find_representative_cluster_test, cluster_two_moves)
{
    std::vector<chopper::layout::Cluster> clusters{chopper::layout::Cluster{0},
                                                   chopper::layout::Cluster{1},
                                                   chopper::layout::Cluster{2}};
    clusters[2].move_to(clusters[1]);
    clusters[1].move_to(clusters[0]);

    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(clusters, clusters[1].id()), clusters[0].id());
    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(clusters, clusters[2].id()), clusters[0].id());
}
