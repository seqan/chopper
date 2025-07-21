#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <chopper/lsh.hpp>

TEST(Cluster_test, ctor_from_id)
{
    size_t user_bin_idx{5};
    chopper::Cluster const cluster{user_bin_idx};

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
    chopper::Cluster cluster1{user_bin_idx1};
    chopper::Cluster cluster2{user_bin_idx2};

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

TEST(Multicluster_test, ctor_from_cluster)
{
    chopper::Cluster const cluster1{5};
    chopper::MultiCluster const multi_cluster1{cluster1};

    EXPECT_EQ(multi_cluster1.id(), cluster1.id());
    EXPECT_FALSE(multi_cluster1.empty());
    EXPECT_EQ(multi_cluster1.size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins().size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins()[0].size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins()[0][0], 5u);
    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
}

TEST(Multicluster_test, ctor_from_moved_cluster)
{
    chopper::Cluster cluster1{5};
    chopper::Cluster cluster2{7};
    cluster2.move_to(cluster1);

    chopper::MultiCluster const multi_cluster2{cluster2};

    EXPECT_EQ(multi_cluster2.id(), cluster2.id());
    EXPECT_TRUE(multi_cluster2.empty());
    EXPECT_TRUE(multi_cluster2.has_been_moved());
    EXPECT_EQ(multi_cluster2.moved_to_cluster_id(), cluster2.moved_to_cluster_id());
    EXPECT_EQ(multi_cluster2.size(), 0u);
    ASSERT_EQ(multi_cluster2.contained_user_bins().size(), 0u);
    EXPECT_TRUE(multi_cluster2.is_valid(cluster2.id()));
}

TEST(Multicluster_test, move_to)
{
    chopper::Cluster cluster1{5};
    chopper::Cluster cluster2{7};
    cluster2.move_to(cluster1);
    ASSERT_EQ(cluster1.size(), 2u);
    chopper::Cluster const cluster3{13};

    chopper::MultiCluster multi_cluster1{cluster1};
    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
    EXPECT_EQ(multi_cluster1.size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins().size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins()[0].size(), 2u);

    chopper::MultiCluster multi_cluster3{cluster3};
    EXPECT_TRUE(multi_cluster3.is_valid(cluster3.id()));
    EXPECT_EQ(multi_cluster3.size(), 1u);

    multi_cluster1.move_to(multi_cluster3);

    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
    EXPECT_TRUE(multi_cluster3.is_valid(cluster3.id()));

    // multi_cluster1 has been moved and is empty now
    EXPECT_TRUE(multi_cluster1.empty());
    EXPECT_EQ(multi_cluster1.size(), 0u);
    EXPECT_TRUE(multi_cluster1.has_been_moved());
    EXPECT_EQ(multi_cluster1.moved_to_cluster_id(), multi_cluster3.id());

    // multi_cluster3 contains 2 clusters now, {13} and {5, 7}
    EXPECT_FALSE(multi_cluster3.empty());
    EXPECT_EQ(multi_cluster3.size(), 2u); // two clusters
    ASSERT_EQ(multi_cluster3.contained_user_bins().size(), 2u);
    ASSERT_EQ(multi_cluster3.contained_user_bins()[0].size(), 1u);
    ASSERT_EQ(multi_cluster3.contained_user_bins()[1].size(), 2u);
    ASSERT_EQ(multi_cluster3.contained_user_bins()[0][0], cluster3.id());
    ASSERT_EQ(multi_cluster3.contained_user_bins()[1][0], cluster1.id());
    ASSERT_EQ(multi_cluster3.contained_user_bins()[1][1], cluster2.id());
}

TEST(LSH_find_representative_cluster_test, cluster_one_move)
{
    std::vector<chopper::Cluster> clusters{chopper::Cluster{0}, chopper::Cluster{1}};
    clusters[1].move_to(clusters[0]);

    EXPECT_EQ(chopper::LSH_find_representative_cluster(clusters, clusters[1].id()), clusters[0].id());
}

TEST(LSH_find_representative_cluster_test, multi_cluster_one_move)
{
    std::vector<chopper::MultiCluster> mclusters{{chopper::Cluster{0}}, {chopper::Cluster{1}}};
    mclusters[1].move_to(mclusters[0]);

    EXPECT_EQ(chopper::LSH_find_representative_cluster(mclusters, mclusters[1].id()), mclusters[0].id());

}

TEST(LSH_find_representative_cluster_test, cluster_two_moves)
{
    std::vector<chopper::Cluster> clusters{chopper::Cluster{0}, chopper::Cluster{1}, chopper::Cluster{2}};
    clusters[2].move_to(clusters[1]);
    clusters[1].move_to(clusters[0]);

    EXPECT_EQ(chopper::LSH_find_representative_cluster(clusters, clusters[1].id()), clusters[0].id());
    EXPECT_EQ(chopper::LSH_find_representative_cluster(clusters, clusters[2].id()), clusters[0].id());
}

TEST(LSH_find_representative_cluster_test, multi_cluster_two_moves)
{
    std::vector<chopper::MultiCluster> mclusters{{chopper::Cluster{0}}, {chopper::Cluster{1}}, {chopper::Cluster{2}}};
    mclusters[2].move_to(mclusters[1]);
    mclusters[1].move_to(mclusters[0]);

    EXPECT_EQ(chopper::LSH_find_representative_cluster(mclusters, mclusters[1].id()), mclusters[0].id());
    EXPECT_EQ(chopper::LSH_find_representative_cluster(mclusters, mclusters[2].id()), mclusters[0].id());
}
