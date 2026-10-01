#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <chopper/layout/fast_layout_cluster.hpp>
#include <chopper/layout/phibf/multi_cluster.hpp>

using chopper::layout::Cluster;
using chopper::layout::phibf::MultiCluster;

TEST(Multicluster_test, ctor_from_cluster)
{
    Cluster const cluster1{5};
    MultiCluster const multi_cluster1{cluster1};

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
    Cluster cluster1{5};
    Cluster cluster2{7};
    cluster2.move_to(cluster1);

    MultiCluster const multi_cluster2{cluster2};

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
    Cluster cluster1{5};
    Cluster cluster2{7};
    cluster2.move_to(cluster1);
    ASSERT_EQ(cluster1.size(), 2u);
    Cluster const cluster3{13};

    MultiCluster multi_cluster1{cluster1};
    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
    EXPECT_EQ(multi_cluster1.size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins().size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins()[0].size(), 2u);

    MultiCluster multi_cluster3{cluster3};
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
    EXPECT_EQ(multi_cluster3.contained_user_bins()[0][0], cluster3.id());
    EXPECT_EQ(multi_cluster3.contained_user_bins()[1][0], cluster1.id());
    EXPECT_EQ(multi_cluster3.contained_user_bins()[1][1], cluster2.id());
}

TEST(Multicluster_test, sort_by_cardinality)
{
    std::vector<size_t> const cardinalities{10, 50, 20, 40, 30};

    Cluster cluster0{0};
    Cluster cluster1{1};
    Cluster cluster2{2};
    Cluster cluster3{3};
    cluster1.move_to(cluster0); // {0, 1}
    cluster3.move_to(cluster2); // {2, 3}
    Cluster cluster4{4};

    MultiCluster multi_cluster0{cluster0}; // {{0, 1}}
    MultiCluster multi_cluster4{cluster4}; // {{4}}
    MultiCluster multi_cluster2{cluster2}; // {{2, 3}}
    multi_cluster0.move_to(multi_cluster4); // {{4}, {0, 1}}
    multi_cluster2.move_to(multi_cluster4); // {{4}, {0, 1}, {2, 3}}

    multi_cluster4.sort_by_cardinality(cardinalities);

    auto const & clusters = multi_cluster4.contained_user_bins();
    ASSERT_EQ(clusters.size(), 3u);
    // inner clusters are sorted by descending cardinality, clusters by descending size
    EXPECT_EQ(clusters[2], (std::vector<size_t>{4}));
    for (size_t i = 0; i < 2; ++i)
    {
        ASSERT_EQ(clusters[i].size(), 2u);
        EXPECT_GE(cardinalities[clusters[i][0]], cardinalities[clusters[i][1]]);
    }
}

TEST(LSH_find_representative_cluster_test, multi_cluster_one_move)
{
    std::vector<MultiCluster> mclusters{{Cluster{0}}, {Cluster{1}}};
    mclusters[1].move_to(mclusters[0]);

    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(mclusters, mclusters[1].id()), mclusters[0].id());
}

TEST(LSH_find_representative_cluster_test, multi_cluster_two_moves)
{
    std::vector<MultiCluster> mclusters{{Cluster{0}}, {Cluster{1}}, {Cluster{2}}};
    mclusters[2].move_to(mclusters[1]);
    mclusters[1].move_to(mclusters[0]);

    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(mclusters, mclusters[1].id()), mclusters[0].id());
    EXPECT_EQ(chopper::layout::LSH_find_representative_cluster(mclusters, mclusters[2].id()), mclusters[0].id());
}
