#include <gtest/gtest.h>

#include <fstream>

#include <seqan/graph_algorithms.h>
#include <seqan/graph_types.h>

#include <chopper/split/map_distance_matrix.hpp>

TEST(map_distance_matrix_test, construction)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    EXPECT_EQ(m.nseq, 8u);
    EXPECT_EQ(m.dummy, 1.0);
    EXPECT_EQ(m.distance_threshold, 0.9);
}

TEST(map_distance_matrix_test, length)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    EXPECT_EQ(length(m), 64u);
}

TEST(map_distance_matrix_test, empty_matrix_has_default_values)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    for (size_t i = 0; i < m.nseq; ++i)
        for (size_t j = 0; j < m.nseq; ++j)
            EXPECT_EQ(m[i*m.nseq+j], 1.0);
}

TEST(map_distance_matrix_test, set_distance_value_distance_score)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    m.set_distance_value(0, 1, distance_score{0.5});
    m.set_distance_value(3, 4, distance_score{0.7});
    m.set_distance_value(3, 5, distance_score{0.95});

    EXPECT_EQ(m[0*m.nseq+1], 0.5);
    EXPECT_EQ(m[3*m.nseq+4], 0.7);
    EXPECT_EQ(m[3*m.nseq+5], 1.0);
}

TEST(map_distance_matrix_test, set_distance_value_similarity_score)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    m.set_distance_value(0, 1, similarity_score{0.5});
    m.set_distance_value(3, 4, similarity_score{0.3});
    m.set_distance_value(3, 5, similarity_score{0.05});

    EXPECT_EQ(m[0*m.nseq+1], 0.5);
    EXPECT_EQ(m[3*m.nseq+4], 0.7);
    EXPECT_EQ(m[3*m.nseq+5], 1.0);
}

TEST(map_distance_matrix_test, proxy_assignment)
{
    map_distance_matrix m{num_seq{8}, dummy_value{1.0}, upper_distance_threshold{0.9}};

    m[0*m.nseq+1] = 0.5;
    m[3*m.nseq+4] = 0.7;
    m[3*m.nseq+5] = 0.95;

    EXPECT_EQ(m[0*m.nseq+1], 0.5);
    EXPECT_EQ(m[3*m.nseq+4], 0.7);
    EXPECT_EQ(m[3*m.nseq+5], 1.0);

    m[0*m.nseq+1] = 0.95;
    m[3*m.nseq+4] = 0.95;

    EXPECT_EQ(m[0*m.nseq+1], 1.0);
    EXPECT_EQ(m[3*m.nseq+4], 1.0);
}
