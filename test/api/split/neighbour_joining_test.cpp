#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <fstream>

#include <seqan3/test/expect_range_eq.hpp>

#include <chopper/split/neighbour_joining.hpp>

TEST(neighbour_joining_test, small_example)
{
    // Create a distance matrix
    seqan::String<double> mat;
    seqan::resize(mat, 8*8, 0);
    seqan::assignValue(mat, 0*8+1, 7); seqan::assignValue(mat, 0*8+2, 8);seqan::assignValue(mat, 0*8+3, 11);seqan::assignValue(mat, 0*8+4, 13);seqan::assignValue(mat, 0*8+5, 16);seqan::assignValue(mat, 0*8+6, 13);seqan::assignValue(mat, 0*8+7, 17);
    seqan::assignValue(mat, 1*8+2, 5); seqan::assignValue(mat, 1*8+3, 8);seqan::assignValue(mat, 1*8+4, 10);seqan::assignValue(mat, 1*8+5, 13);seqan::assignValue(mat, 1*8+6, 10);seqan::assignValue(mat, 1*8+7, 14);
    seqan::assignValue(mat, 2*8+3, 5); seqan::assignValue(mat, 2*8+4, 7);seqan::assignValue(mat, 2*8+5, 10);seqan::assignValue(mat, 2*8+6, 7);seqan::assignValue(mat, 2*8+7, 11);
    seqan::assignValue(mat, 3*8+4, 8); seqan::assignValue(mat, 3*8+5, 11);seqan::assignValue(mat, 3*8+6, 8);seqan::assignValue(mat, 3*8+7, 12);
    seqan::assignValue(mat, 4*8+5, 5); seqan::assignValue(mat, 4*8+6, 6);seqan::assignValue(mat, 4*8+7, 10);
    seqan::assignValue(mat, 5*8+6, 9); seqan::assignValue(mat, 5*8+7, 13);
    seqan::assignValue(mat, 6*8+7, 8);

    auto guideTreeOut = neighbour_joining(mat);
    //std::cout << guideTreeOut << std::endl;

    EXPECT_EQ(seqan::numVertices(guideTreeOut), 15);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 8, 1) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 8, 0) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 9, 5) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 9, 4) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 2) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 8) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 2) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 8) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 11, 3) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 11, 10) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 12, 9) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 12, 11) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 13, 12) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 13, 6) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 14, 13) != 0);
    EXPECT_TRUE(seqan::findEdge(guideTreeOut, 14, 7) != 0);
    EXPECT_TRUE(seqan::getRoot(guideTreeOut) == 14);
}

TEST(neighbour_joining_test, wikipedia_example)
{ // see https://de.wikipedia.org/wiki/Neighbor-Joining-Algorithmus

    // Create a distance matrix
    seqan::String<double> mat;
    seqan::resize(mat, 4*4, 0);
    seqan::assignValue(mat, 0*4+1, 3); seqan::assignValue(mat, 0*4+2, 14);seqan::assignValue(mat, 0*4+3, 12);
    seqan::assignValue(mat, 1*4+2, 13); seqan::assignValue(mat, 1*4+3, 11);
    seqan::assignValue(mat, 2*4+3, 4);

    auto guideTreeOut = neighbour_joining(mat);
    std::cout << guideTreeOut << std::endl;

    std::ofstream dotFile("/tmp/nj.dot");
    writeRecords(dotFile, guideTreeOut, seqan::DotDrawing());

    EXPECT_EQ(seqan::numVertices(guideTreeOut), 7);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 8, 1) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 8, 0) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 9, 5) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 9, 4) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 2) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 8) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 2) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 10, 8) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 11, 3) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 11, 10) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 12, 9) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 12, 11) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 13, 12) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 13, 6) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 14, 13) != 0);
    //EXPECT_TRUE(seqan::findEdge(guideTreeOut, 14, 7) != 0);
    //EXPECT_TRUE(seqan::getRoot(guideTreeOut) == 14);
}
