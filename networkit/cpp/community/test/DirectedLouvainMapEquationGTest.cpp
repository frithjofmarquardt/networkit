/*
 * DirectedLouvainMapEquationGTest.cpp
 *
 * Created on: 2021-07-04
 * Author: Frithjof Marquardt
 */

// networkit-format

#include <gtest/gtest.h>

#include <networkit/community/DirectedLouvainMapEquation.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class DirectedMapEquationGTest : public testing::Test {
public:
    void SetUp() { Aux::Random::setSeed(435913, false); }
};

void addDirectedClique(Graph &graph, Partition &groundTruth, node lowestId, node highestId) {
    index subsetId = groundTruth.upperBound();
    groundTruth.setUpperBound(subsetId + 1);
    for (node i = lowestId; i <= highestId; ++i) {
        groundTruth.addToSubset(subsetId, i);
        for (node j = i + 1; j <= highestId; ++j) {
            graph.addEdge(i, j, 1);
            graph.addEdge(j, i, 1);
        }
    }
}

TEST_F(DirectedMapEquationGTest, testLocalMoveSmall) {
    Aux::Random::setSeed(2342556, false);
    Graph G(10, true, true);
    Partition groundTruth(10);
    addDirectedClique(G, groundTruth, 0, 4);
    addDirectedClique(G, groundTruth, 5, 9);
    G.addEdge(0, 9, 1);
    G.addEdge(9, 0, 1);
    G.addEdge(7, 2, 1);
    G.addEdge(2, 7, 1);
    DirectedLouvainMapEquation mapequation(G, 256);
    mapequation.run();
    auto partition = mapequation.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

} // namespace NetworKit