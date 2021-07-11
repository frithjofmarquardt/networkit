/*
 * LouvainMapEquationBenchmark.cpp
 *
 * Created on: 2019-10-31
 * Author: Armin Wiebigke
 */
// networkit-format

#include <gtest/gtest.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/DirectedLouvainMapEquation.hpp>
#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

class DirectedMapEquationBenchmark : public testing::Test {
public:
    void SetUp() { Aux::Random::setSeed(435913, false); }
};

TEST_F(DirectedMapEquationBenchmark, benchLarge) {
    int numNodes = 1000;
    int numClusters = 20;
    int maxIter = 32;

    ClusteredRandomGraphGenerator generator(numNodes, numClusters, 0.5, 0.002);
    Graph G = generator.generate();
    Partition groundTruth = generator.getCommunities();


    Aux::Timer timer{};
    timer.start();

    LouvainMapEquation mapequation(G, false, maxIter, "none");
    mapequation.run();
    auto partition = mapequation.getPartition();

    timer.stop();

    Aux::Log::setLogLevel("INFO");
    INFO(mapequation.toString(), " took ", timer.elapsedMilliseconds(), "ms");
    INFO("Number Subsets = ", partition.numberOfSubsets(), " (Groundtruth = ", groundTruth.numberOfSubsets(), ")");
    

    Graph Gdirected(numNodes, true, true);

    G.forEdges([&](node u, node v, edgeweight weight, index edgeid) {
        Gdirected.addEdge(u,v,weight);
        Gdirected.addEdge(v,u,weight);  
    });

    timer.start();

    DirectedLouvainMapEquation dmapequation(Gdirected, maxIter);
    dmapequation.run();
    partition = dmapequation.getPartition();

    timer.stop();

    INFO(dmapequation.toString(), " took ", timer.elapsedMilliseconds(), "ms");
    INFO("Number Subsets = ", partition.numberOfSubsets(), " (Groundtruth = ", groundTruth.numberOfSubsets(), ")");
    
}

} // namespace NetworKit