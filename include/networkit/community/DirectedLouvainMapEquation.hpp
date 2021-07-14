/*
 * DirectedLouvainMapEquation.hpp
 *
 * Created on: 2021-07-04
 * Author: Frithjof Marquardt
 *
 * Based on LouvainMapEquation.hpp by: 
 *         Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */
// networkit-format

#ifndef NETWORKIT_COMMUNITY_DIRECTED_LOUVAIN_MAP_EQUATION_HPP_
#define NETWORKIT_COMMUNITY_DIRECTED_LOUVAIN_MAP_EQUATION_HPP_

#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

class DirectedLouvainMapEquation : public CommunityDetectionAlgorithm {
public:
    /**
     * @param[in] G input graph
     * @param[in] maxIterations maximum number of iterations for move phase
     *
     */
    explicit DirectedLouvainMapEquation(const Graph &graph, count maxIterations = 32);

    void run() override;

    std::string toString() const override;

private:
    count maxIterations;

    std::vector<double> clusterCut, clusterVolume;
    std::vector<count> clusterSize;
    std::vector<double> nodeFrequencies;
    std::vector<double> weightedOutDegrees;
    std::vector<SparseVector<double>> ets_neighborClusterWeights;

    double totalCut, totalVolume;
    double tau;
    double fitness;

    count localMoving(std::vector<node> &nodes, count iteration);

    bool tryLocalMove(node u, SparseVector<double> &neighborClusterWeights);

    bool performMove(node u, double frequency, node currentCluster, node targetCluster, 
                     double currentCutDelta, double targetCutDelta, double fitnessDelta);

    void calculateNodeFrequencies();

    double computeAbsoluteDifference(std::vector<double> v, std::vector<double> w);

    void calculateInitialClusterCutAndVolume();

    /**
     * Calculate the change in the map equation if the node is moved from its current cluster to the
     * target cluster.
     */
    double fitnessChange(node u, double frequency, double currentCutDelta, double targetCutDelta, 
                         node currentCluster, node targetCluster);

    double plogp(double w);

    
#ifndef NDEBUG
    long double mapEquation();

    long double oldfitness;
    long double newfitness;
    double lastChange;
#endif

};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_DIRECTED_LOUVAIN_MAP_EQUATION_HPP_
