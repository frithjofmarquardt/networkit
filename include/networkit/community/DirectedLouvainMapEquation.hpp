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
    struct Move {
        node movedNode = none;
        double volume = 0.0;
        index originCluster = none, targetCluster = none;
        double cutUpdateToOriginCluster = 0.0, cutUpdateToTargetCluster = 0.0;

        Move(node n, double vol, index cc, index tc, double cuptoc, double cupttc)
            : movedNode(n), volume(vol), originCluster(cc), targetCluster(tc),
              cutUpdateToOriginCluster(cuptoc), cutUpdateToTargetCluster(cupttc) {}
    };

    static_assert(std::is_trivially_destructible<Move>::value,
                  "DirectedLouvainMapEquation::Move struct is not trivially destructible");

    count maxIterations;

    std::vector<double> clusterCut, clusterVolume;
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
