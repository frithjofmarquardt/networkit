/*
 * DirectedLouvainMapEquation.cpp
 *
 * Created on: 2021-07-04
 * Author: Frithjof Marquardt
 * 
 * Based on by LouvainMapEquation.cpp:
 *         Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */
// networkit-format

#include <algorithm>
#include <cassert>
#include <cmath>

#include <omp.h>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/DirectedLouvainMapEquation.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

// cannot use stadard CommunityDetectionAlgorithm constructor as in LouvainMapEquation since it does not support directed graphs
DirectedLouvainMapEquation::DirectedLouvainMapEquation(const Graph &graph, count maxIterations)
    : CommunityDetectionAlgorithm(graph, Partition(graph.upperNodeIdBound())),
      maxIterations(maxIterations), clusterCut(graph.upperNodeIdBound()),
      clusterVolume(graph.upperNodeIdBound()), nodeFrequencies(graph.upperNodeIdBound(), 0),
      weightedOutDegrees(G->upperNodeIdBound(), 0), ets_neighborClusterWeights(1) {
        tau = 0.15;
        fitness = -1;;
      }

void DirectedLouvainMapEquation::run() {
#ifndef NDEBUG
    Aux::Log::setLogLevel("DEBUG");
#endif
    if (hasRun)
        throw std::runtime_error("Algorithm was already run!");

    if (!G->isDirected())
        throw std::runtime_error("Graph has to be directed to use this algorithm!");

    Aux::SignalHandler handler;

    result.allToSingletons();

    // remove unused nodes from clustering
    if (G->numberOfNodes() != G->upperNodeIdBound()) {
        for (node u = 0; u < G->upperNodeIdBound(); ++u) {
            if (!G->hasNode(u)) {
                result.remove(u);
            }
        }
    }
    handler.assureRunning();
    
    // calculate out degrees for relative edgeweights
    for (node u = 0; u < G->upperNodeIdBound(); u++) {
        if (G->hasNode(u))
        {
            weightedOutDegrees[u] = G->weightedDegree(u, false);
        }
    }
    // initial setup for algorithm
    calculateNodeFrequencies();
    calculateInitialClusterCutAndVolume();

    std::vector<node> nodes{G->nodeRange().begin(), G->nodeRange().end()};
    count numberOfNodesMoved = 1;
    
    for (count iteration = 0; iteration < maxIterations && numberOfNodesMoved > 0; ++iteration) {
        handler.assureRunning();
        
        numberOfNodesMoved = localMoving(nodes, iteration);
    }
    

    hasRun = true;
}

count DirectedLouvainMapEquation::localMoving(std::vector<node> &nodes, count iteration) {
    count nodesMoved = 0;

    SparseVector<double> &neighborClusterWeights = ets_neighborClusterWeights[0];
    if (iteration == 0) {
        neighborClusterWeights.resize(G->upperNodeIdBound(), -1.0);
    }
    for (node u : nodes) {
        if (tryLocalMove(u, neighborClusterWeights)) {
            nodesMoved += 1;
        }
    }

    return nodesMoved;
}

// try to move node u to clusters of neighboring nodes (all in and outgoing neighbors)
bool DirectedLouvainMapEquation::tryLocalMove(node u, SparseVector<double> &neighborClusterWeights) {
    
    const index currentCluster = result[u];
    const double frequency = nodeFrequencies[u];
    const double weightedDegree = weightedOutDegrees[u];

    index targetCluster;

    double weightChangeCurrent = 0;

    G->forEdgesOf(u, [&](node, node v, edgeweight weight) {
        if (u != v) {
            targetCluster = result[v];
            if (currentCluster != targetCluster) {
                if (!neighborClusterWeights.indexIsUsed(targetCluster))
                    neighborClusterWeights.insert(targetCluster, 0);
                // outgoing to other cluster from u wont affect cut of current cluster anymore
                weightChangeCurrent -= frequency * weight / weightedDegree;

                // sum up portion of outgoing edges to new target cluster
                neighborClusterWeights[targetCluster] += weight;
            }
        }
    });

    // calculate change of target cluster cut caused by outgoing edges as 
    // relative weighted frequency of all outgoing edges of u minus outgoing edges to target cluster
    if (neighborClusterWeights.size() > 0) {
        neighborClusterWeights.forElements([&](index targetCluster,
                                               double weightDelta) {
            neighborClusterWeights[targetCluster] = (weightedDegree - weightDelta) / weightedDegree * frequency;
        });
    } 

    G->forInEdgesOf(u, [&](node, node v, edgeweight weight) {
        if (u != v) {
            targetCluster = result[v];
            if (currentCluster != targetCluster) {
                if (!neighborClusterWeights.indexIsUsed(targetCluster)) {
                    // if no outgoing edge from u to target cluster, relative weighted frequency to other clusters equals node frequency of u
                    neighborClusterWeights.insert(targetCluster, frequency);
                }

                neighborClusterWeights[targetCluster] -= nodeFrequencies[v] * weight / weightedOutDegrees[v];

            } else {
                // edges to u in current cluster are now outgoing edges
                weightChangeCurrent += nodeFrequencies[v] * weight / weightedOutDegrees[v];
            }
        }
    });

    double bestChange, change;
    index bestTarget = currentCluster;
    double bestCutDelta;
    double currentVolumeNew = clusterVolume[currentCluster] - frequency;
    double targetCutDelta;


    const auto subsetSizeMap = result.subsetSizeMap(); // sizes of clusters
    const count numberNodes = G->numberOfNodes();   

    // change for current cluster cut if u moves
    const double currentCutDelta = tau * (currentVolumeNew / (double) numberNodes - (numberNodes - subsetSizeMap.at(currentCluster)) / (double) numberNodes * frequency) + (1 - tau) * weightChangeCurrent;
    
    // determine best target cluster
    if (neighborClusterWeights.size() > 0) {

        bestChange = 0; // do not move u 

        neighborClusterWeights.forElements([&](index targetCluster,
                                               double weightDelta) {
            targetCutDelta = tau * ((numberNodes - subsetSizeMap.at(targetCluster) - 1) / (double) numberNodes * frequency - clusterVolume[targetCluster] / numberNodes) + (1 - tau) * weightDelta;

            change = fitnessChange(u, frequency, currentCutDelta, targetCutDelta, currentCluster, targetCluster);

            if (change < bestChange || (change == bestChange && targetCluster < bestTarget && bestTarget != currentCluster)) {
                bestChange = change;
                bestTarget = targetCluster;
                bestCutDelta = targetCutDelta;
            }
        });
    }

    neighborClusterWeights.reset();

    // move u to best target cluster
    if (bestTarget != currentCluster) {

#ifndef NDEBUG
        oldfitness = mapEquation();
        lastChange = bestChange;
#endif
        bool res = performMove(u, frequency, currentCluster, bestTarget, currentCutDelta, bestCutDelta, bestChange);
#ifndef NDEBUG
        newfitness = mapEquation();
        assert(abs(newfitness - oldfitness - lastChange) < pow(10, -12));
#endif  
        return res;
    }

    return false;
}

double DirectedLouvainMapEquation::fitnessChange(node u, double frequency, double currentDelta, double targetDelta, node currentCluster, node targetCluster) {
    const double currentCutOld = clusterCut[currentCluster];
    const double targetCutOld = clusterCut[targetCluster];

    const double currentVolumeNew = clusterVolume[currentCluster] - frequency;
    const double targetVolumeNew = clusterVolume[targetCluster] + frequency;

    // calculate deltas for three terms in map equation that change when moving node u
    const double totalCutDelta = plogp(totalCut + currentDelta + targetDelta) - plogp(totalCut);
    const double clusterCutDelta = 2 * (plogp(currentCutOld) + plogp(targetCutOld) - plogp(currentCutOld + currentDelta) - plogp(targetCutOld + targetDelta));
    const double cutPlusVolumeDelta = plogp(currentCutOld + currentDelta + currentVolumeNew) + plogp(targetCutOld + targetDelta + targetVolumeNew) - plogp(currentCutOld + clusterVolume[currentCluster]) - plogp(targetCutOld + clusterVolume[targetCluster]);

    return totalCutDelta + clusterCutDelta + cutPlusVolumeDelta;
}


bool DirectedLouvainMapEquation::performMove(node u, double frequency, node currentCluster, node targetCluster, 
                                             double currentCutDelta, double targetCutDelta, double fitnessDelta) {
    bool moved = true;

    // node moves from current to target -> adjust sum of frequencies in both clusters
    clusterVolume[currentCluster] -= frequency;
    clusterVolume[targetCluster] += frequency;

    // adjust sum of all q_i
    totalCut += currentCutDelta + targetCutDelta;

    // set new q_i for current and target cluster
    clusterCut[currentCluster] += currentCutDelta;
    clusterCut[targetCluster] += targetCutDelta;

    // adjust result partition with move
    result.moveToSubset(targetCluster, u);

    fitness += fitnessDelta;

    return moved;
}


// calculate node visit frequencies using page rank algorithm
void DirectedLouvainMapEquation::calculateNodeFrequencies() {
    // initial distribution of uniform frequencies
    std::vector<double> nodeFrequenciesNew(G->upperNodeIdBound(), 1 / (double) G->numberOfNodes());

    int iter = 0;

    double distributeNeighborPortion, distributeGraphPortion;

    while (computeAbsoluteDifference(nodeFrequencies, nodeFrequenciesNew) > pow(10, -12)) {
        // store current iterative in nodeFrequencies and reset nodeFrequenciesNew to zeros
        std::copy(nodeFrequenciesNew.begin(), nodeFrequenciesNew.end(), nodeFrequencies.begin());
        std::fill(nodeFrequenciesNew.begin(), nodeFrequenciesNew.end(), 0);

        for (node u = 0; u < G->upperNodeIdBound(); ++u) {
            if (G->hasNode(u)) {
                distributeNeighborPortion = (1 - tau) * nodeFrequencies[u];
                
                distributeGraphPortion = nodeFrequencies[u] / (double) G->numberOfNodes();
                // if no outgoing edges set tau to 1 otherwise scale with tau
                if (G->degreeOut(u) > 0) {
                    distributeGraphPortion *= tau;
                }

                // distribute portion of frequency to neighbors
                G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                     nodeFrequenciesNew[v] += distributeNeighborPortion * ew / weightedOutDegrees[u];
                });

                // distribute rest uniformly to whole graph
                for (node x = 0; x < G->upperNodeIdBound(); x++) {
                    if (G->hasNode(x))
                        nodeFrequenciesNew[x] += distributeGraphPortion;
                }
            }
        }

        if (iter++ > 500) {
            throw std::runtime_error("Power iteration for calculating node frequencies needed too many iterations to finish!");
        }
    }

    // final copy to store last iterative in nodeFrequencies
    std::copy(nodeFrequenciesNew.begin(), nodeFrequenciesNew.end(), nodeFrequencies.begin());

#ifndef NDEBUG
    double sum = 0;
    for (int i=0; i< nodeFrequencies.size();i++) {
        sum+= nodeFrequencies[i];
    }
    assert(abs(sum - 1)<pow(10,-12)) ;
#endif
}

double DirectedLouvainMapEquation::computeAbsoluteDifference(std::vector<double> v, std::vector<double> w) {
    double sum = 0;

    assert(v.size() == w.size());

    for (int i = 0; i < v.size(); i++) {
        sum += abs(v[i] - w[i]);
    }

    return sum;
}

void DirectedLouvainMapEquation::calculateInitialClusterCutAndVolume() {
    totalCut = 0.0;
    totalVolume = 0.0;


    // sums for calculation of initial value of map equation
    double sumCut = 0; // sums plogp(q_i) over all clusters i
    double sumCutVolume = 0; // sums plogp(q_i + volume_i) over all clusters i

    int numberNodes = G->numberOfNodes();

    for (node u = 0; u < G->upperNodeIdBound(); ++u) {
        clusterCut[u] = 0;
        if (G->hasNode(u)) {
            G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                // ignore loops
                if (u != v)
                    clusterCut[u] += ew * nodeFrequencies[u];
            });
        }
        // clusterCut contains q_i for cluster i (here cluster is one node)
        clusterCut[u] *= 1 - tau;
        clusterCut[u] /= weightedOutDegrees[u];
        clusterCut[u] += tau * (numberNodes - 1) / numberNodes * nodeFrequencies[u];
        // clusterVolume contains sum of node frequencies of nodes in cluster
        clusterVolume[u] = nodeFrequencies[u];

        // total cut is sum of q_i over all clusters; volume is sum of all node frequencies
        totalCut += clusterCut[u];
        totalVolume += clusterVolume[u];

        sumCut += plogp(clusterCut[u]);
        sumCutVolume += plogp(clusterCut[u] + clusterVolume[u]);
    }
    double sumFrequencies = 0;
    for (node u = 0; u<G->upperNodeIdBound(); u++) {
        if(G->hasNode(u))
            sumFrequencies += plogp(nodeFrequencies[u]);
    }

    fitness = plogp(totalCut) - 2 * sumCut + sumCutVolume - sumFrequencies;
#ifndef NDEBUG
    assert(abs(fitness - mapEquation())< pow(10, -12));
#endif
}

std::string DirectedLouvainMapEquation::toString() const {
    return "DirectedLouvainMapEquation";
}


double DirectedLouvainMapEquation::plogp(double w) {
    if (w > 0) {
        return w * log(w);
    }
    return 0;
}

#ifndef NDEBUG

long double DirectedLouvainMapEquation::mapEquation() {

    std::set<std::set<index>> subsets = result.getSubsets();

    std::vector<double> q(subsets.size());
    std::vector<double> p(subsets.size()); 
    double sumQ = 0;

    int idx = 0;
    for (auto subset: subsets) {
        q[idx] = 0;
        p[idx] = 0;
        for (std::set<index>::iterator it = subset.begin(); it != subset.end(); it++) {
            node u = *it;
            if (G->hasNode(u)) {
                p[idx] += nodeFrequencies[u];

                G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                    if (result[u] != result[v])
                        q[idx] += ew * nodeFrequencies[u] / weightedOutDegrees[u];
                });
            }
        }

        q[idx] *= 1 - tau;
        q[idx] += tau * (G->numberOfNodes() - subset.size()) / G->numberOfNodes() * p[idx];
        sumQ += q[idx];
        idx++;
    }


    double sumPLogPClusterCutPlusVol = 0;
    double sumPLogPClusterCut = 0;

    for (index i = 0; i < q.size(); ++i) {
        sumPLogPClusterCutPlusVol += plogp(q[i] + p[i]);
        sumPLogPClusterCut += plogp(q[i]);
    }

    double sumPLogPwAlpha = 0;
    for (node u = 0; u<G->upperNodeIdBound(); u++) {
        if(G->hasNode(u))
            sumPLogPwAlpha += plogp(nodeFrequencies[u]);
    }

    return plogp(sumQ) - 2 * sumPLogPClusterCut + sumPLogPClusterCutPlusVol - sumPLogPwAlpha;
}
#endif

} // namespace NetworKit

