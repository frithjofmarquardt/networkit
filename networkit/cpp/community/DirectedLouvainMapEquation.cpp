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
      weightedOutDegrees(G->upperNodeIdBound(), 0) {
        tau = 0.15;
        fitness = -1;;
      }

void DirectedLouvainMapEquation::run() {

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
        
        numberOfNodesMoved = localMoving(nodes);
    }
    

    hasRun = true;
}

count DirectedLouvainMapEquation::localMoving(std::vector<node> &nodes) {
    count nodesMoved = 0;

    for (node u : nodes) {
        if (tryLocalMove(u)) {
            nodesMoved += 1;
        }
#ifndef NDEBUG
        newfitness = mapEquation();
        assert(abs(newfitness - oldfitness - lastChange) < pow(10, -12));
#endif
    }

    return nodesMoved;
}

// try to move node u to clusters of neighboring nodes (all in and outgoing neighbors)
bool DirectedLouvainMapEquation::tryLocalMove(node u) {
    
    const index currentCluster = result[u];
    index targetCluster;
    const double frequency = nodeFrequencies[u];
    double newCurrentCut, newTargetCut;
    double newCurrentVolume, newTargetVolume;

    double change;
    double bestChange;
    index bestTarget = currentCluster;
    double bestCurrentCut, bestTargetCut;

    // for all neighbors of u, test if moving u to neighbors cluster improves map equation
    auto lambda = [&](node, node v, edgeweight weight) {
        targetCluster = result[v];
        if (currentCluster != targetCluster && targetCluster != bestTarget) {
            calculateNewCutAndVolume(u, frequency, currentCluster, targetCluster, 
                newCurrentVolume, newTargetVolume, newCurrentCut, newTargetCut);
            change = fitnessChange(u, frequency, newCurrentVolume, newTargetVolume, 
                newCurrentCut, newTargetCut, currentCluster, targetCluster);

            if (!bestChange || change < bestChange) {
                bestTarget = targetCluster;
                bestCurrentCut = newCurrentCut;
                bestTargetCut = newTargetCut;
                bestChange = change;
            }
        }
    };

    G->forEdgesOf(u, lambda);
    G->forInEdgesOf(u, lambda);

    // move u to best target cluster
    if (bestTarget != currentCluster) {

#ifndef NDEBUG
        oldfitness = mapEquation();
        lastChange = bestChange;
#endif
        return performMove(u, frequency, currentCluster, bestTarget, bestCurrentCut, bestTargetCut, bestChange);


    }

    return false;
}

double DirectedLouvainMapEquation::fitnessChange(node u, double frequency, double newCurrentVolume, double newTargetVolume, 
                                                 double newCurrentCut, double newTargetCut, node currentCluster, node targetCluster) {
    
    const double currentCutOld = clusterCut[currentCluster];
    const double targetCutOld = clusterCut[targetCluster];
    
    // calculate deltas for three terms in map equation that change when moving node u
    const double totalCutDelta = plogp(totalCut - currentCutOld - targetCutOld + newCurrentCut + newTargetCut) - plogp(totalCut);
    double clusterCutDelta = 2 * (plogp(currentCutOld) + plogp(targetCutOld) - plogp(newTargetCut) - plogp(newCurrentCut));
    double cutPlusVolumeDelta = plogp(newCurrentCut + newCurrentVolume) + plogp(newTargetCut + newTargetVolume) - plogp(currentCutOld + clusterVolume[currentCluster]) - plogp(targetCutOld + clusterVolume[targetCluster]);

    return totalCutDelta + clusterCutDelta + cutPlusVolumeDelta;
}

void DirectedLouvainMapEquation::calculateNewCutAndVolume(node u, double frequency, node currentCluster, 
                                                          node targetCluster, double &newCurrentVolume, double &newTargetVolume,
                                                          double &newCurrentCut, double &newTargetCut) {
    const count numberNodes = G->numberOfNodes();

    const auto subsetSizeMap = result.subsetSizeMap();

    double sumOut = 0;
    double sumIn = 0;

    if (subsetSizeMap.at(currentCluster) == 1) {
        // if u is last node in cluster set q_i and volume to 0
        newCurrentVolume = 0;
        newCurrentCut = 0;
    } else {
        newCurrentVolume = clusterVolume[currentCluster] - frequency;

        G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
            if (result[v] != currentCluster)
                sumOut += frequency * ew;
        });

        sumOut /= weightedOutDegrees[u];

        G->forInEdgesOf(u, [&](node, node v, edgeweight ew) {
            if (result[v] == currentCluster)
                sumIn += nodeFrequencies[v] * ew / weightedOutDegrees[v];
        });

        // newCut = oldCut + tau*(newVolume/n - (n-n_i)*frequency[u]/n) + (1-tau)*(sum_in-sum_out)
        newCurrentCut = clusterCut[currentCluster] + tau * (newCurrentVolume / (double) numberNodes - (numberNodes - subsetSizeMap.at(currentCluster)) / (double) numberNodes * frequency) + (1 - tau) * (sumIn - sumOut);
    }
    
    if (subsetSizeMap.at(targetCluster) == numberNodes - 1) {
        // if nodes are all in 1 cluster set volume to 1, q_i to 0
        newTargetVolume = 1;
        newTargetCut = 0;
    } else {
        newTargetVolume = clusterVolume[targetCluster] + frequency;
        sumOut = 0;
        G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
            if (result[v] != targetCluster)
                sumOut += frequency * ew;
        });

        sumOut /= weightedOutDegrees[u];

        sumIn = 0;
        G->forInEdgesOf(u, [&](node, node v, edgeweight ew) {
            if (result[v] == targetCluster)
                sumIn += nodeFrequencies[v] * ew / weightedOutDegrees[v];
        });
        // newCut = oldCut + tau*(oldVolume/n - (n-n_i-1)*frequency[u]/n) + (1-tau)*(sum_out-sum_in)
        newTargetCut = clusterCut[targetCluster] + tau * ((numberNodes - subsetSizeMap.at(targetCluster) - 1) / (double) numberNodes * frequency - clusterVolume[targetCluster] / numberNodes) + (1 - tau) * (sumOut - sumIn);
    }
}

bool DirectedLouvainMapEquation::performMove(node u, double frequency, node currentCluster, node targetCluster, 
                                             double currentCut, double targetCut, double fitnessDelta) {
    bool moved = true;

    // node moves from current to target -> adjust sum of frequencies in both clusters
    clusterVolume[currentCluster] -= frequency;
    clusterVolume[targetCluster] += frequency;

    // adjust sum of all q_i
    totalCut += currentCut + targetCut - clusterCut[currentCluster] - clusterCut[targetCluster];

    // set new q_i for current and target cluster
    clusterCut[currentCluster] = currentCut;
    clusterCut[targetCluster] = targetCut;

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

