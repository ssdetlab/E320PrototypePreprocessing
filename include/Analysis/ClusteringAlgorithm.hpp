#pragma once

#include <cstddef>

#include "Types/TypeDefinitions.hpp"

/// @brief algorithm performing BFS clustering 
class ClusteringAlgorithm {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Maximum cluster size
    std::size_t maxClusterSize;
  };

  ClusteringAlgorithm() = delete;

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit ClusteringAlgorithm(const Config &cfg) : m_cfg(cfg) {}

  /// @brief Build clusters from a hit collection
  ///
  /// @param hits collection of hits
  ///
  /// @return vector of clusters
  std::vector<Cluster> getClusters(HitSet hits) const;

 private:
  /// Configuration
  Config m_cfg;

  /// @brief Perform BFS clustering on the hit collection
  ///
  /// @param clusterPixels collection of hits in a found cluster
  /// @param pivot cluster to start the search from
  /// @param pixels hit collection to cluster
  void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                     HitSet &pixels) const;
};
