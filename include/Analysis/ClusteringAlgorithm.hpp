#pragma once

#include <cstddef>

#include "Types/DetectorEvent.hpp"
#include "Types/TypeDefinitions.hpp"

class ClusteringAlgorithm {
 public:
  /// Hit-related data types
  using Cluster = E320Io::ChipEvent::Cluster;

  struct Config {
    /// Maximum cluster size
    std::size_t maxClusterSize;
  };

  ClusteringAlgorithm() = delete;

  explicit ClusteringAlgorithm(const Config &cfg) : m_cfg(cfg) {}

  std::vector<Cluster> getClusters(HitSet &hits) const;

 private:
  Config m_cfg;

  void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                     HitSet &pixels) const;
};
