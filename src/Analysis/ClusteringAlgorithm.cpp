#include "Analysis/ClusteringAlgorithm.hpp"

#include <deque>
#include <limits>

std::vector<Cluster> ClusteringAlgorithm::getClusters(HitSet hits) const {
  std::vector<Cluster> clusters;

  while (!hits.empty()) {
    auto pixel = *hits.begin();
    std::vector<Hit> clusterPixels;
    bfsClustering(clusterPixels, pixel, hits);
    std::size_t clSize = clusterPixels.size();
    if (clSize > m_cfg.maxClusterSize) {
      continue;
    }

    std::size_t xMin = std::numeric_limits<std::size_t>::max();
    std::size_t xMax = 0;
    std::size_t yMin = std::numeric_limits<std::size_t>::max();
    std::size_t yMax = 0;
    double xCenter = 0;
    double yCenter = 0;
    for (auto pix : clusterPixels) {
      xCenter += pix.first;
      yCenter += pix.second;

      if (pix.first <= xMin) {
        xMin = pix.first;
      }
      if (pix.first >= xMax) {
        xMax = pix.first;
      }
      if (pix.second <= yMin) {
        yMin = pix.second;
      }
      if (pix.second >= yMax) {
        yMax = pix.second;
      }
    }
    std::size_t dx = xMax - xMin + 1;
    std::size_t dy = yMax - yMin + 1;
    xCenter /= clusterPixels.size();
    yCenter /= clusterPixels.size();

    clusters.push_back(
        Cluster{xCenter, yCenter, dx, dy, clSize, clusterPixels});
  }

  return clusters;
}

void ClusteringAlgorithm::bfsClustering(std::vector<Hit> &clusterPixels,
                                        Hit &pivot, HitSet &pixels) const {
  std::deque<Hit> queue{pivot};
  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop_front();
    auto currentPtr = pixels.find(current);

    if (currentPtr != pixels.end()) {
      // Add this pixel to the cluster
      clusterPixels.push_back(current);
      // Remove pixel from live pixels list
      pixels.erase(currentPtr);

      std::vector<Hit> neighbours{{current.first + 1, current.second},
                                  {current.first - 1, current.second},
                                  {current.first, current.second + 1},
                                  {current.first, current.second - 1}};
      for (const auto &n : neighbours) {
        auto nptr = pixels.find(n);
        if (nptr != pixels.end()) {
          queue.push_back(*nptr);
        }
      }
    }
  }
}
