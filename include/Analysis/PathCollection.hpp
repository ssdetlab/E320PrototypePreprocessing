#pragma once

#include <string>
#include <vector>

/// @brief run data paths container
struct PathCollection {
  /// Prefix of the run files
  std::string runPrefix = "";
  /// Parent path of the run files
  std::string parentPath = "";
  /// Collection of the run files
  std::vector<std::string> inPaths{};
  /// Collection of the run orbit reconstruction files
  std::vector<std::string> inOrbitPaths{};
  /// Output path
  std::string outPath = "";
  /// Number of events to skip in the beginning
  std::size_t skip = 0;
  /// Tree entry number to stop at (exclusive),
  /// the run entries are [skip, stop)
  std::size_t stop = 0;
};
