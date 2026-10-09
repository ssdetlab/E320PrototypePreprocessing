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
  /// Output path
  std::string outPath = "";
  /// Number of events to skip in the beginning
  std::size_t skip = 0;
  /// Number of events to skip in the end
  std::size_t stop = 0;
};
