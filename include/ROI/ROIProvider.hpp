#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "ROI/ROI.hpp"

/// @brief class providing the per-chip ROIs of the runs
class ROIProvider {
 public:
  /// Chip ID to ROI map of a run
  using ROIContainer = std::map<std::size_t, std::shared_ptr<ROI>>;

  /// @brief supported ROI types enumeration
  enum ROIType : int { Ellipse = 0, Butterfly = 1, Unbounded = 2 };

  /// @brief Nested configuration struct
  struct Config {
    /// ROI CSV data directory path
    std::string inROIDir;
    /// ROI type
    ROIType roiType;
    /// Chip IDs of the dataset
    std::set<std::size_t> chipIds;
    /// Run prefix to ROI run prefix map, the ROIs of a run
    /// are read from <inROIDir>/<ROI run prefix>.csv
    /// (e.g. ROIs estimated on a signal run applied
    /// to a background run)
    std::map<std::string, std::string> roiRunMap;
  };

  ROIProvider() = delete;

  /// @brief Constructor, reads the ROIs
  /// of the runs in the ROI run map
  ///
  /// @param cfg configuration struct
  explicit ROIProvider(const Config &cfg);

  /// @brief Get the ROIs of a run
  ///
  /// @param runPrefix unique run identifier
  ///
  /// @return chip ID to ROI map for the run
  const ROIContainer &getROIs(const std::string &runPrefix) const;

 private:
  /// Configuration
  Config m_cfg;

  /// Run prefix to ROI collection map
  std::map<std::string, ROIContainer> m_rois;

  /// ROI collection shared by all the runs
  /// for the unbounded ROI type
  ROIContainer m_unboundedROIs;
};
