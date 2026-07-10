#pragma once

#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "ROI/ROI.hpp"
#include "TF1.h"
#include "TGraph.h"
#include "TH2.h"

/// @brief Buttrfly-shaped ROI
class ButterflyROI : public ROI {
 public:
  /// @brief constructor
  ///
  /// @param x0 butterfly center x
  /// @param y0 butterfly center y
  /// @param w butterfly minor axis size
  /// @param theta butterfly rotation angle
  /// @param curv butterfly curvature
  ButterflyROI(double x0, double y0, double w, double theta, double curv);

  /// @brief check if the point is inside the ROI
  ///
  /// @param x point X
  /// @param y point Y
  ///
  /// @return true if (x,y) is inside the ROI, false otherwise
  bool isInside(double x, double y) const override;

  /// Parameters getters
  double centerX() const { return m_x0; };
  double centerY() const { return m_y0; };
  double width() const { return m_w; };
  double theta() const { return m_theta; };
  double curvature() const { return m_curv; };

 private:
  /// Butterfly center x
  double m_x0;
  /// Butterfly center y
  double m_y0;
  /// Butterfly minor axis size
  double m_w;
  /// Butterfly rotation angle
  double m_theta;
  /// Butterfly curvature
  double m_curv;
};

/// @brief class handling the butterfly ROI construction
class ButterflyROIFactory {
 public:
  using ROIContainer = std::map<std::size_t, std::shared_ptr<ButterflyROI>>;

  ButterflyROIFactory() = default;
  ~ButterflyROIFactory() = default;

  /// @brief Output diagnostic plots struct
  struct Stats {
    /// Raw hit map
    TH2D rawHist;
    /// Raw bin counts
    TH1D rawBinCounts;
    /// Bin counts with quantile cut applied
    TH1D trimmedBinCounts;
    /// Butterfly slope graph
    TGraph slopeGr;
    /// Butterfly slope fit function
    TF1 slopeFitFunc;
    /// Decorrelated hit map
    TH2D decorrHist;
    /// Collection of 1D slices of
    /// the decorrelated hit map
    std::vector<std::tuple<double, TH1D, TF1>> slices;
    /// BUtterfly curvature graph
    TGraph curvGr;
    /// BUtterfly curvature fit function
    TF1 curvFitFunc;
    /// Hit map with the ROI applied
    TH2D testHist;
  };

  /// @brief ROI calculation options struct
  struct Options {
    /// Number of steps along the butterfly
    /// taken to  construct the slope, curvature,
    /// and slices in number of bins
    int nSteps;
    /// Width of the cummulation tunnel around
    /// the step in number of bins
    int stepWidth;
    /// Size of the step in number of bins
    int stepSize;
    /// Number of std devs defining the butterfly width
    double nStdDevsWidth;
    /// Application of the noisy pixel filter flag
    bool applyNoiseFilter;
    /// Noisy pixel low quantile
    double noiseQuantileLow;
    /// Noisy pixel high quantile
    double noiseQuantileHigh;
    /// Noisy pixel high quantile
    double nStdDevsNoise;
    /// Chip IDs to analyze
    std::set<std::size_t> chipIds;
  };

  /// @brief Construct the butterfly for a spefic run data
  ///
  /// @param pathToDir directory with the run data
  /// @param treeName name of the data tree
  /// @param prefix unique run identifier
  /// @opt Options struct
  ///
  /// @return map of chip IDs to ROIs for the run
  ROIContainer buildROIForRun(const std::string &pathToDir,
                              const std::string &treeName,
                              const std::string &prefix,
                              const Options &opt) const;

  /// @brief Read the butterfly for a spefic run data
  ///
  /// @param pathToCsv CSV file containing the parameters
  ///
  /// @return map of chip IDs to ROIs for the run
  ROIContainer buildROIForRun(const std::string &pathToCsv) const;

  /// @brief Dump run ROI collection into a CSV file
  ///
  /// @param path CSV file containing the parameters
  /// @param rois ROIs to dump
  void dumpToCsvForRun(const std::string &path, const ROIContainer &rois) const;
};
