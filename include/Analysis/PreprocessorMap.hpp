#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Eigen>

#include "Analysis/AcceleratorStateFilter.hpp"
#include "Analysis/ClusteringAlgorithm.hpp"
#include "Analysis/PathCollection.hpp"
#include "DetectorEvent.h"
#include "ROI/ROIProvider.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "Types/DetectorEvent.hpp"
#include "detail/PairHash.hpp"

/// @brief class handling the E320 data preprocessing
class Preprocessor {
 public:
  /// Shorthand for the ROOT file handles
  using FileHandle = std::tuple<TFile *, TChain *, TTree *>;

  /// Hit-related data types
  using Cluster = E320Io::ChipEvent::Cluster;
  using Hit = std::pair<std::uint16_t, std::uint16_t>;
  using HitSet = std::unordered_set<Hit, PairHash>;

  /// @brief preprocessing options struct
  struct Options {
    /// Data directory path
    std::string inDataDir;
    /// Orbit reconstruction data directory path
    std::string inOrbitDir;
    /// Data tree name
    std::string inDataTreeName;
    /// Orbit reconstruction tree name
    std::string inOrbitTreeName;
    /// Output directory path
    std::string outDir;
    /// Output data tree name
    std::string outDataTreeName;
    /// Output magnet tree name
    std::string outMagnetTreeName;
    /// Chip IDs of the dataset
    std::set<std::size_t> chipIds;
    /// Magnets parameters
    double dipoleLength;      /// in meters
    double xCorrectorLength;  /// in meters
    double deflectionAngle;   /// in rad
    /// Maximum number of hits
    std::size_t maxHits;
    /// Maximum number of hits in the ROI
    std::size_t maxRoiHits;
    /// Exclude shuttered events flag
    bool excludeShuttered;
    /// Clustering algo
    std::shared_ptr<ClusteringAlgorithm> clusteringAlgo;
    /// ROI provider
    std::shared_ptr<ROIProvider> roiProvider;
    /// Accelerator state filter
    std::shared_ptr<AcceleratorStateFilter> acceleratorStateFilter;
  };

  /// @brief Output diagnostic plots struct
  struct Stats {
    Stats() = default;

    /// @brief initialize the diagnostic plots
    ///
    /// @param prefix unique run identifier
    /// @param ids chip IDs under analysis
    /// @param pvNames accelerator state PV labels
    void init(const std::string &prefix, const std::set<std::size_t> &ids,
              const std::vector<std::string> &pvNames);

    /// @brief reset the diagnostics plots
    void reset();

    /// @brief dump the plots into a file
    void print();

    /// Chip ids
    std::set<std::size_t> chipIds;
    /// Run prefix
    std::string runPrefix = "";

    /// N total events
    std::size_t nEventsTotal = 0;
    /// N accepted events
    std::size_t nEventsAccepted = 0;
    /// N shuttered events
    std::size_t nEventsUnshuttered = 0;
    /// N even PID events
    std::size_t nEventsEvenPID = 0;
    /// N odd PID events
    std::size_t nEventsOddPID = 0;
    /// N events with full detector occupancy
    /// and even PID
    std::size_t nEventsFullDetEvenPID = 0;
    /// N events with full detector occupancy
    /// and odd PID
    std::size_t nEventsFullDetOddPID = 0;

    /// Event category hist
    TH1I eventCategotyHist;
    /// Raw hit map hists
    std::map<int, TH2D> chipHitMapRawHists;
    /// Filtered hit map hists
    std::map<int, TH2D> chipHitMapFilteredHists;
    /// Raw occupancy graphs
    std::map<int, TGraph> chipOccRawGrs;
    /// Filtered occupancy graphs
    std::map<int, TGraph> chipOccFilteredGrs;

    /// Accelerator state PV labels
    std::vector<std::string> accPvNames;
    /// Raw accelerator state PV graphs
    std::vector<TGraph> accRawGrs;
    /// Filtered accelerator state PV graphs
    std::vector<TGraph> accFilteredGrs;
  };

  /// @brief Run data preprocessing
  ///
  /// @param opt preprocessing options
  void processFiles(const Options &opt) const;

 private:
  /// @brief collect run data paths for a collection of runs
  ///
  /// @param opt preprocessing options
  ///
  /// @return map of run prefix to its path collection
  std::map<std::string, PathCollection> collectDataPaths(
      const Options &opt) const;

  /// @brief get file handle for a given ROOT tree
  ///
  /// @param treeName name of the data tree
  /// @param paths run data paths
  ///
  /// @return file handle for run data reading
  FileHandle getFileHandle(const std::string &treeName,
                           const std::vector<std::string> &paths) const;
};
