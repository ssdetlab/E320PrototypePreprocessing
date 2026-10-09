#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Eigen>

#include "AcceleratorState.hpp"
#include "Analysis/PathCollection.hpp"
#include "DetectorEvent.h"
#include "ROI/ROI.hpp"
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

  /// @brief supported ROI types enumeration
  enum ROIType : int { Ellipse = 0, Butterfly = 1, Unbounded = 2 };

  /// @brief preprocessing options struct
  struct Options {
    /// Data directory path
    std::string inDataDir;
    /// ROI CSV data directory path
    std::string inROIDir;
    /// Data tree name
    std::string inDataTreeName;
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
    /// Accelerator data noise low quantile
    double noiseQuantileLowAcc;
    /// Accelerator data noise high quantile
    double noiseQuantileHighAcc;
    /// Max number of std devs for the
    /// accelerator data
    int nStdDevsAcc;
    /// ROI type
    ROIType roiType;
    /// Maximum number of hits
    std::size_t maxHits;
    /// Maximum number of hits in the ROI
    std::size_t maxRoiHits;
    /// Maximum cluster size
    std::size_t maxClusterSize;
    /// Exclude shuttered events flag
    bool excludeShuttered;
  };

  /// @brief Output diagnostic plots struct
  struct Stats {
    Stats() = default;

    /// @brief initialize the diagnostic plots
    ///
    /// @param prefix unique run identifier
    /// @param ids chip IDs under analysis
    void init(const std::string &prefix, const std::set<std::size_t> &ids);

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

    /// Raw PMTs
    TGraph gr_PMT_S20_3060_Raw;
    TGraph gr_PMT_S20_3070_Raw;
    TGraph gr_PMT_S20_3179_Raw;
    TGraph gr_PMT_S20_3350_Raw;
    TGraph gr_PMT_S20_3360_Raw;

    /// Filtered PMTs
    TGraph gr_PMT_S20_3060_Filtered;
    TGraph gr_PMT_S20_3070_Filtered;
    TGraph gr_PMT_S20_3179_Filtered;
    TGraph gr_PMT_S20_3350_Filtered;
    TGraph gr_PMT_S20_3360_Filtered;

    /// Raw toroids
    TGraph gr_TORO_LI20_1988_TMIT_PC_Raw;
    TGraph gr_TORO_LI20_2040_TMIT_PC_Raw;
    TGraph gr_TORO_LI20_2452_TMIT_PC_Raw;
    TGraph gr_TORO_LI20_3163_TMIT_PC_Raw;
    TGraph gr_TORO_LI20_3255_TMIT_PC_Raw;

    /// Filtered toroids
    TGraph gr_TORO_LI20_1988_TMIT_PC_Filtered;
    TGraph gr_TORO_LI20_2040_TMIT_PC_Filtered;
    TGraph gr_TORO_LI20_2452_TMIT_PC_Filtered;
    TGraph gr_TORO_LI20_3163_TMIT_PC_Filtered;
    TGraph gr_TORO_LI20_3255_TMIT_PC_Filtered;
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

  /// @brief get file handle for a given ROOT tree
  ///
  /// @param pathToROIData path to directory with ROI CSV files
  /// @param opt option struct
  ///
  /// @return chip ID to ROI map for the run
  std::map<std::size_t, std::shared_ptr<ROI>> getROIs(
      const std::string &pathToROIData, const Options &opt) const;

  /// @brief perform BFS clustering on the collection of hits
  ///
  /// @param clusterPixels collection cluster hits
  /// @param pivot hit to start the search from
  /// @param pixels collection of active pixels to search the clusters in
  void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                     HitSet &pixels) const;

  /// @brief find active pixel clusters in the collection of hits
  ///
  /// @param hits active pixels collection
  /// @param maxCluterSize maximum allowed size of a pixel cluster
  ///
  /// @return vector of active pixel clusters
  std::vector<Cluster> getClusters(HitSet &hits,
                                   std::size_t maxCluterSize) const;

  /// @brief get accelerator state means and std devs for a run
  ///
  /// @param tree tree with the accelerator data
  /// @param skip number of BXs to skip in the beginning
  /// @param quantileLow low quantile trim for the diagnostic distributions
  /// @param quantileHigh high quantile trim for the diagnostic distributions
  ///
  /// @return means and std devs of the accelerator diagnostics
  std::pair<AcceleratorState, AcceleratorState> getAcceleratorStateMeanStd(
      TTree *tree, std::size_t skip, const Options &opt) const;

  /// @brief tag anomalous BXs by analysing the accelerator state
  ///
  /// @param tree tree with the accelerator data
  /// @param skip number of BXs to skip in the beginning
  /// @param nStdDevs number of the diagnostics std devs to consider a BX anomalous
  ///
  /// @return collection of "bad" BXs' tree entry numbers
  std::unordered_set<std::size_t> getBadEvents(TTree *tree, std::size_t skip,
                                               const Options &opt) const;

  /// @brief convert EPICS dipole PV to T
  ///
  /// @param val PV value -- beam energy in eV
  /// @param theta deflection angle in rad
  /// @param bz dipole length in m
  ///
  /// @return field strength in T
  double convertEpicsDipole(double val, double theta, double bz) const;

  /// @brief convert EPICS x-corrector PV to T
  ///
  /// @param val PV value -- kG * m
  /// @param bz dipole length in m
  ///
  /// @return field strength in T
  double convertEpicsXCorrector(double val, double bz) const;

  /// @brief convert EPICS quadrupole PV to T
  ///
  /// @param val PV value -- kG
  ///
  /// @return field gradient in T/m
  double convertEpicsQuad(double val) const;

  /// @brief extract EPICS PV value as a double
  ///
  /// @param detEvent detector event pointer
  /// @param pvname name of the PV to extract
  ///
  /// @return value of the PV
  double getPv(const detector_event *detEvent, const std::string &pvname) const;
};
