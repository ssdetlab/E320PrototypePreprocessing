#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Eigen>

#include "AcceleratorState.hpp"
#include "Analysis/ROI.hpp"
#include "Analysis/Stats.hpp"
#include "DetectorEvent.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "Types/DetectorEvent.hpp"
#include "detail/PairHash.hpp"

class Preprocessor {
 public:
  using FileHandle = std::tuple<TFile *, TChain *, TTree *>;

  using Cluster = E320Io::ChipEvent::Cluster;
  using Hit = std::pair<std::uint16_t, std::uint16_t>;
  using HitSet = std::unordered_set<Hit, PairHash>;

  enum ROIType : int { Ellipse = 0, Butterfly = 1, Unbounded = 2 };

  struct Options {
    /// Data directory path
    std::string pathToDataDir;
    /// ROI data directory path
    std::string pathToROIDir;
    /// Data tree name
    std::string dataTreeName;
    /// ROI data tree name
    std::string roiDataTreeName;
    /// Output directory path
    std::string outPathDir;
    /// Output data tree name
    std::string outDataTreeName;
    /// Output magnet tree name
    std::string outMagnetTreeName;
    /// Chip IDs of the dataset
    std::set<int> chipIds;
    /// Magnets parameters
    double dipoleLength;      // in meters
    double xCorrectorLength;  // in meters
    double deflectionAngle;   // in rad
    /// Max number of StdDevs for the
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
  };

  struct PathCollection {
    /// Prefix of the run files
    std::string runPrefix = "";
    /// Collection of the run files
    std::vector<std::string> inPaths{};
    /// Output path
    std::string outPath = "";
    /// Number of events to skip in the beginning
    std::size_t skip = 0;
    /// Number of events to skip in the end
    std::size_t stop = 0;
  };

  /// Run data preprocessing
  void processFiles(const Options &opt);

 private:
  ///
  std::map<std::string, PathCollection> collectDataPaths(const Options &opt);

  std::unordered_map<int, std::shared_ptr<ROI>> computeROIs(const Options &opt,
                                                            Stats &stats);

  FileHandle getFileHandle(const std::string &treeName,
                           const std::vector<std::string> &paths);

  void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                     HitSet &pixels);

  std::vector<Cluster> getClusters(HitSet &hits, std::size_t maxCluterSize);

  std::pair<AcceleratorState, AcceleratorState> getAcceleratorStateMeanStd(
      TTree *tree, std::size_t skip);

  std::unordered_set<std::size_t> getBadEvents(TTree *tree, std::size_t skip,
                                               int nStdDevs);

  // @param val PV value -- beam energy in eV
  // @param theta deflection angle in rad
  // @param bz dipole length in m
  //
  // @return field strength in T
  double convertEpicsDipole(double val, double theta, double bz);

  // @param val PV value -- kG * m
  // @param bz dipole length in m
  //
  // @return field strength in T
  double convertEpicsXCorrector(double val, double bz);

  // @param val PV value -- kG * m
  //
  // @return field strength in T
  double convertEpicsQuad(double val);

  double getPv(const detector_event *detEvent, const std::string &pvname);
};
