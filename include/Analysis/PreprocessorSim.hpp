#pragma once

#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"

/// @brief class handling the E320 data preprocessing
class PreprocessorSim {
 public:
  /// Shorthand for the ROOT file handles
  using FileHandle = std::tuple<TFile*, TChain*, TTree*>;

  /// @brief preprocessing options struct
  struct Options {
    /// Tracker data directory path
    std::string inTrackerDir;
    /// Beam orbit data directory path
    std::string inBeamDir;
    /// Tracker data tree name
    std::string inTrackerTreeName;
    /// Beam orbit data tree name
    std::string inBeamTreeName;
    /// Output directory path
    std::string outPath;
    /// Output tree name
    std::string outTreeName;
  };

  /// @brief Run data preprocessing
  ///
  /// @param opt preprocessing options
  void processFiles(const Options& opt) const;

 private:
  /// @brief get file handle for a given ROOT tree
  ///
  /// @param treeName name of the data tree
  /// @param paths run data paths
  ///
  /// @return file handle for run data reading
  FileHandle getFileHandle(const std::string& treeName,
                           const std::vector<std::string>& paths) const;

 protected:
};
