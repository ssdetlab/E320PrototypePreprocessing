#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DetectorEvent.h"

/// @brief algorithm tagging anomalous BXs by analysing the accelerator state
class AcceleratorStateFilter {
 public:
  /// @brief supported PV cut types enumeration
  enum CutType : int { Upper = 0, Central = 1 };

  /// @brief accelerator state PV configuration
  struct PvConfig {
    /// PV name in the EPICS frame
    std::string key;
    /// Label of the PV in the histograms and diagnostic plots
    std::string name;
    /// Cut applied to the PV
    CutType cutType;
    /// Number of bins of the PV distribution
    int nBins;
    /// Lower edge of the PV distribution
    double low;
    /// Upper edge of the PV distribution
    double high;
  };

  /// @brief accelerator data of a BX
  struct BxData {
    /// Tree entry number
    std::size_t entry = 0;
    /// Trigger number
    std::uint64_t trgN = 0;
    /// PID PV presence flag
    bool hasPid = false;
    /// Shuttered BX flag
    bool isShuttered = false;
    /// Values of the accelerator state PVs, in the
    /// configuration order up to the first missing one
    std::vector<double> values;
    /// All accelerator state PVs presence flag
    bool hasPvs = false;
  };

  /// @brief Nested configuration struct
  struct Config {
    /// Accelerator data noise low quantile
    double noiseQuantileLow;
    /// Accelerator data noise high quantile
    double noiseQuantileHigh;
    /// Max number of std devs for the
    /// accelerator data
    int nStdDevs;
    /// Accelerator state PVs, read and cut in this order
    std::vector<PvConfig> pvs;
  };

  AcceleratorStateFilter() = delete;

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit AcceleratorStateFilter(const Config &cfg) : m_cfg(cfg) {}

  /// @brief Read the accelerator data of a BX
  ///
  /// @param detEvent detector event pointer
  /// @param entry tree entry number of the BX
  ///
  /// @return accelerator data of the BX
  BxData getBxData(const detector_event *detEvent, std::size_t entry) const;

  /// @brief Tag anomalous BXs by analysing the accelerator state
  ///
  /// @param bxs accelerator data of the run BXs
  ///
  /// @return collection of "bad" BXs' tree entry numbers
  std::unordered_set<std::size_t> getBadEvents(
      const std::vector<BxData> &bxs) const;

  /// @brief Get the accelerator state PV configurations
  ///
  /// @return PV configurations in the reading order
  const std::vector<PvConfig> &getPvs() const { return m_cfg.pvs; }

 private:
  /// Configuration
  Config m_cfg;

  /// @brief Get the accelerator state means and std devs
  ///
  /// @param bxs accelerator data of the run BXs
  ///
  /// @return means and std devs of the PVs in the configuration order
  std::pair<std::vector<double>, std::vector<double>>
  getAcceleratorStateMeanStd(const std::vector<BxData> &bxs) const;
};
