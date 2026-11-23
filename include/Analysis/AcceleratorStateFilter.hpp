#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <utility>
#include <vector>

#include "AcceleratorState.hpp"
#include "DetectorEvent.h"

/// @brief algorithm tagging anomalous BXs by analysing the accelerator state
class AcceleratorStateFilter {
 public:
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
    /// Number of the accelerator state PVs read,
    /// in the reading order up to the first missing one
    std::size_t nPvs = 0;
    /// All accelerator state PVs presence flag
    bool hasPvs = false;
    /// Accelerator state
    AcceleratorState state;
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

  /// @brief Get the accelerator state values in the PV reading order
  ///
  /// @param state accelerator state
  ///
  /// @return accelerator state values
  static std::vector<double> getStateValues(const AcceleratorState &state);

 private:
  /// Configuration
  Config m_cfg;

  /// @brief Get the accelerator state means and std devs
  ///
  /// @param bxs accelerator data of the run BXs
  ///
  /// @return means and std devs of the accelerator diagnostics
  std::pair<AcceleratorState, AcceleratorState> getAcceleratorStateMeanStd(
      const std::vector<BxData> &bxs) const;
};
