#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>

#include <Eigen/Eigen>
#include <TCanvas.h>

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

class Stats {
 public:
  Stats() = delete;

  explicit Stats(const std::set<int> &ids);

  void init(const std::string &prefix);

  void reset();

  void print();

  /// Chip ids
  std::set<int> chipIds;
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
  /// Raw hit map ROI hists
  std::map<int, TH2D> chipHitMapROIRawHists;
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
