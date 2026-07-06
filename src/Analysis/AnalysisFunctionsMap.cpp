#include <cmath>
#include <cstddef>
#include <deque>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>

#include "Analysis/AnalysisFuncitons.hpp"
#include "Analysis/PVList.hpp"
#include "DetectorEvent.h"

double convertEpicsDipole(double val, double theta, double bz) {
  const double c_light = 2.99792458E8;  // speed of light in m/s
  const double me = 0.5109989461e6;     // electron mass in eV

  return std::sqrt(val * val - me * me) * std::sin(theta) / (c_light * bz);
}

double convertEpicsXCorrector(double val, double bz) {
  return val / bz * 0.1;
}

double convertEpicsQuad(double val) {
  return val * 0.1;
}

double getPv(const detector_event *detEvent, const std::string &pvname) {
  try {
    return std::stod(detEvent->ev_epics_frame.pv_map.at(pvname).value);
  } catch (std::out_of_range &err) {
    std::cout << "MISSING PV " << pvname << "\n";
    throw err;
  }
}

std::pair<AcceleratorState, AcceleratorState> getAcceleratorStateMeanStd(
    TTree *tree, std::size_t skip) {
  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);

  // PMT_S20_3060
  TH1D pmt_s20_3060_hist("pmt_s20_3060", "pmt_s20_3060", 10000, 0, 50000);
  // PMT_S20_3070
  TH1D pmt_s20_3070_hist("pmt_s20_3070", "pmt_s20_3070", 10000, 0, 50000);
  // PMT_S20_3179
  TH1D pmt_s20_3179_hist("pmt_s20_3179", "pmt_s20_3179", 10000, 0, 50000);
  // PMT_S20_3350
  TH1D pmt_s20_3350_hist("pmt_s20_3350", "pmt_s20_3350", 10000, 0, 50000);
  // PMT_S20_3360
  TH1D pmt_s20_3360_hist("pmt_s20_3360", "pmt_s20_3360", 10000, 0, 50000);
  // TORO:LI20:1988:TMIT_PC
  TH1D toro_li20_1988_tmit_pc_hist(
      "toro_li20_1988_tmit_pc", "toro_li20_1988_tmit_pc", 10000, -5000, 5000);
  // TORO:LI20:2040:TMIT_PC
  TH1D toro_li20_2040_tmit_pc_hist(
      "toro_li20_2040_tmit_pc", "toro_li20_2040_tmit_pc", 10000, -5000, 5000);
  // TORO:LI20:2452:TMIT_PC
  TH1D toro_li20_2452_tmit_pc_hist(
      "toro_li20_2452_tmit_pc", "toro_li20_2452_tmit_pc", 10000, -5000, 5000);
  // // TORO:LI20:3163:TMIT_PC
  // TH1D toro_li20_3163_tmit_pc_hist(
  //     "toro_li20_3163_tmit_pc", "toro_li20_3163_tmit_pc", 10000, -5000,
  //     5000);
  // TORO:LI20:3163:TMIT_PC
  TH1D toro_li20_3255_tmit_pc_hist(
      "toro_li20_3255_tmit_pc", "toro_li20_3255_tmit_pc", 10000, -5000, 5000);

  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    std::size_t epicsPID = 0;
    try {
      epicsPID = getPv(detEvent, epics_pid_key);
    } catch (std::out_of_range &err) {
      continue;
    }
    bool isEven = static_cast<int>(std::round((epicsPID - 10) / 36)) % 2 == 0;
    bool isShuttered =
        isEven && static_cast<int>(std::round((epicsPID - 10) / 36)) % 20 == 0;
    if (isShuttered) {
      continue;
    }

    try {
      // PMT_S20_3060
      double pmt_s20_3060_val = getPv(detEvent, pmt_s20_3060_key);
      pmt_s20_3060_hist.Fill(pmt_s20_3060_val);
      // PMT_S20_3070
      double pmt_s20_3070_val = getPv(detEvent, pmt_s20_3070_key);
      pmt_s20_3070_hist.Fill(pmt_s20_3070_val);
      // PMT_S20_3179
      double pmt_s20_3179_val = getPv(detEvent, pmt_s20_3179_key);
      pmt_s20_3179_hist.Fill(pmt_s20_3179_val);
      // PMT_S20_3350
      double pmt_s20_3350_val = getPv(detEvent, pmt_s20_3350_key);
      pmt_s20_3350_hist.Fill(pmt_s20_3350_val);
      // PMT_S20_3360
      double pmt_s20_3360_val = getPv(detEvent, pmt_s20_3360_key);
      pmt_s20_3360_hist.Fill(pmt_s20_3360_val);

      // TORO:LI20:1988:TMIT_PC
      double toro_li20_1988_tmit_pc_val =
          getPv(detEvent, toro_li20_1988_tmit_pc_key);
      toro_li20_1988_tmit_pc_hist.Fill(toro_li20_1988_tmit_pc_val);
      // TORO:LI20:2040:TMIT_PC
      double toro_li20_2040_tmit_pc_val =
          getPv(detEvent, toro_li20_2040_tmit_pc_key);
      toro_li20_2040_tmit_pc_hist.Fill(toro_li20_2040_tmit_pc_val);
      // TORO:LI20:2452:TMIT_PC
      double toro_li20_2452_tmit_pc_val =
          getPv(detEvent, toro_li20_2452_tmit_pc_key);
      toro_li20_2452_tmit_pc_hist.Fill(toro_li20_2452_tmit_pc_val);
      // // TORO:LI20:3163:TMIT_PC
      // double toro_li20_3163_tmit_pc_val =
      //     getPv(detEvent, toro_li20_3163_tmit_pc_key);
      // toro_li20_3163_tmit_pc_hist.Fill(toro_li20_3163_tmit_pc_val);
      // TORO:LI20:3255:TMIT_PC
      double toro_li20_3255_tmit_pc_val =
          getPv(detEvent, toro_li20_3255_tmit_pc_key);
      toro_li20_3255_tmit_pc_hist.Fill(toro_li20_3255_tmit_pc_val);
    } catch (std::out_of_range &err) {
      std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
      continue;
    }
  }

  AcceleratorState accStateMin;
  AcceleratorState accStateMax;
  double quantilePos[2] = {0.01, 0.99};

  double pmt_s20_3060_qs[2];
  pmt_s20_3060_hist.GetQuantiles(2, pmt_s20_3060_qs, quantilePos);
  accStateMin.pmt_s20_3060 = pmt_s20_3060_qs[0];
  accStateMax.pmt_s20_3060 = pmt_s20_3060_qs[1];

  double pmt_s20_3070_qs[2];
  pmt_s20_3070_hist.GetQuantiles(2, pmt_s20_3070_qs, quantilePos);
  accStateMin.pmt_s20_3070 = pmt_s20_3070_qs[0];
  accStateMax.pmt_s20_3070 = pmt_s20_3070_qs[1];

  double pmt_s20_3179_qs[2];
  pmt_s20_3179_hist.GetQuantiles(2, pmt_s20_3179_qs, quantilePos);
  accStateMin.pmt_s20_3179 = pmt_s20_3179_qs[0];
  accStateMax.pmt_s20_3179 = pmt_s20_3179_qs[1];

  double pmt_s20_3350_qs[2];
  pmt_s20_3350_hist.GetQuantiles(2, pmt_s20_3350_qs, quantilePos);
  accStateMin.pmt_s20_3350 = pmt_s20_3350_qs[0];
  accStateMax.pmt_s20_3350 = pmt_s20_3350_qs[1];

  double pmt_s20_3360_qs[2];
  pmt_s20_3360_hist.GetQuantiles(2, pmt_s20_3360_qs, quantilePos);
  accStateMin.pmt_s20_3360 = pmt_s20_3360_qs[0];
  accStateMax.pmt_s20_3360 = pmt_s20_3360_qs[1];

  double toro_li20_1988_tmit_pc_qs[2];
  toro_li20_1988_tmit_pc_hist.GetQuantiles(2, toro_li20_1988_tmit_pc_qs,
                                           quantilePos);
  accStateMin.toro_li20_1988_tmit_pc = toro_li20_1988_tmit_pc_qs[0];
  accStateMax.toro_li20_1988_tmit_pc = toro_li20_1988_tmit_pc_qs[1];

  double toro_li20_2040_tmit_pc_qs[2];
  toro_li20_2040_tmit_pc_hist.GetQuantiles(2, toro_li20_2040_tmit_pc_qs,
                                           quantilePos);
  accStateMin.toro_li20_2040_tmit_pc = toro_li20_2040_tmit_pc_qs[0];
  accStateMax.toro_li20_2040_tmit_pc = toro_li20_2040_tmit_pc_qs[1];

  double toro_li20_2452_tmit_pc_qs[2];
  toro_li20_2452_tmit_pc_hist.GetQuantiles(2, toro_li20_2452_tmit_pc_qs,
                                           quantilePos);
  accStateMin.toro_li20_2452_tmit_pc = toro_li20_2452_tmit_pc_qs[0];
  accStateMax.toro_li20_2452_tmit_pc = toro_li20_2452_tmit_pc_qs[1];

  // double toro_li20_3163_tmit_pc_qs[2];
  // toro_li20_3163_tmit_pc_hist.GetQuantiles(2, toro_li20_3163_tmit_pc_qs,
  //                                          quantilePos);
  // accStateMin.toro_li20_3163_tmit_pc = toro_li20_3163_tmit_pc_qs[0];
  // accStateMax.toro_li20_3163_tmit_pc = toro_li20_3163_tmit_pc_qs[1];

  double toro_li20_3255_tmit_pc_qs[2];
  toro_li20_3255_tmit_pc_hist.GetQuantiles(2, toro_li20_3255_tmit_pc_qs,
                                           quantilePos);
  accStateMin.toro_li20_3255_tmit_pc = toro_li20_3255_tmit_pc_qs[0];
  accStateMax.toro_li20_3255_tmit_pc = toro_li20_3255_tmit_pc_qs[1];

  // PMT_S20_3060
  pmt_s20_3060_hist.Reset();
  // PMT_S20_3070
  pmt_s20_3070_hist.Reset();
  // PMT_S20_3179
  pmt_s20_3179_hist.Reset();
  // PMT_S20_3350
  pmt_s20_3350_hist.Reset();
  // PMT_S20_3360
  pmt_s20_3360_hist.Reset();
  // TORO:LI20:1988:TMIT_PC
  toro_li20_1988_tmit_pc_hist.Reset();
  // TORO:LI20:2040:TMIT_PC
  toro_li20_2040_tmit_pc_hist.Reset();
  // TORO:LI20:2452:TMIT_PC
  toro_li20_2452_tmit_pc_hist.Reset();
  // // TORO:LI20:3163:TMIT_PC
  // toro_li20_3163_tmit_pc_hist.Reset();
  // TORO:LI20:3255:TMIT_PC
  toro_li20_3255_tmit_pc_hist.Reset();

  AcceleratorState accStateMean;
  AcceleratorState accStateStdDev;
  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    double pmt_s20_3060_val = 0;
    double pmt_s20_3070_val = 0;
    double pmt_s20_3179_val = 0;
    double pmt_s20_3350_val = 0;
    double pmt_s20_3360_val = 0;

    double toro_li20_1988_tmit_pc_val = 0;
    double toro_li20_2040_tmit_pc_val = 0;
    double toro_li20_2452_tmit_pc_val = 0;
    // double toro_li20_3163_tmit_pc_val = 0;
    double toro_li20_3255_tmit_pc_val = 0;
    try {
      pmt_s20_3060_val = getPv(detEvent, pmt_s20_3060_key);
      pmt_s20_3070_val = getPv(detEvent, pmt_s20_3070_key);
      pmt_s20_3179_val = getPv(detEvent, pmt_s20_3179_key);
      pmt_s20_3350_val = getPv(detEvent, pmt_s20_3350_key);
      pmt_s20_3360_val = getPv(detEvent, pmt_s20_3360_key);

      toro_li20_1988_tmit_pc_val = getPv(detEvent, toro_li20_1988_tmit_pc_key);
      toro_li20_2040_tmit_pc_val = getPv(detEvent, toro_li20_2040_tmit_pc_key);
      toro_li20_2452_tmit_pc_val = getPv(detEvent, toro_li20_2452_tmit_pc_key);
      // toro_li20_3163_tmit_pc_val =
      //     getPv(detEvent, toro_li20_3163_tmit_pc_key);
      toro_li20_3255_tmit_pc_val = getPv(detEvent, toro_li20_3255_tmit_pc_key);
    } catch (std::out_of_range &err) {
      std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
      continue;
    }

    if (pmt_s20_3060_val < accStateMin.pmt_s20_3060 ||
        pmt_s20_3060_val > accStateMax.pmt_s20_3060) {
      continue;
    }
    if (pmt_s20_3070_val < accStateMin.pmt_s20_3070 ||
        pmt_s20_3070_val > accStateMax.pmt_s20_3070) {
      continue;
    }
    if (pmt_s20_3179_val < accStateMin.pmt_s20_3179 ||
        pmt_s20_3179_val > accStateMax.pmt_s20_3179) {
      continue;
    }
    if (pmt_s20_3350_val < accStateMin.pmt_s20_3350 ||
        pmt_s20_3350_val > accStateMax.pmt_s20_3350) {
      continue;
    }
    if (pmt_s20_3360_val < accStateMin.pmt_s20_3360 ||
        pmt_s20_3360_val > accStateMax.pmt_s20_3360) {
      continue;
    }

    if (toro_li20_1988_tmit_pc_val < accStateMin.toro_li20_1988_tmit_pc ||
        toro_li20_1988_tmit_pc_val > accStateMax.toro_li20_1988_tmit_pc) {
      continue;
    }
    if (toro_li20_2040_tmit_pc_val < accStateMin.toro_li20_2040_tmit_pc ||
        toro_li20_2040_tmit_pc_val > accStateMax.toro_li20_2040_tmit_pc) {
      continue;
    }
    if (toro_li20_2452_tmit_pc_val < accStateMin.toro_li20_2452_tmit_pc ||
        toro_li20_2452_tmit_pc_val > accStateMax.toro_li20_2452_tmit_pc) {
      continue;
    }
    // if (toro_li20_3163_tmit_pc_val < accStateMin.toro_li20_3163_tmit_pc ||
    //     toro_li20_3163_tmit_pc_val > accStateMax.toro_li20_3163_tmit_pc) {
    //   continue;
    // }
    if (toro_li20_3255_tmit_pc_val < accStateMin.toro_li20_3255_tmit_pc ||
        toro_li20_3255_tmit_pc_val > accStateMax.toro_li20_3255_tmit_pc) {
      continue;
    }

    // PMT_S20_3060
    pmt_s20_3060_hist.Fill(pmt_s20_3060_val);
    // PMT_S20_3070
    pmt_s20_3070_hist.Fill(pmt_s20_3070_val);
    // PMT_S20_3179
    pmt_s20_3179_hist.Fill(pmt_s20_3179_val);
    // PMT_S20_3350
    pmt_s20_3350_hist.Fill(pmt_s20_3350_val);
    // PMT_S20_3360
    pmt_s20_3360_hist.Fill(pmt_s20_3360_val);
    // TORO:LI20:1988:TMIT_PC
    toro_li20_1988_tmit_pc_hist.Fill(toro_li20_1988_tmit_pc_val);
    // TORO:LI20:2040:TMIT_PC
    toro_li20_2040_tmit_pc_hist.Fill(toro_li20_2040_tmit_pc_val);
    // TORO:LI20:2452:TMIT_PC
    toro_li20_2452_tmit_pc_hist.Fill(toro_li20_2452_tmit_pc_val);
    // // TORO:LI20:3163:TMIT_PC
    // toro_li20_3163_tmit_pc_hist.Fill(toro_li20_3163_tmit_pc_val);
    // TORO:LI20:3255:TMIT_PC
    toro_li20_3255_tmit_pc_hist.Fill(toro_li20_3255_tmit_pc_val);
  }

  // PMT_S20_3060
  accStateMean.pmt_s20_3060 = pmt_s20_3060_hist.GetMean();
  accStateStdDev.pmt_s20_3060 = pmt_s20_3060_hist.GetStdDev();
  // PMT_S20_3070
  accStateMean.pmt_s20_3070 = pmt_s20_3070_hist.GetMean();
  accStateStdDev.pmt_s20_3070 = pmt_s20_3070_hist.GetStdDev();
  // PMT_S20_3179
  accStateMean.pmt_s20_3179 = pmt_s20_3179_hist.GetMean();
  accStateStdDev.pmt_s20_3179 = pmt_s20_3179_hist.GetStdDev();
  // PMT_S20_3350
  accStateMean.pmt_s20_3350 = pmt_s20_3350_hist.GetMean();
  accStateStdDev.pmt_s20_3350 = pmt_s20_3350_hist.GetStdDev();
  // PMT_S20_3360
  accStateMean.pmt_s20_3360 = pmt_s20_3360_hist.GetMean();
  accStateStdDev.pmt_s20_3360 = pmt_s20_3360_hist.GetStdDev();

  // TORO:LI20:1988:TMIT_PC
  accStateMean.toro_li20_1988_tmit_pc = toro_li20_1988_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_1988_tmit_pc =
      toro_li20_1988_tmit_pc_hist.GetStdDev();
  // TORO:LI20:2040:TMIT_PC
  accStateMean.toro_li20_2040_tmit_pc = toro_li20_2040_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_2040_tmit_pc =
      toro_li20_2040_tmit_pc_hist.GetStdDev();
  // TORO:LI20:2452:TMIT_PC
  accStateMean.toro_li20_2452_tmit_pc = toro_li20_2452_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_2452_tmit_pc =
      toro_li20_2452_tmit_pc_hist.GetStdDev();
  // // TORO:LI20:3163:TMIT_PC
  // accStateMean.toro_li20_3163_tmit_pc =
  // toro_li20_3163_tmit_pc_hist.GetMean();
  // accStateStdDev.toro_li20_3163_tmit_pc =
  //     toro_li20_3163_tmit_pc_hist.GetStdDev();
  // TORO:LI20:3255:TMIT_PC
  accStateMean.toro_li20_3255_tmit_pc = toro_li20_3255_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_3255_tmit_pc =
      toro_li20_3255_tmit_pc_hist.GetStdDev();

  return {accStateMean, accStateStdDev};
}

std::unordered_set<std::size_t> getBadEvents(TTree *tree, std::size_t skip,
                                             int nStdDevs) {
  // Get accelerator state parameters
  auto [accStateMean, accStateStdDev] = getAcceleratorStateMeanStd(tree, skip);

  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);

  std::unordered_set<std::size_t> badEvents;
  auto isBadTwoSided = [&nStdDevs](double mean, double std, double val) {
    return (std::abs(mean - val) / std > nStdDevs);
  };
  auto isBadUpper = [&nStdDevs](double mean, double std, double val) {
    return ((val - mean) / std > nStdDevs);
  };

  std::size_t nMissingPV = 0;
  std::size_t nShuttered = 0;
  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    std::size_t epicsPID = 0;
    try {
      epicsPID = getPv(detEvent, epics_pid_key);
    } catch (std::out_of_range &err) {
      badEvents.insert(i);
      nMissingPV++;
      continue;
    }
    bool isEven = static_cast<int>(std::round((epicsPID - 10) / 36)) % 2 == 0;
    bool isShuttered =
        isEven && static_cast<int>(std::round((epicsPID - 10) / 36)) % 20 == 0;
    if (isShuttered) {
      badEvents.insert(i);
      nShuttered++;
      continue;
    }

    double pmt_s20_3060_val = 0;
    double pmt_s20_3070_val = 0;
    double pmt_s20_3179_val = 0;
    double pmt_s20_3350_val = 0;
    double pmt_s20_3360_val = 0;

    double toro_li20_1988_tmit_pc_val = 0;
    double toro_li20_2040_tmit_pc_val = 0;
    double toro_li20_2452_tmit_pc_val = 0;
    // double toro_li20_3163_tmit_pc_val =
    //     0;
    double toro_li20_3255_tmit_pc_val = 0;
    try {
      pmt_s20_3060_val = getPv(detEvent, pmt_s20_3060_key);
      pmt_s20_3070_val = getPv(detEvent, pmt_s20_3070_key);
      pmt_s20_3179_val = getPv(detEvent, pmt_s20_3179_key);
      pmt_s20_3350_val = getPv(detEvent, pmt_s20_3350_key);
      pmt_s20_3360_val = getPv(detEvent, pmt_s20_3360_key);

      toro_li20_1988_tmit_pc_val = getPv(detEvent, toro_li20_1988_tmit_pc_key);
      toro_li20_2040_tmit_pc_val = getPv(detEvent, toro_li20_2040_tmit_pc_key);
      toro_li20_2452_tmit_pc_val = getPv(detEvent, toro_li20_2452_tmit_pc_key);
      // toro_li20_3163_tmit_pc_val =
      //     getPv(detEvent, toro_li20_3163_tmit_pc_key);
      toro_li20_3255_tmit_pc_val = getPv(detEvent, toro_li20_3255_tmit_pc_key);
    } catch (std::out_of_range &err) {
      std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
      badEvents.insert(i);
      nMissingPV++;
      continue;
    }

    // PMT_S20_3060
    if (isBadUpper(accStateMean.pmt_s20_3060, accStateStdDev.pmt_s20_3060,
                   pmt_s20_3060_val)) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3070
    if (isBadUpper(accStateMean.pmt_s20_3070, accStateStdDev.pmt_s20_3070,
              pmt_s20_3070_val)) {
      badEvents.insert(i);
      continue;
    }
    if (pmt_s20_3070_val > 580) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3179
    if (isBadUpper(accStateMean.pmt_s20_3179, accStateStdDev.pmt_s20_3179,
              pmt_s20_3179_val)) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3350
    if (isBadUpper(accStateMean.pmt_s20_3350, accStateStdDev.pmt_s20_3350,
              pmt_s20_3350_val)) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3360
    if (isBadUpper(accStateMean.pmt_s20_3360, accStateStdDev.pmt_s20_3360,
              pmt_s20_3360_val)) {
      badEvents.insert(i);
      continue;
    }

    // TORO:LI20:1988:TMIT_PC
    if (isBadUpper(accStateMean.toro_li20_1988_tmit_pc,
              accStateStdDev.toro_li20_1988_tmit_pc,
              toro_li20_1988_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2040:TMIT_PC
    if (isBadUpper(accStateMean.toro_li20_2040_tmit_pc,
              accStateStdDev.toro_li20_2040_tmit_pc,
              toro_li20_2040_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2452:TMIT_PC
    if (isBadUpper(accStateMean.toro_li20_2452_tmit_pc,
              accStateStdDev.toro_li20_2452_tmit_pc,
              toro_li20_2452_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // // TORO:LI20:3163:TMIT_PC
    // if (isBadUpper(accStateMean.toro_li20_3163_tmit_pc,
    //           accStateStdDev.toro_li20_3163_tmit_pc,
    //           toro_li20_3163_tmit_pc_val)) {
    //   badEvents.insert(i);
    //   continue;
    // }
    // TORO:LI20:3255:TMIT_PC
    if (isBadUpper(accStateMean.toro_li20_3255_tmit_pc,
              accStateStdDev.toro_li20_3255_tmit_pc,
              toro_li20_3255_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
  }
  std::cout << "MISSING PVS " << nMissingPV << "\n";
  std::cout << "SHUTTERED " << nShuttered << "\n";
  return badEvents;
}

void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                   HitSet &pixels) {
  std::deque<Hit> queue{pivot};
  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop_front();
    auto currentPtr = pixels.find(current);

    if (currentPtr != pixels.end()) {
      // add this pixel to the cluster
      clusterPixels.push_back(current);
      // remove pixel from live pixels list
      pixels.erase(currentPtr);

      std::vector<Hit> neighbours{{current.first + 1, current.second},
                                  {current.first - 1, current.second},
                                  {current.first, current.second + 1},
                                  {current.first, current.second - 1}};
      for (const auto &n : neighbours) {
        auto nptr = pixels.find(n);
        if (nptr != pixels.end()) {
          queue.push_back(*nptr);
        }
      }
    }
  }
}

std::vector<Cluster> getClusters(HitSet &hits, std::size_t maxCluterSize,
                                 const std::shared_ptr<ROI> &roi) {
  std::vector<Cluster> clusters;

  while (!hits.empty()) {
    auto pixel = *hits.begin();
    std::vector<Hit> clusterPixels;
    bfsClustering(clusterPixels, pixel, hits);
    std::size_t clSize = clusterPixels.size();
    if (clSize > maxCluterSize) {
      continue;
    }

    std::size_t xMin = 1024;
    std::size_t xMax = 0;
    std::size_t yMin = 512;
    std::size_t yMax = 0;
    double xCenter = 0;
    double yCenter = 0;
    for (auto pix : clusterPixels) {
      xCenter += pix.first;
      yCenter += pix.second;

      if (pix.first <= xMin) {
        xMin = pix.first;
      }
      if (pix.first >= xMax) {
        xMax = pix.first;
      }
      if (pix.second <= yMin) {
        yMin = pix.second;
      }
      if (pix.second >= yMax) {
        yMax = pix.second;
      }
    }
    std::size_t dx = xMax - xMin + 1;
    std::size_t dy = yMax - yMin + 1;
    xCenter /= clusterPixels.size();
    yCenter /= clusterPixels.size();

    if (!roi->isInside(xCenter, yCenter)) {
      continue;
    }

    clusters.push_back(Cluster{xCenter, yCenter, dx, dy, clSize});
  }

  return clusters;
}

std::vector<Cluster> getClusters(HitSet &hits, std::size_t maxCluterSize) {
  std::vector<Cluster> clusters;

  while (!hits.empty()) {
    auto pixel = *hits.begin();
    std::vector<Hit> clusterPixels;
    bfsClustering(clusterPixels, pixel, hits);
    std::size_t clSize = clusterPixels.size();
    if (clSize > maxCluterSize) {
      continue;
    }

    std::size_t xMin = 1024;
    std::size_t xMax = 0;
    std::size_t yMin = 512;
    std::size_t yMax = 0;
    double xCenter = 0;
    double yCenter = 0;
    for (auto pix : clusterPixels) {
      xCenter += pix.first;
      yCenter += pix.second;

      if (pix.first <= xMin) {
        xMin = pix.first;
      }
      if (pix.first >= xMax) {
        xMax = pix.first;
      }
      if (pix.second <= yMin) {
        yMin = pix.second;
      }
      if (pix.second >= yMax) {
        yMax = pix.second;
      }
    }
    std::size_t dx = xMax - xMin + 1;
    std::size_t dy = yMax - yMin + 1;
    xCenter /= clusterPixels.size();
    yCenter /= clusterPixels.size();

    clusters.push_back(Cluster{xCenter, yCenter, dx, dy, clSize});
  }

  return clusters;
}
