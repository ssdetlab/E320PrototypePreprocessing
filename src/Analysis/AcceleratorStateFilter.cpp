#include "Analysis/AcceleratorStateFilter.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Analysis/PVList.hpp"
#include "TH1.h"
#include "detail/EpicsUtils.hpp"

AcceleratorStateFilter::BxData AcceleratorStateFilter::getBxData(
    const detector_event *detEvent, std::size_t entry) const {
  using namespace EPICS;

  BxData bx;
  bx.entry = entry;
  bx.trgN = detEvent->trg_n;

  // PID and shutter state
  try {
    std::size_t epicsPID = getPv(detEvent, epics_pid_key);
    bx.hasPid = true;
    bx.isShuttered =
        static_cast<int>(std::round((epicsPID - 10) / 36)) % 20 == 0;
  } catch (std::out_of_range &err) {
    bx.hasPid = false;
  }

  // Accelerator state PVs, read in order up to the first missing one
  try {
    // PMT_S20_3060
    bx.state.pmt_s20_3060 = getPv(detEvent, pmt_s20_3060_key);
    bx.nPvs++;
    // PMT_S20_3070
    bx.state.pmt_s20_3070 = getPv(detEvent, pmt_s20_3070_key);
    bx.nPvs++;
    // PMT_S20_3179
    bx.state.pmt_s20_3179 = getPv(detEvent, pmt_s20_3179_key);
    bx.nPvs++;
    // PMT_S20_3350
    bx.state.pmt_s20_3350 = getPv(detEvent, pmt_s20_3350_key);
    bx.nPvs++;
    // PMT_S20_3360
    bx.state.pmt_s20_3360 = getPv(detEvent, pmt_s20_3360_key);
    bx.nPvs++;
    // TORO:LI20:1988:TMIT_PC
    bx.state.toro_li20_1988_tmit_pc =
        getPv(detEvent, toro_li20_1988_tmit_pc_key);
    bx.nPvs++;
    // TORO:LI20:2040:TMIT_PC
    bx.state.toro_li20_2040_tmit_pc =
        getPv(detEvent, toro_li20_2040_tmit_pc_key);
    bx.nPvs++;
    // TORO:LI20:2452:TMIT_PC
    bx.state.toro_li20_2452_tmit_pc =
        getPv(detEvent, toro_li20_2452_tmit_pc_key);
    bx.nPvs++;
    // TORO:LI20:3163:TMIT_PC
    bx.state.toro_li20_3163_tmit_pc =
        getPv(detEvent, toro_li20_3163_tmit_pc_key);
    bx.nPvs++;
    // TORO:LI20:3255:TMIT_PC
    bx.state.toro_li20_3255_tmit_pc =
        getPv(detEvent, toro_li20_3255_tmit_pc_key);
    bx.nPvs++;
    bx.hasPvs = true;
  } catch (std::out_of_range &err) {
    std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
  }

  return bx;
}

std::vector<double> AcceleratorStateFilter::getStateValues(
    const AcceleratorState &state) {
  return {state.pmt_s20_3060,           state.pmt_s20_3070,
          state.pmt_s20_3179,           state.pmt_s20_3350,
          state.pmt_s20_3360,           state.toro_li20_1988_tmit_pc,
          state.toro_li20_2040_tmit_pc, state.toro_li20_2452_tmit_pc,
          state.toro_li20_3163_tmit_pc, state.toro_li20_3255_tmit_pc};
}

std::pair<AcceleratorState, AcceleratorState>
AcceleratorStateFilter::getAcceleratorStateMeanStd(
    const std::vector<BxData> &bxs) const {
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
  // TORO:LI20:3163:TMIT_PC
  TH1D toro_li20_3163_tmit_pc_hist(
      "toro_li20_3163_tmit_pc", "toro_li20_3163_tmit_pc", 10000, -5000, 5000);
  // TORO:LI20:3255:TMIT_PC
  TH1D toro_li20_3255_tmit_pc_hist(
      "toro_li20_3255_tmit_pc", "toro_li20_3255_tmit_pc", 10000, -5000, 5000);
  std::vector<TH1D *> hists = {
      &pmt_s20_3060_hist,           &pmt_s20_3070_hist,
      &pmt_s20_3179_hist,           &pmt_s20_3350_hist,
      &pmt_s20_3360_hist,           &toro_li20_1988_tmit_pc_hist,
      &toro_li20_2040_tmit_pc_hist, &toro_li20_2452_tmit_pc_hist,
      &toro_li20_3163_tmit_pc_hist, &toro_li20_3255_tmit_pc_hist};

  for (const auto &bx : bxs) {
    if (!bx.hasPid || bx.isShuttered) {
      continue;
    }

    // Fill the PVs read before the first missing one
    auto values = getStateValues(bx.state);
    for (std::size_t j = 0; j < bx.nPvs; j++) {
      hists.at(j)->Fill(values.at(j));
    }
  }

  AcceleratorState accStateMin;
  AcceleratorState accStateMax;
  double quantilePos[2] = {m_cfg.noiseQuantileLow, m_cfg.noiseQuantileHigh};

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

  double toro_li20_3163_tmit_pc_qs[2];
  toro_li20_3163_tmit_pc_hist.GetQuantiles(2, toro_li20_3163_tmit_pc_qs,
                                           quantilePos);
  accStateMin.toro_li20_3163_tmit_pc = toro_li20_3163_tmit_pc_qs[0];
  accStateMax.toro_li20_3163_tmit_pc = toro_li20_3163_tmit_pc_qs[1];

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
  // TORO:LI20:3163:TMIT_PC
  toro_li20_3163_tmit_pc_hist.Reset();
  // TORO:LI20:3255:TMIT_PC
  toro_li20_3255_tmit_pc_hist.Reset();

  AcceleratorState accStateMean;
  AcceleratorState accStateStdDev;
  for (const auto &bx : bxs) {
    if (!bx.hasPid || bx.isShuttered) {
      continue;
    }
    if (!bx.hasPvs) {
      continue;
    }

    if (bx.state.pmt_s20_3060 < accStateMin.pmt_s20_3060 ||
        bx.state.pmt_s20_3060 > accStateMax.pmt_s20_3060) {
      continue;
    }
    if (bx.state.pmt_s20_3070 < accStateMin.pmt_s20_3070 ||
        bx.state.pmt_s20_3070 > accStateMax.pmt_s20_3070) {
      continue;
    }
    if (bx.state.pmt_s20_3179 < accStateMin.pmt_s20_3179 ||
        bx.state.pmt_s20_3179 > accStateMax.pmt_s20_3179) {
      continue;
    }
    if (bx.state.pmt_s20_3350 < accStateMin.pmt_s20_3350 ||
        bx.state.pmt_s20_3350 > accStateMax.pmt_s20_3350) {
      continue;
    }
    if (bx.state.pmt_s20_3360 < accStateMin.pmt_s20_3360 ||
        bx.state.pmt_s20_3360 > accStateMax.pmt_s20_3360) {
      continue;
    }
    if (bx.state.toro_li20_1988_tmit_pc < accStateMin.toro_li20_1988_tmit_pc ||
        bx.state.toro_li20_1988_tmit_pc > accStateMax.toro_li20_1988_tmit_pc) {
      continue;
    }
    if (bx.state.toro_li20_2040_tmit_pc < accStateMin.toro_li20_2040_tmit_pc ||
        bx.state.toro_li20_2040_tmit_pc > accStateMax.toro_li20_2040_tmit_pc) {
      continue;
    }
    if (bx.state.toro_li20_2452_tmit_pc < accStateMin.toro_li20_2452_tmit_pc ||
        bx.state.toro_li20_2452_tmit_pc > accStateMax.toro_li20_2452_tmit_pc) {
      continue;
    }
    if (bx.state.toro_li20_3163_tmit_pc < accStateMin.toro_li20_3163_tmit_pc ||
        bx.state.toro_li20_3163_tmit_pc > accStateMax.toro_li20_3163_tmit_pc) {
      continue;
    }
    if (bx.state.toro_li20_3255_tmit_pc < accStateMin.toro_li20_3255_tmit_pc ||
        bx.state.toro_li20_3255_tmit_pc > accStateMax.toro_li20_3255_tmit_pc) {
      continue;
    }

    // PMT_S20_3060
    pmt_s20_3060_hist.Fill(bx.state.pmt_s20_3060);
    // PMT_S20_3070
    pmt_s20_3070_hist.Fill(bx.state.pmt_s20_3070);
    // PMT_S20_3179
    pmt_s20_3179_hist.Fill(bx.state.pmt_s20_3179);
    // PMT_S20_3350
    pmt_s20_3350_hist.Fill(bx.state.pmt_s20_3350);
    // PMT_S20_3360
    pmt_s20_3360_hist.Fill(bx.state.pmt_s20_3360);
    // TORO:LI20:1988:TMIT_PC
    toro_li20_1988_tmit_pc_hist.Fill(bx.state.toro_li20_1988_tmit_pc);
    // TORO:LI20:2040:TMIT_PC
    toro_li20_2040_tmit_pc_hist.Fill(bx.state.toro_li20_2040_tmit_pc);
    // TORO:LI20:2452:TMIT_PC
    toro_li20_2452_tmit_pc_hist.Fill(bx.state.toro_li20_2452_tmit_pc);
    // TORO:LI20:3163:TMIT_PC
    toro_li20_3163_tmit_pc_hist.Fill(bx.state.toro_li20_3163_tmit_pc);
    // TORO:LI20:3255:TMIT_PC
    toro_li20_3255_tmit_pc_hist.Fill(bx.state.toro_li20_3255_tmit_pc);
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
  // TORO:LI20:3163:TMIT_PC
  accStateMean.toro_li20_3163_tmit_pc = toro_li20_3163_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_3163_tmit_pc =
      toro_li20_3163_tmit_pc_hist.GetStdDev();
  // TORO:LI20:3255:TMIT_PC
  accStateMean.toro_li20_3255_tmit_pc = toro_li20_3255_tmit_pc_hist.GetMean();
  accStateStdDev.toro_li20_3255_tmit_pc =
      toro_li20_3255_tmit_pc_hist.GetStdDev();

  return {accStateMean, accStateStdDev};
}

std::unordered_set<std::size_t> AcceleratorStateFilter::getBadEvents(
    const std::vector<BxData> &bxs) const {
  std::cout << "Getting bad events with AccStdDev = " << m_cfg.nStdDevs << "\n";
  // Get accelerator state parameters
  auto [accStateMean, accStateStdDev] = getAcceleratorStateMeanStd(bxs);

  std::unordered_set<std::size_t> badEvents;
  auto isBadCentral = [this](double mean, double std, double val) {
    return (std::abs(mean - val) / std > m_cfg.nStdDevs);
  };
  auto isBadUpper = [this](double mean, double std, double val) {
    return ((val - mean) / std > m_cfg.nStdDevs);
  };

  for (const auto &bx : bxs) {
    if (!bx.hasPid) {
      badEvents.insert(bx.entry);
      continue;
    }

    // Shuttered BXs are not classified
    if (bx.isShuttered) {
      continue;
    }

    if (!bx.hasPvs) {
      badEvents.insert(bx.entry);
      continue;
    }
    // PMT_S20_3060
    if (isBadUpper(accStateMean.pmt_s20_3060, accStateStdDev.pmt_s20_3060,
                   bx.state.pmt_s20_3060)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // PMT_S20_3070
    if (isBadUpper(accStateMean.pmt_s20_3070, accStateStdDev.pmt_s20_3070,
                   bx.state.pmt_s20_3070)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // PMT_S20_3179
    if (isBadUpper(accStateMean.pmt_s20_3179, accStateStdDev.pmt_s20_3179,
                   bx.state.pmt_s20_3179)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // PMT_S20_3350
    if (isBadUpper(accStateMean.pmt_s20_3350, accStateStdDev.pmt_s20_3350,
                   bx.state.pmt_s20_3350)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // PMT_S20_3360
    if (isBadUpper(accStateMean.pmt_s20_3360, accStateStdDev.pmt_s20_3360,
                   bx.state.pmt_s20_3360)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // TORO:LI20:1988:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_1988_tmit_pc,
                     accStateStdDev.toro_li20_1988_tmit_pc,
                     bx.state.toro_li20_1988_tmit_pc)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // TORO:LI20:2040:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_2040_tmit_pc,
                     accStateStdDev.toro_li20_2040_tmit_pc,
                     bx.state.toro_li20_2040_tmit_pc)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // TORO:LI20:2452:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_2452_tmit_pc,
                     accStateStdDev.toro_li20_2452_tmit_pc,
                     bx.state.toro_li20_2452_tmit_pc)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // TORO:LI20:3163:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_3163_tmit_pc,
                     accStateStdDev.toro_li20_3163_tmit_pc,
                     bx.state.toro_li20_3163_tmit_pc)) {
      badEvents.insert(bx.entry);
      continue;
    }
    // TORO:LI20:3255:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_3255_tmit_pc,
                     accStateStdDev.toro_li20_3255_tmit_pc,
                     bx.state.toro_li20_3255_tmit_pc)) {
      badEvents.insert(bx.entry);
      continue;
    }
  }
  std::cout << "Found bad events " << badEvents.size() << "\n";
  std::cout << "\n\n\n\n";
  return badEvents;
}
