#include "DetectorEvent.h"
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <cmath>
#include <cstddef>
#include <deque>
#include <iostream>
#include <unordered_set>

#include "Analysis/AnalysisFuncitons.hpp"
#include "detail/Painters.hpp"

std::unordered_set<std::size_t>
getBadEvents(TTree *tree, const AcceleratorState &accStateMean,
             const AcceleratorState &accStateStdDev, int nStdDevs,
             std::size_t skip) {
  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);

  std::unordered_set<std::size_t> badEvents;
  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    // PMT_S20_3060
    if (std::abs(accStateMean.pmt_s20_3060 -
                 detEvent->epics_frame.pmt_s20_3060) /
            accStateStdDev.pmt_s20_3060 >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3070
    if (std::abs(accStateMean.pmt_s20_3070 -
                 detEvent->epics_frame.pmt_s20_3070) /
            accStateStdDev.pmt_s20_3070 >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3179
    if (std::abs(accStateMean.pmt_s20_3179 -
                 detEvent->epics_frame.pmt_s20_3179) /
            accStateStdDev.pmt_s20_3179 >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3350
    if (std::abs(accStateMean.pmt_s20_3350 -
                 detEvent->epics_frame.pmt_s20_3350) /
            accStateStdDev.pmt_s20_3350 >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // PMT_S20_3360
    if (std::abs(accStateMean.pmt_s20_3360 -
                 detEvent->epics_frame.pmt_s20_3360) /
            accStateStdDev.pmt_s20_3360 >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:1988:TMIT_PC
    if (std::abs(accStateMean.toro_li20_1988_tmit_pc -
                 detEvent->epics_frame.toro_li20_1988_tmit_pc) /
            accStateStdDev.toro_li20_1988_tmit_pc >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2040:TMIT_PC
    if (std::abs(accStateMean.toro_li20_2040_tmit_pc -
                 detEvent->epics_frame.toro_li20_2040_tmit_pc) /
            accStateStdDev.toro_li20_2040_tmit_pc >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2452:TMIT_PC
    if (std::abs(accStateMean.toro_li20_2452_tmit_pc -
                 detEvent->epics_frame.toro_li20_2452_tmit_pc) /
            accStateStdDev.toro_li20_2452_tmit_pc >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:3163:TMIT_PC
    if (std::abs(accStateMean.toro_li20_3163_tmit_pc -
                 detEvent->epics_frame.toro_li20_3163_tmit_pc) /
            accStateStdDev.toro_li20_3163_tmit_pc >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:3255:TMIT_PC
    if (std::abs(accStateMean.toro_li20_3255_tmit_pc -
                 detEvent->epics_frame.toro_li20_3255_tmit_pc) /
            accStateStdDev.toro_li20_3255_tmit_pc >
        nStdDevs) {
      badEvents.insert(i);
      continue;
    }
    // RADM:LI20:1:CH01:MEAS
    if (std::abs(accStateMean.radm_li20_1_ch01_meas -
                 detEvent->epics_frame.radm_li20_1_ch01_meas) /
            accStateStdDev.radm_li20_1_ch01_meas >
        nStdDevs) {

      badEvents.insert(i);
      continue;
    }
  }
  return badEvents;
}

std::pair<AcceleratorState, AcceleratorState>
getAcceleratorStateMeanStd(TTree *tree, std::size_t skip) {
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
  // TORO:LI20:3163:TMIT_PC
  TH1D toro_li20_3163_tmit_pc_hist(
      "toro_li20_3163_tmit_pc", "toro_li20_3163_tmit_pc", 10000, -5000, 5000);
  // TORO:LI20:3255:TMIT_PC
  TH1D toro_li20_3255_tmit_pc_hist(
      "toro_li20_3255_tmit_pc", "toro_li20_3255_tmit_pc", 10000, -5000, 5000);
  // RADM:LI20:1:CH01:MEAS
  TH1D radm_li20_1_ch01_meas_hist("radm_li20_1_ch01_meas",
                                  "radm_li20_1_ch01_meas", 1000, -5000, 5000);

  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    // PMT_S20_3060
    pmt_s20_3060_hist.Fill(detEvent->epics_frame.pmt_s20_3060);
    // PMT_S20_3070
    pmt_s20_3070_hist.Fill(detEvent->epics_frame.pmt_s20_3070);
    // PMT_S20_3179
    pmt_s20_3179_hist.Fill(detEvent->epics_frame.pmt_s20_3179);
    // PMT_S20_3350
    pmt_s20_3350_hist.Fill(detEvent->epics_frame.pmt_s20_3350);
    // PMT_S20_3360
    pmt_s20_3360_hist.Fill(detEvent->epics_frame.pmt_s20_3360);
    // TORO:LI20:1988:TMIT_PC
    toro_li20_1988_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_1988_tmit_pc);
    // TORO:LI20:2040:TMIT_PC
    toro_li20_2040_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_2040_tmit_pc);
    // TORO:LI20:2452:TMIT_PC
    toro_li20_2452_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_2452_tmit_pc);
    // TORO:LI20:3163:TMIT_PC
    toro_li20_3163_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_3163_tmit_pc);
    // TORO:LI20:3255:TMIT_PC
    toro_li20_3255_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_3255_tmit_pc);
    // RADM:LI20:1:CH01:MEAS
    radm_li20_1_ch01_meas_hist.Fill(
        detEvent->epics_frame.radm_li20_1_ch01_meas);
  }

  AcceleratorState accStateMin;
  AcceleratorState accStateMax;
  double quantilePos[2] = {0.1, 0.9};

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

  double radm_li20_1_ch01_meas_qs[2];
  radm_li20_1_ch01_meas_hist.GetQuantiles(2, radm_li20_1_ch01_meas_qs,
                                          quantilePos);
  accStateMin.radm_li20_1_ch01_meas = radm_li20_1_ch01_meas_qs[0];
  accStateMax.radm_li20_1_ch01_meas = radm_li20_1_ch01_meas_qs[1];

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
  // RADM:LI20:1:CH01:MEAS
  radm_li20_1_ch01_meas_hist.Reset();

  AcceleratorState accStateMean;
  AcceleratorState accStateStdDev;
  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (detEvent->epics_frame.pmt_s20_3060 < accStateMin.pmt_s20_3060 ||
        detEvent->epics_frame.pmt_s20_3060 > accStateMax.pmt_s20_3060) {
      continue;
    }
    if (detEvent->epics_frame.pmt_s20_3070 < accStateMin.pmt_s20_3070 ||
        detEvent->epics_frame.pmt_s20_3070 > accStateMax.pmt_s20_3070) {
      continue;
    }
    if (detEvent->epics_frame.pmt_s20_3179 < accStateMin.pmt_s20_3179 ||
        detEvent->epics_frame.pmt_s20_3179 > accStateMax.pmt_s20_3179) {
      continue;
    }
    if (detEvent->epics_frame.pmt_s20_3350 < accStateMin.pmt_s20_3350 ||
        detEvent->epics_frame.pmt_s20_3350 > accStateMax.pmt_s20_3350) {
      continue;
    }
    if (detEvent->epics_frame.pmt_s20_3360 < accStateMin.pmt_s20_3360 ||
        detEvent->epics_frame.pmt_s20_3360 > accStateMax.pmt_s20_3360) {
      continue;
    }
    if (detEvent->epics_frame.toro_li20_1988_tmit_pc <
            accStateMin.toro_li20_1988_tmit_pc ||
        detEvent->epics_frame.toro_li20_1988_tmit_pc >
            accStateMax.toro_li20_1988_tmit_pc) {
      continue;
    }
    if (detEvent->epics_frame.toro_li20_2040_tmit_pc <
            accStateMin.toro_li20_2040_tmit_pc ||
        detEvent->epics_frame.toro_li20_2040_tmit_pc >
            accStateMax.toro_li20_2040_tmit_pc) {
      continue;
    }
    if (detEvent->epics_frame.toro_li20_2452_tmit_pc <
            accStateMin.toro_li20_2452_tmit_pc ||
        detEvent->epics_frame.toro_li20_2452_tmit_pc >
            accStateMax.toro_li20_2452_tmit_pc) {
      continue;
    }
    if (detEvent->epics_frame.toro_li20_3163_tmit_pc <
            accStateMin.toro_li20_3163_tmit_pc ||
        detEvent->epics_frame.toro_li20_3163_tmit_pc >
            accStateMax.toro_li20_3163_tmit_pc) {
      continue;
    }
    if (detEvent->epics_frame.toro_li20_3255_tmit_pc <
            accStateMin.toro_li20_3255_tmit_pc ||
        detEvent->epics_frame.toro_li20_3255_tmit_pc >
            accStateMax.toro_li20_3255_tmit_pc) {
      continue;
    }
    if (detEvent->epics_frame.radm_li20_1_ch01_meas <
            accStateMin.radm_li20_1_ch01_meas ||
        detEvent->epics_frame.radm_li20_1_ch01_meas >
            accStateMax.radm_li20_1_ch01_meas) {
      continue;
    }

    // PMT_S20_3060
    pmt_s20_3060_hist.Fill(detEvent->epics_frame.pmt_s20_3060);
    // PMT_S20_3070
    pmt_s20_3070_hist.Fill(detEvent->epics_frame.pmt_s20_3070);
    // PMT_S20_3179
    pmt_s20_3179_hist.Fill(detEvent->epics_frame.pmt_s20_3179);
    // PMT_S20_3350
    pmt_s20_3350_hist.Fill(detEvent->epics_frame.pmt_s20_3350);
    // PMT_S20_3360
    pmt_s20_3360_hist.Fill(detEvent->epics_frame.pmt_s20_3360);
    // TORO:LI20:1988:TMIT_PC
    toro_li20_1988_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_1988_tmit_pc);
    // TORO:LI20:2040:TMIT_PC
    toro_li20_2040_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_2040_tmit_pc);
    // TORO:LI20:2452:TMIT_PC
    toro_li20_2452_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_2452_tmit_pc);
    // TORO:LI20:3163:TMIT_PC
    toro_li20_3163_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_3163_tmit_pc);
    // TORO:LI20:3255:TMIT_PC
    toro_li20_3255_tmit_pc_hist.Fill(
        detEvent->epics_frame.toro_li20_3255_tmit_pc);
    // RADM:LI20:1:CH01:MEAS
    radm_li20_1_ch01_meas_hist.Fill(
        detEvent->epics_frame.radm_li20_1_ch01_meas);
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
  // RADM:LI20:1:CH01:MEAS
  accStateMean.radm_li20_1_ch01_meas = radm_li20_1_ch01_meas_hist.GetMean();
  accStateStdDev.radm_li20_1_ch01_meas = radm_li20_1_ch01_meas_hist.GetStdDev();

  return {accStateMean, accStateStdDev};
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

void saveGraphsHists(TTree *tree, std::size_t skip, std::size_t stop,
                     std::unordered_set<std::size_t> badEvents) {

  TH2D bpm3156Hist("bpm3156Hist", "bpm3156Hist", 200, -5, 5, 200, -5, 5);
  TH2D bpm3218Hist("bpm3218Hist", "bpm3218Hist", 200, -5, 5, 200, -5, 5);
  TH2D bpm3265Hist("bpm3265Hist", "bpm3265Hist", 200, -5, 5, 200, -5, 5);
  TH2D bpm3315Hist("bpm3315Hist", "bpm3315Hist", 200, -5, 5, 200, -5, 5);

  TGraph gr_PMT_S20_3060_Raw;
  TGraph gr_PMT_S20_3070_Raw;
  TGraph gr_PMT_S20_3179_Raw;
  TGraph gr_PMT_S20_3350_Raw;
  TGraph gr_PMT_S20_3360_Raw;
  TGraph gr_TORO_LI20_1988_TMIT_PC_Raw;
  TGraph gr_TORO_LI20_2040_TMIT_PC_Raw;
  TGraph gr_TORO_LI20_2452_TMIT_PC_Raw;
  TGraph gr_TORO_LI20_3163_TMIT_PC_Raw;

  TGraph gr_PMT_S20_3060_Filtered;
  TGraph gr_PMT_S20_3070_Filtered;
  TGraph gr_PMT_S20_3179_Filtered;
  TGraph gr_PMT_S20_3350_Filtered;
  TGraph gr_PMT_S20_3360_Filtered;
  TGraph gr_TORO_LI20_1988_TMIT_PC_Filtered;
  TGraph gr_TORO_LI20_2040_TMIT_PC_Filtered;
  TGraph gr_TORO_LI20_2452_TMIT_PC_Filtered;
  TGraph gr_TORO_LI20_3163_TMIT_PC_Filtered;

  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);
  std::size_t nEntries =
      std::min(static_cast<std::size_t>(tree->GetEntries()), stop);
  for (std::size_t i = skip; i < nEntries; i++) {
    tree->GetEntry(i);
    gr_PMT_S20_3060_Raw.SetPoint(i - skip, detEvent->trg_n,
                                 detEvent->epics_frame.pmt_s20_3060);
    gr_PMT_S20_3070_Raw.SetPoint(i - skip, detEvent->trg_n,
                                 detEvent->epics_frame.pmt_s20_3070);
    gr_PMT_S20_3179_Raw.SetPoint(i - skip, detEvent->trg_n,
                                 detEvent->epics_frame.pmt_s20_3179);
    gr_PMT_S20_3350_Raw.SetPoint(i - skip, detEvent->trg_n,
                                 detEvent->epics_frame.pmt_s20_3350);
    gr_PMT_S20_3360_Raw.SetPoint(i - skip, detEvent->trg_n,
                                 detEvent->epics_frame.pmt_s20_3360);
    gr_TORO_LI20_1988_TMIT_PC_Raw.SetPoint(
        i - skip, detEvent->trg_n,
        detEvent->epics_frame.toro_li20_1988_tmit_pc);
    gr_TORO_LI20_2040_TMIT_PC_Raw.SetPoint(
        i - skip, detEvent->trg_n,
        detEvent->epics_frame.toro_li20_2040_tmit_pc);
    gr_TORO_LI20_2452_TMIT_PC_Raw.SetPoint(
        i - skip, detEvent->trg_n,
        detEvent->epics_frame.toro_li20_2452_tmit_pc);
    gr_TORO_LI20_3163_TMIT_PC_Raw.SetPoint(
        i - skip, detEvent->trg_n,
        detEvent->epics_frame.toro_li20_3163_tmit_pc);

    if (!badEvents.contains(i)) {
      bpm3156Hist.Fill(detEvent->epics_frame.bpm_pb_3156_x,
                       detEvent->epics_frame.bpm_pb_3156_y);
      bpm3218Hist.Fill(detEvent->epics_frame.bpm_quad0_3218_x,
                       detEvent->epics_frame.bpm_quad0_3218_y);
      bpm3265Hist.Fill(detEvent->epics_frame.bpm_quad1_3265_x,
                       detEvent->epics_frame.bpm_quad1_3265_y);
      bpm3315Hist.Fill(detEvent->epics_frame.bpm_quad2_3315_x,
                       detEvent->epics_frame.bpm_quad2_3315_y);

      gr_PMT_S20_3060_Filtered.SetPoint(i - skip, detEvent->trg_n,
                                        detEvent->epics_frame.pmt_s20_3060);
      gr_PMT_S20_3070_Filtered.SetPoint(i - skip, detEvent->trg_n,
                                        detEvent->epics_frame.pmt_s20_3070);
      gr_PMT_S20_3179_Filtered.SetPoint(i - skip, detEvent->trg_n,
                                        detEvent->epics_frame.pmt_s20_3179);
      gr_PMT_S20_3350_Filtered.SetPoint(i - skip, detEvent->trg_n,
                                        detEvent->epics_frame.pmt_s20_3350);
      gr_PMT_S20_3360_Filtered.SetPoint(i - skip, detEvent->trg_n,
                                        detEvent->epics_frame.pmt_s20_3360);
      gr_TORO_LI20_1988_TMIT_PC_Filtered.SetPoint(
          i - skip, detEvent->trg_n,
          detEvent->epics_frame.toro_li20_1988_tmit_pc);
      gr_TORO_LI20_2040_TMIT_PC_Filtered.SetPoint(
          i - skip, detEvent->trg_n,
          detEvent->epics_frame.toro_li20_2040_tmit_pc);
      gr_TORO_LI20_2452_TMIT_PC_Filtered.SetPoint(
          i - skip, detEvent->trg_n,
          detEvent->epics_frame.toro_li20_2452_tmit_pc);
      gr_TORO_LI20_3163_TMIT_PC_Filtered.SetPoint(
          i - skip, detEvent->trg_n,
          detEvent->epics_frame.toro_li20_3163_tmit_pc);
    } else {
      gr_PMT_S20_3060_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_PMT_S20_3070_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_PMT_S20_3179_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_PMT_S20_3350_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_PMT_S20_3360_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_TORO_LI20_1988_TMIT_PC_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_TORO_LI20_2040_TMIT_PC_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_TORO_LI20_2452_TMIT_PC_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
      gr_TORO_LI20_3163_TMIT_PC_Filtered.SetPoint(i - skip, detEvent->trg_n, 0);
    }
  }

  LegendPainterConfig legCfg;
  legCfg.draw = true;
  legCfg.option = "p";

  GraphPainterConfig gr_PMT_S20_3060_Raw_Cfg;
  gr_PMT_S20_3060_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_PMT_S20_3060_Filtered_Cfg;
  gr_PMT_S20_3060_Filtered_Cfg.option = "P";
  gr_PMT_S20_3060_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3060_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_PMT_S20_3060_Raw, gr_PMT_S20_3060_Raw_Cfg},
             {gr_PMT_S20_3060_Filtered, gr_PMT_S20_3060_Filtered_Cfg}},
            "PMT_S20_3060.pdf", legCfg);

  GraphPainterConfig gr_PMT_S20_3070_Raw_Cfg;
  gr_PMT_S20_3070_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_PMT_S20_3070_Filtered_Cfg;
  gr_PMT_S20_3070_Filtered_Cfg.option = "P";
  gr_PMT_S20_3070_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3070_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_PMT_S20_3070_Raw, gr_PMT_S20_3070_Raw_Cfg},
             {gr_PMT_S20_3070_Filtered, gr_PMT_S20_3070_Filtered_Cfg}},
            "PMT_S20_3070.pdf", legCfg);

  GraphPainterConfig gr_PMT_S20_3179_Raw_Cfg;
  gr_PMT_S20_3179_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_PMT_S20_3179_Filtered_Cfg;
  gr_PMT_S20_3179_Filtered_Cfg.option = "P";
  gr_PMT_S20_3179_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3179_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_PMT_S20_3179_Raw, gr_PMT_S20_3179_Raw_Cfg},
             {gr_PMT_S20_3179_Filtered, gr_PMT_S20_3179_Filtered_Cfg}},
            "PMT_S20_3179.pdf", legCfg);

  GraphPainterConfig gr_PMT_S20_3350_Raw_Cfg;
  gr_PMT_S20_3350_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_PMT_S20_3350_Filtered_Cfg;
  gr_PMT_S20_3350_Filtered_Cfg.option = "P";
  gr_PMT_S20_3350_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3350_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_PMT_S20_3350_Raw, gr_PMT_S20_3350_Raw_Cfg},
             {gr_PMT_S20_3350_Filtered, gr_PMT_S20_3350_Filtered_Cfg}},
            "PMT_S20_3350.pdf", legCfg);

  GraphPainterConfig gr_PMT_S20_3360_Raw_Cfg;
  gr_PMT_S20_3070_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_PMT_S20_3360_Filtered_Cfg;
  gr_PMT_S20_3360_Filtered_Cfg.option = "P";
  gr_PMT_S20_3360_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3360_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_PMT_S20_3360_Raw, gr_PMT_S20_3360_Raw_Cfg},
             {gr_PMT_S20_3360_Filtered, gr_PMT_S20_3360_Filtered_Cfg}},
            "PMT_S20_3360.pdf", legCfg);

  legCfg.xBottomLeft = 0.2;
  legCfg.yBottomLeft = 0.6;
  legCfg.xTopRight = 0.4;
  legCfg.yTopRight = 0.8;

  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.xTitle = "Trig. N";
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.yTitle = "pC";
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_TORO_LI20_1988_TMIT_PC_Raw, gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_1988_TMIT_PC_Filtered,
              gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg}},
            "TORO_LI20_1988_TMIT_PC.pdf", legCfg);

  legCfg.draw = false;
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.xTitle = "Trig. N";
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.yTitle = "pC";
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_TORO_LI20_2040_TMIT_PC_Raw, gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2040_TMIT_PC_Filtered,
              gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg}},
            "TORO_LI20_2040_TMIT_PC.pdf", legCfg);

  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.xTitle = "Trig. N";
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.yTitle = "pC";
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_TORO_LI20_2452_TMIT_PC_Raw, gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2452_TMIT_PC_Filtered,
              gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg}},
            "TORO_LI20_2452_TMIT_PC.pdf", legCfg);

  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.xTitle = "Trig. N";
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.yTitle = "pC";
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph({{gr_TORO_LI20_3163_TMIT_PC_Raw, gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_3163_TMIT_PC_Filtered,
              gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg}},
            "TORO_LI20_3163_TMIT_PC.pdf", legCfg);

  TCanvas c23("c23", "c23", 800, 800);
  c23.SetGridx();
  c23.SetGridy();
  c23.SetCanvasSize(2000, 2000);
  bpm3156Hist.Draw("colz");
  bpm3156Hist.SetTitle("BPM_3156;Trg N;");
  c23.SaveAs("BPM_3156.pdf");

  TCanvas c20("c20", "c20", 800, 800);
  c20.SetGridx();
  c20.SetGridy();
  c20.SetCanvasSize(2000, 2000);
  bpm3218Hist.Draw("colz");
  bpm3218Hist.SetTitle("BPM_3218;Trg N;");
  c20.SaveAs("BPM_3218.pdf");

  TCanvas c21("c21", "c21", 800, 800);
  c21.SetGridx();
  c21.SetGridy();
  c21.SetCanvasSize(2000, 2000);
  bpm3265Hist.Draw("colz");
  bpm3265Hist.SetTitle("BPM_3265;Trg N;");
  c21.SaveAs("BPM_3265.pdf");

  TCanvas c22("c22", "c22", 800, 800);
  c22.SetGridx();
  c22.SetGridy();
  c22.SetCanvasSize(2000, 2000);
  bpm3315Hist.Draw("colz");
  bpm3315Hist.SetTitle("BPM_3315;Trg N;");
  c22.SaveAs("BPM_3315.pdf");
}
