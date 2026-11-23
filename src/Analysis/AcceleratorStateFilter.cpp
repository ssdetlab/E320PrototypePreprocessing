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
  bx.values.reserve(m_cfg.pvs.size());
  try {
    for (const auto &pv : m_cfg.pvs) {
      bx.values.push_back(getPv(detEvent, pv.key));
    }
    bx.hasPvs = true;
  } catch (std::out_of_range &err) {
    std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
  }

  return bx;
}

std::pair<std::vector<double>, std::vector<double>>
AcceleratorStateFilter::getAcceleratorStateMeanStd(
    const std::vector<BxData> &bxs) const {
  std::vector<TH1D> hists;
  hists.reserve(m_cfg.pvs.size());
  for (const auto &pv : m_cfg.pvs) {
    hists.emplace_back(pv.name.c_str(), pv.name.c_str(), pv.nBins, pv.low,
                       pv.high);
  }

  // Quantiles of the unshuttered BXs
  for (const auto &bx : bxs) {
    if (!bx.hasPid || bx.isShuttered) {
      continue;
    }

    // Fill the PVs read before the first missing one
    for (std::size_t j = 0; j < bx.values.size(); j++) {
      hists.at(j).Fill(bx.values.at(j));
    }
  }

  double quantilePos[2] = {m_cfg.noiseQuantileLow, m_cfg.noiseQuantileHigh};
  std::vector<double> accStateMin(m_cfg.pvs.size());
  std::vector<double> accStateMax(m_cfg.pvs.size());
  for (std::size_t j = 0; j < hists.size(); j++) {
    double qs[2];
    hists.at(j).GetQuantiles(2, qs, quantilePos);
    accStateMin.at(j) = qs[0];
    accStateMax.at(j) = qs[1];
    hists.at(j).Reset();
  }

  // Means and std devs of the unshuttered BXs
  // inside the quantile range of every PV
  for (const auto &bx : bxs) {
    if (!bx.hasPid || bx.isShuttered) {
      continue;
    }
    if (!bx.hasPvs) {
      continue;
    }

    bool isInside = true;
    for (std::size_t j = 0; j < bx.values.size(); j++) {
      if (bx.values.at(j) < accStateMin.at(j) ||
          bx.values.at(j) > accStateMax.at(j)) {
        isInside = false;
        break;
      }
    }
    if (!isInside) {
      continue;
    }

    for (std::size_t j = 0; j < bx.values.size(); j++) {
      hists.at(j).Fill(bx.values.at(j));
    }
  }

  std::vector<double> accStateMean(m_cfg.pvs.size());
  std::vector<double> accStateStdDev(m_cfg.pvs.size());
  for (std::size_t j = 0; j < hists.size(); j++) {
    accStateMean.at(j) = hists.at(j).GetMean();
    accStateStdDev.at(j) = hists.at(j).GetStdDev();
  }

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

    for (std::size_t j = 0; j < bx.values.size(); j++) {
      double mean = accStateMean.at(j);
      double std = accStateStdDev.at(j);
      double val = bx.values.at(j);

      bool isBad = false;
      switch (m_cfg.pvs.at(j).cutType) {
        case Upper: {
          isBad = isBadUpper(mean, std, val);
          break;
        }
        case Central: {
          isBad = isBadCentral(mean, std, val);
          break;
        }
        default: {
          throw std::runtime_error("Unknown PV cut type");
          break;
        }
      }
      if (isBad) {
        badEvents.insert(bx.entry);
        break;
      }
    }
  }
  std::cout << "Found bad events " << badEvents.size() << "\n";
  std::cout << "\n\n\n\n";
  return badEvents;
}
