#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <unordered_set>
#include <vector>

#include "AcceleratorState.hpp"
#include "Analysis/ROI.hpp"
#include "DetectorEvent.h"
#include "TTree.h"
#include "Types/DetectorEvent.hpp"
#include "detail/PairHash.hpp"

using Cluster = E320Io::ChipEvent::Cluster;
using Hit = std::pair<std::uint16_t, std::uint16_t>;
using HitSet = std::unordered_set<Hit, PairHash>;

void bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot, HitSet &pixels);

std::vector<Cluster> getClusters(HitSet &hits, std::size_t maxCluterSize);

std::vector<Cluster> getClusters(HitSet &hits, std::size_t maxCluterSize,
                                 const std::shared_ptr<ROI> &roi);

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
