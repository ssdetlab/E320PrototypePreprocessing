#include "Analysis/PreprocessorMap.hpp"

#include <cmath>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>

#include "Analysis/PVList.hpp"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "detail/EpicsUtils.hpp"
#include "detail/Painters.hpp"

std::map<std::string, PathCollection> Preprocessor::collectDataPaths(
    const Options &opt) const {
  std::cout << "Collecting data paths from\n" << opt.inDataDir << "\n";

  // Collect directories with run data
  std::vector<std::filesystem::path> inDirs;
  for (const auto &entry : std::filesystem::directory_iterator(opt.inDataDir)) {
    if (!entry.is_directory()) {
      continue;
    }
    inDirs.push_back(entry.path());
  }
  std::sort(inDirs.begin(), inDirs.end());
  std::cout << "Runs found: " << inDirs.size() << "\n";

  // Collect paths
  std::map<std::string, PathCollection> paths;
  std::size_t idx = 0;
  for (const auto &inPath : inDirs) {
    std::cout << "Collecting data paths from\n" << inPath << "\n";
    auto beginIt = inPath.begin();
    std::size_t nPathSteps = std::distance(beginIt, inPath.end());
    std::advance(beginIt, nPathSteps - 1);
    std::string runString = beginIt->string();

    paths[runString].parentPath = inPath.string();
    paths[runString].runPrefix = runString;
    paths[runString].outPath = opt.outDir + "/" + runString + ".root";
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(inPath)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      paths[runString].inPaths.push_back(entry.path().string());
    }
    paths[runString].skip = 0;
    paths[runString].stop = std::numeric_limits<std::size_t>::max();
    std::sort(paths[runString].inPaths.begin(), paths[runString].inPaths.end());

    // Collect the orbit reconstruction paths of the run
    std::filesystem::path orbitDir = opt.inOrbitDir + "/" + runString;
    if (!std::filesystem::is_directory(orbitDir)) {
      throw std::runtime_error("Missing orbit reconstruction directory " +
                               orbitDir.string());
    }
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(orbitDir)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      paths[runString].inOrbitPaths.push_back(entry.path().string());
    }
    if (paths[runString].inOrbitPaths.empty()) {
      throw std::runtime_error("No orbit reconstruction files in " +
                               orbitDir.string());
    }
    std::sort(paths[runString].inOrbitPaths.begin(),
              paths[runString].inOrbitPaths.end());
    idx++;
    if (idx > 2) {
      // break;
    }
  }

  // Print debug info
  std::cout << "All paths:\n";
  for (const auto &[runPrefix, collection] : paths) {
    std::cout << "---------------------------------------------\n";
    std::cout << "Run Prefix: " << collection.runPrefix << "\n";
    std::cout << "Out path: " << collection.outPath << "\n";
    std::cout << "In paths:\n";
    for (const auto &inPath : collection.inPaths) {
      std::cout << inPath << "\n";
    }
    std::cout << "In orbit paths:\n";
    for (const auto &inPath : collection.inOrbitPaths) {
      std::cout << inPath << "\n";
    }
    std::cout << "---------------------------------------------\n";
  }
  std::cout << "\n\n\n\n";
  return paths;
};

Preprocessor::FileHandle Preprocessor::getFileHandle(
    const std::string &treeName, const std::vector<std::string> &paths) const {
  TFile *file = nullptr;
  TTree *tree = nullptr;
  TChain *chainOwner = new TChain(treeName.c_str(), treeName.c_str());
  if (paths.size() == 1) {
    file = new TFile(paths.at(0).c_str());
    tree = file->Get<TTree>(treeName.c_str());
  } else {
    chainOwner = new TChain(treeName.c_str());
    for (const auto &path : paths) {
      chainOwner->Add(path.c_str());
    }
    tree = dynamic_cast<TTree *>(chainOwner);
  }

  return {file, chainOwner, tree};
}

void Preprocessor::processFiles(const Options &opt) const {
  using namespace EPICS;

  // Accelerator state PV labels
  std::vector<std::string> accPvNames;
  for (const auto &pv : opt.acceleratorStateFilter->getPvs()) {
    accPvNames.push_back(pv.name);
  }

  Stats globalStats;
  globalStats.init("Total", opt.chipIds, accPvNames);

  auto paths = collectDataPaths(opt);

  ULong64_t eventId = 0;
  for (const auto &[runPrefix, collection] : paths) {
    std::cout << "Processing " << runPrefix << "\n";

    const auto &rois = opt.roiProvider->getROIs(runPrefix);

    // --------------------------------------------------
    // Initialize stats
    Stats localStats;
    localStats.init(runPrefix, opt.chipIds, accPvNames);

    // --------------------------------------------------
    // Get tracker data tree
    auto [inTrackerFile, inTrackerChain, inTrackerTree] =
        getFileHandle(opt.inDataTreeName, collection.inPaths);

    // --------------------------------------------------
    // Initialize tracker branches
    detector_event *detEvent = nullptr;
    inTrackerTree->SetBranchAddress("event", &detEvent);

    // --------------------------------------------------
    // Get orbit reconstruction tree
    auto [inOrbitFile, inOrbitChain, inOrbitTree] =
        getFileHandle(opt.inOrbitTreeName, collection.inOrbitPaths);

    // --------------------------------------------------
    // Initialize orbit reconstruction branches

    // EUDAQ event identifiers
    UInt_t orbitRunNumber = 0;
    ULong64_t orbitTrgN = 0;

    // Reconstructed beam position at the reference
    TVector3 *referencePositionEst = nullptr;

    // Reconstructed beam momentum at the reference
    TLorentzVector *referenceMomentumEst = nullptr;

    // Covariance of the orbit fit at the reference
    TMatrixD *referenceCov = nullptr;

    // Covariance of the beam bunch parameters at the reference
    TMatrixD *referenceBeamBunchCov = nullptr;

    // Covariance at the reference due to the alignment uncertainty
    TMatrixD *referenceAlignmentCov = nullptr;

    // Covariance at the reference of the gradient dependent systematic
    TMatrixD *referenceGradientSystematicCov = nullptr;

    // Shifts at the reference for the alignment samples
    std::vector<TVectorD> *referenceAlignmentSampleShifts = nullptr;

    // Covariance at the reference from the alignment samples
    TMatrixD *referenceAlignmentSampleCov = nullptr;

    // chi2 of the orbit fit
    double referenceFitChi2Est = 0;

    inOrbitTree->SetBranchAddress("runNumber", &orbitRunNumber);
    inOrbitTree->SetBranchAddress("trgN", &orbitTrgN);
    inOrbitTree->SetBranchAddress("referencePositionEst",
                                  &referencePositionEst);
    inOrbitTree->SetBranchAddress("referenceMomentumEst",
                                  &referenceMomentumEst);
    inOrbitTree->SetBranchAddress("referenceCov", &referenceCov);
    inOrbitTree->SetBranchAddress("referenceBeamBunchCov",
                                  &referenceBeamBunchCov);
    inOrbitTree->SetBranchAddress("referenceAlignmentCov",
                                  &referenceAlignmentCov);
    inOrbitTree->SetBranchAddress("referenceGradientSystematicCov",
                                  &referenceGradientSystematicCov);
    inOrbitTree->SetBranchAddress("referenceAlignmentSampleShifts",
                                  &referenceAlignmentSampleShifts);
    inOrbitTree->SetBranchAddress("referenceAlignmentSampleCov",
                                  &referenceAlignmentSampleCov);
    inOrbitTree->SetBranchAddress("referenceFitChi2Est", &referenceFitChi2Est);

    // Map the EUDAQ event identifiers to the orbit entries
    std::map<std::pair<std::uint32_t, std::uint64_t>, std::size_t> orbitEntries;
    for (std::size_t i = 0; i < inOrbitTree->GetEntries(); i++) {
      inOrbitTree->GetEntry(i);
      auto [it, inserted] =
          orbitEntries.insert({{orbitRunNumber, orbitTrgN}, i});
      if (!inserted) {
        throw std::runtime_error("Duplicate orbit reconstruction for run " +
                                 std::to_string(orbitRunNumber) + " trigger " +
                                 std::to_string(orbitTrgN));
      }
    }
    std::size_t nEventsNoOrbit = 0;

    // --------------------------------------------------
    // Output files
    if (auto parentPath =
            std::filesystem::path(collection.outPath).parent_path();
        !std::filesystem::exists(parentPath)) {
      std::filesystem::create_directory(parentPath);
    }
    TFile *outFile = new TFile(collection.outPath.c_str(), "RECREATE");
    TTree *outDataTree =
        new TTree(opt.outDataTreeName.c_str(), opt.outDataTreeName.c_str());
    TTree *outMagFieldTree =
        new TTree(opt.outMagnetTreeName.c_str(), opt.outMagnetTreeName.c_str());

    // Output data tree branches
    E320Io::DetectorEvent outDetEvent;
    outDataTree->Branch("event", &outDetEvent);
    outDataTree->Branch("eventId", &eventId);

    // Output orbit reconstruction branches, the orbit prefix marks the
    // quantities of the reconstructed beam orbit
    TVector3 outOrbitReferencePositionEst;
    TLorentzVector outOrbitReferenceMomentumEst;
    TMatrixD outOrbitReferenceCov;
    TMatrixD outOrbitReferenceBeamBunchCov;
    TMatrixD outOrbitReferenceAlignmentCov;
    TMatrixD outOrbitReferenceGradientSystematicCov;
    std::vector<TVectorD> outOrbitReferenceAlignmentSampleShifts;
    TMatrixD outOrbitReferenceAlignmentSampleCov;
    double outOrbitReferenceFitChi2Est = 0;

    outDataTree->Branch("orbitReferencePositionEst",
                        &outOrbitReferencePositionEst);
    outDataTree->Branch("orbitReferenceMomentumEst",
                        &outOrbitReferenceMomentumEst);
    outDataTree->Branch("orbitReferenceCov", &outOrbitReferenceCov);
    outDataTree->Branch("orbitReferenceBeamBunchCov",
                        &outOrbitReferenceBeamBunchCov);
    outDataTree->Branch("orbitReferenceAlignmentCov",
                        &outOrbitReferenceAlignmentCov);
    outDataTree->Branch("orbitReferenceGradientSystematicCov",
                        &outOrbitReferenceGradientSystematicCov);
    outDataTree->Branch("orbitReferenceAlignmentSampleShifts",
                        &outOrbitReferenceAlignmentSampleShifts);
    outDataTree->Branch("orbitReferenceAlignmentSampleCov",
                        &outOrbitReferenceAlignmentSampleCov);
    outDataTree->Branch("orbitReferenceFitChi2Est",
                        &outOrbitReferenceFitChi2Est);

    // Output magnet tree branches
    double quadTripletFocusingEnergy = 0;
    double quadTripletObjectZ = 0;
    double quadTripletImageZ = 0;
    double quadTripletM12 = 0;
    double quadTripletM34 = 0;
    double quad1Grad = 0;
    double quad2Grad = 0;
    double quad3Grad = 0;
    double xCorrectorStrength = 0;
    double dipoleStrength = 0;

    outMagFieldTree->Branch("quadTripletFocusingEnergy",
                            &quadTripletFocusingEnergy);
    outMagFieldTree->Branch("quadTripletObjectZ", &quadTripletObjectZ);
    outMagFieldTree->Branch("quadTripletImageZ", &quadTripletImageZ);
    outMagFieldTree->Branch("quadTripletM12", &quadTripletM12);
    outMagFieldTree->Branch("quadTripletM34", &quadTripletM34);
    outMagFieldTree->Branch("quad1Grad", &quad1Grad);
    outMagFieldTree->Branch("quad2Grad", &quad2Grad);
    outMagFieldTree->Branch("quad3Grad", &quad3Grad);
    outMagFieldTree->Branch("xCorrectorStrength", &xCorrectorStrength);
    outMagFieldTree->Branch("dipoleStrength", &dipoleStrength);
    outMagFieldTree->Branch("eventId", &eventId);

    // --------------------------------------------------
    // Collect the accelerator data of the run BXs
    std::size_t nEntries = std::min(
        static_cast<std::size_t>(inTrackerTree->GetEntries()), collection.stop);
    std::cout << "Collecting accelerator data\n";
    std::vector<AcceleratorStateFilter::BxData> bxs;
    for (std::size_t i = collection.skip; i < nEntries; i++) {
      inTrackerTree->GetEntry(i);
      bxs.push_back(opt.acceleratorStateFilter->getBxData(detEvent, i));
    }

    // --------------------------------------------------
    // Collect bad events
    std::unordered_set<std::size_t> badEvents =
        opt.acceleratorStateFilter->getBadEvents(bxs);

    // --------------------------------------------------
    // Fill in accelerator events
    std::cout << "Filling in accelerator data\n";
    for (const auto &bx : bxs) {
      // Shuttered BXs are excluded from the accelerator diagnostics
      if (!bx.hasPid || bx.isShuttered) {
        continue;
      }

      // Raw values of the PVs read before the first missing one
      for (std::size_t j = 0; j < bx.values.size(); j++) {
        auto &gr = localStats.accRawGrs.at(j);
        gr.SetPoint(gr.GetN(), bx.trgN, bx.values.at(j));
      }
      if (!bx.hasPvs) {
        continue;
      }

      if (badEvents.contains(bx.entry)) {
        continue;
      }
      for (std::size_t j = 0; j < bx.values.size(); j++) {
        auto &gr = localStats.accFilteredGrs.at(j);
        gr.SetPoint(gr.GetN(), bx.trgN, bx.values.at(j));
      }
    }

    // --------------------------------------------------
    // Fill in the tracker data
    std::cout << "Filling in data\n";
    for (std::size_t i = collection.skip; i < nEntries; i++) {
      // Total events
      localStats.nEventsTotal++;
      globalStats.nEventsTotal++;

      // Bad and accepted events
      if (badEvents.contains(i)) {
        continue;
      }
      localStats.nEventsAccepted++;
      globalStats.nEventsAccepted++;
      inTrackerTree->GetEntry(i);

      std::size_t epicsPID = 0;
      try {
        epicsPID = getPv(detEvent, epics_pid_key);
      } catch (std::out_of_range &err) {
        continue;
      }

      // Shuttered events
      std::size_t pid =
          static_cast<std::size_t>(std::round((epicsPID - 10) / 36));
      bool isShuttered = pid % 20 == 0;
      if (isShuttered && opt.excludeShuttered) {
        continue;
      }
      localStats.nEventsUnshuttered++;
      globalStats.nEventsUnshuttered++;

      // Even/odd PID events
      bool isEven = pid % 2 == 0;
      if (isEven) {
        localStats.nEventsEvenPID++;
        globalStats.nEventsEvenPID++;
      } else {
        localStats.nEventsOddPID++;
        globalStats.nEventsOddPID++;
      }

      outDetEvent.trgN = detEvent->trg_n;
      outDetEvent.tsBegin = detEvent->ts_begin;
      outDetEvent.tsEnd = detEvent->ts_end;
      outDetEvent.runNumber = detEvent->run_number;
      outDetEvent.epicsParity = isEven ? 0 : 1;
      outDetEvent.epicsPID = std::round((epicsPID - 10) / 36) * 36 + 10;
      outDetEvent.epicsDAQNumber = getPv(detEvent, epics_daq_number_key);

      std::vector<E320Io::StaveEvent> outStaveEvents;
      outStaveEvents.reserve(detEvent->st_ev_buffer.size());

      quad1Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3141_bact));
      quad2Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3261_bact));
      quad3Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3091_bact));
      dipoleStrength =
          -convertEpicsDipole(getPv(detEvent, li20_lgps_3330_bact) * 1e9,
                              opt.deflectionAngle, opt.dipoleLength);
      xCorrectorStrength = convertEpicsXCorrector(
          getPv(detEvent, li20_xcor_3276_bact), opt.xCorrectorLength);
      quadTripletFocusingEnergy = getPv(detEvent, quad_focusing_energy_key);
      quadTripletObjectZ = getPv(detEvent, quad_object_z_key);
      quadTripletImageZ = getPv(detEvent, quad_image_z_key);
      quadTripletM12 = getPv(detEvent, quad_M12_key);
      quadTripletM34 = getPv(detEvent, quad_M34_key);

      for (const auto &staveEv : detEvent->st_ev_buffer) {
        E320Io::StaveEvent outStaveEvent;

        outStaveEvent.staveId = staveEv.stave_id;
        std::vector<E320Io::ChipEvent> outChipEvents;
        outChipEvents.reserve(staveEv.ch_ev_buffer.size());
        for (const auto &chipEv : staveEv.ch_ev_buffer) {
          int id = 100 * detEvent->run_number + 10 * staveEv.stave_id +
                   chipEv.chip_id;
          E320Io::ChipEvent outChipEvent;

          // ALPIDE readout flags
          outChipEvent.isBusyViolation = chipEv.is_busy_violation;
          outChipEvent.isFlushedIncomplete = chipEv.is_flushed_incomplete;
          outChipEvent.isStrobeExtended = chipEv.is_strobe_extended;
          outChipEvent.isBusyTransition = chipEv.is_busy_transition;

          // MOSAIC readout flags
          outChipEvent.endOfRun = chipEv.end_of_run;
          outChipEvent.overflow = chipEv.overflow;
          outChipEvent.timeout = chipEv.timeout;
          outChipEvent.headerError = chipEv.header_error;
          outChipEvent.decoder10b8bError = chipEv.decoder_10b8b_error;
          outChipEvent.eventOversizeError = chipEv.event_oversize_error;

          outChipEvent.chipId = chipEv.chip_id;
          outChipEvent.channel = chipEv.channel;

          // Collect hits
          std::unordered_set<Hit, PairHash> hits;
          for (const auto &hit : chipEv.hits) {
            hits.insert(hit);
          }

          // Get clusters
          auto clusters = opt.clusteringAlgo->getClusters(hits);

          // Continue if the chip is overpopulated
          if (clusters.size() > opt.maxHits) {
            continue;
          }

          // Fill non-ROId hit map
          localStats.chipOccRawGrs.at(outChipEvent.chipId)
              .SetPoint(localStats.chipOccRawGrs.at(outChipEvent.chipId).GetN(),
                        eventId, clusters.size());
          for (const auto &hit : clusters) {
            localStats.chipHitMapRawHists.at(chipEv.chip_id)
                .Fill(hit.xCenter, hit.yCenter);
          }

          // Pass hits through the ROI
          outChipEvent.hits.reserve(clusters.size());
          for (const auto &cl : clusters) {
            if (rois.at(outChipEvent.chipId)
                    ->isInside(cl.xCenter, cl.yCenter)) {
              outChipEvent.hits.push_back(cl);
            }
          }
          outChipEvent.hits.shrink_to_fit();

          // Continue if the ROI is overpopulated
          if (outChipEvent.hits.size() > opt.maxRoiHits) {
            continue;
          }
          // Continue if the hits are empty
          if (outChipEvent.hits.empty()) {
            continue;
          }

          // Fill ROId hit map
          localStats.chipOccFilteredGrs.at(outChipEvent.chipId)
              .SetPoint(
                  localStats.chipOccFilteredGrs.at(outChipEvent.chipId).GetN(),
                  eventId, outChipEvent.hits.size());
          for (const auto &hit : outChipEvent.hits) {
            localStats.chipHitMapFilteredHists.at(outChipEvent.chipId)
                .Fill(hit.xCenter, hit.yCenter);
          }

          outChipEvents.push_back(std::move(outChipEvent));
        }

        // Continue if the stave is not full
        if (outChipEvents.size() != opt.chipIds.size()) {
          continue;
        }
        outStaveEvent.chipEventBuffer = outChipEvents;
        outStaveEvents.push_back(std::move(outStaveEvent));
      }

      // Continue if the detector is not full
      if (!outStaveEvents.empty()) {
        if (isEven) {
          localStats.nEventsFullDetEvenPID++;
          globalStats.nEventsFullDetEvenPID++;
        } else {
          localStats.nEventsFullDetOddPID++;
          globalStats.nEventsFullDetOddPID++;
        }
      } else {
        continue;
      }

      // Continue if the event has no orbit reconstruction
      auto orbitIt = orbitEntries.find({detEvent->run_number, detEvent->trg_n});
      if (orbitIt == orbitEntries.end()) {
        nEventsNoOrbit++;
        continue;
      }
      inOrbitTree->GetEntry(orbitIt->second);

      outOrbitReferencePositionEst = *referencePositionEst;
      outOrbitReferenceMomentumEst = *referenceMomentumEst;
      outOrbitReferenceCov.ResizeTo(*referenceCov);
      outOrbitReferenceCov = *referenceCov;
      outOrbitReferenceBeamBunchCov.ResizeTo(*referenceBeamBunchCov);
      outOrbitReferenceBeamBunchCov = *referenceBeamBunchCov;
      outOrbitReferenceAlignmentCov.ResizeTo(*referenceAlignmentCov);
      outOrbitReferenceAlignmentCov = *referenceAlignmentCov;
      outOrbitReferenceGradientSystematicCov.ResizeTo(
          *referenceGradientSystematicCov);
      outOrbitReferenceGradientSystematicCov = *referenceGradientSystematicCov;
      outOrbitReferenceAlignmentSampleShifts.clear();
      outOrbitReferenceAlignmentSampleShifts = *referenceAlignmentSampleShifts;
      outOrbitReferenceAlignmentSampleCov.ResizeTo(
          *referenceAlignmentSampleCov);
      outOrbitReferenceAlignmentSampleCov = *referenceAlignmentSampleCov;
      outOrbitReferenceFitChi2Est = referenceFitChi2Est;

      outDetEvent.staveEventBuffer = outStaveEvents;
      outDataTree->Fill();
      outMagFieldTree->Fill();

      eventId++;
      if (i % 100 == 0) {
        std::cout << i << "/" << nEntries << "\n";
      }
    }

    outDataTree->Write();
    outMagFieldTree->Write();
    outFile->Close();

    std::cout << "Events without orbit reconstruction: " << nEventsNoOrbit
              << "\n";

    localStats.print();
    localStats.reset();
    std::cout << "Finished " << runPrefix << "\n";
    std::cout << "\n\n\n\n";
  }

  std::cout << "Finished processing\n";
  globalStats.print();
}

void Preprocessor::Stats::init(const std::string &prefix,
                               const std::set<std::size_t> &ids,
                               const std::vector<std::string> &pvNames) {
  runPrefix = prefix;
  chipIds = ids;
  accPvNames = pvNames;

  std::string nameEventSplit = runPrefix + "EventSplit";
  eventCategotyHist =
      TH1I(nameEventSplit.c_str(), nameEventSplit.c_str(), 6, 0, 6);
  eventCategotyHist.Fill("nEventsTotal", 0);
  eventCategotyHist.Fill("nEventsAccepted", 0);
  eventCategotyHist.Fill("nEventsUnshuttered", 0);
  eventCategotyHist.Fill("nEventsEvenPID", 0);
  eventCategotyHist.Fill("nEventsOddPID", 0);
  eventCategotyHist.Fill("nEventsFullDetEvenPID", 0);
  eventCategotyHist.Fill("nEventsFullDetOddPID", 0);

  for (int id : chipIds) {
    // Raw hit map hists
    std::string nameChipHitMapRaw =
        runPrefix + "Chip" + std::to_string(id) + "HitMapRaw";
    chipHitMapRawHists.insert(
        {id, TH2D(nameChipHitMapRaw.c_str(), nameChipHitMapRaw.c_str(), 256, 0,
                  1024, 128, 0, 512)});
    // Filtered hit map hists
    std::string nameChipHitMapFiltered =
        runPrefix + "Chip" + std::to_string(id) + "HitMapFiltered";
    chipHitMapFilteredHists.insert(
        {id, TH2D(nameChipHitMapFiltered.c_str(),
                  nameChipHitMapFiltered.c_str(), 256, 0, 1024, 128, 0, 512)});
    // Raw occupancy graphs
    std::string nameChipOccR =
        runPrefix + "Chip" + std::to_string(id) + "OccRaw";
    chipOccRawGrs.insert({id, TGraph()});
    // Filtered occupancy graphs
    std::string nameChipOccFilter =
        runPrefix + "Chip" + std::to_string(id) + "OccFiltered";
    chipOccFilteredGrs.insert({id, TGraph()});
  }

  // Accelerator state PV graphs
  accRawGrs.assign(accPvNames.size(), TGraph());
  accFilteredGrs.assign(accPvNames.size(), TGraph());
}

void Preprocessor::Stats::reset() {
  // N total events
  nEventsTotal = 0;
  // N accepted events
  nEventsAccepted = 0;
  // N shuttered events
  nEventsUnshuttered = 0;
  // N even PID events
  nEventsEvenPID = 0;
  // N odd PID events
  nEventsOddPID = 0;
  // N events with full detector occupancy
  // and even PID
  nEventsFullDetEvenPID = 0;
  // N events with full detector occupancy
  // and odd PID
  nEventsFullDetOddPID = 0;

  // Raw hit map hists
  chipHitMapRawHists.clear();
  // Filtered hit map hists
  chipHitMapFilteredHists.clear();
  // Raw occupancy graphs
  chipOccRawGrs.clear();
  // Filtered occupancy graphs
  chipOccFilteredGrs.clear();
}

void Preprocessor::Stats::print() {
  TCanvas c = TCanvas(runPrefix.c_str(), "", 800, 800);
  c.GetPad(1);
  c.SetMargin(0.1, 0.1, 0.1, 0.1);

  c.SetGridx(1);
  c.SetGridy(1);

  c.SetLogx(0);
  c.SetLogy(0);

  LegendPainterConfig legCfg;
  legCfg.draw = true;
  legCfg.option = "p";

  std::string cOptionFirst = runPrefix + ".pdf(";
  std::string cOptionLast = runPrefix + ".pdf)";
  std::string cOption = runPrefix + ".pdf";

  HistPainterConfig eventSplitCfg;
  eventSplitCfg.yTitle = "Count";
  eventSplitCfg.option = "hist text";

  eventCategotyHist.Fill("nEventsTotal", nEventsTotal);
  eventCategotyHist.Fill("nEventsAccepted", nEventsAccepted);
  eventCategotyHist.Fill("nEventsUnshuttered", nEventsUnshuttered);
  eventCategotyHist.Fill("nEventsEvenPID", nEventsEvenPID);
  eventCategotyHist.Fill("nEventsOddPID", nEventsOddPID);
  eventCategotyHist.Fill("nEventsFullDetEvenPID", nEventsFullDetEvenPID);
  eventCategotyHist.Fill("nEventsFullDetOddPID", nEventsFullDetOddPID);
  saveHist(c, cOptionFirst, {{eventCategotyHist, eventSplitCfg}});

  // Accelerator state PVs
  for (std::size_t j = 0; j < accPvNames.size(); j++) {
    std::string grName = "gr_" + accPvNames.at(j);
    GraphPainterConfig grRawCfg;
    grRawCfg.legTitle = "Raw";
    grRawCfg.printTitle = true;
    accRawGrs.at(j).SetNameTitle(grName.c_str(), grName.c_str());
    GraphPainterConfig grFilteredCfg;
    grFilteredCfg.option = "P";
    grFilteredCfg.color = kRed;
    grFilteredCfg.legTitle = "Filtered";
    saveGraph(
        c, cOption,
        {{accRawGrs.at(j), grRawCfg}, {accFilteredGrs.at(j), grFilteredCfg}},
        legCfg);
  }

  GraphPainterConfig grOccRawCfg;
  grOccRawCfg.xTitle = "EventId";
  grOccRawCfg.yTitle = "Clusters";
  grOccRawCfg.legTitle = "Raw";
  grOccRawCfg.printTitle = true;
  GraphPainterConfig grOccFiltered;
  grOccFiltered.option = "P";
  grOccFiltered.color = kRed;
  grOccFiltered.legTitle = "Filtered";
  for (int id : chipIds) {
    auto &grRaw = chipOccRawGrs.at(id);
    std::string grName = runPrefix + "ChipOccRaw" + std::to_string(id);
    grRaw.SetNameTitle(grName.c_str(), grName.c_str());
    saveGraph(
        c, cOption,
        {{grRaw, grOccRawCfg}, {chipOccFilteredGrs.at(id), grOccFiltered}},
        legCfg);
  }

  HistPainterConfig chipHitMapCfg;
  chipHitMapCfg.xTitle = "Pix X";
  chipHitMapCfg.yTitle = "Pix Y";
  chipHitMapCfg.stats = true;
  chipHitMapCfg.option = "colz";
  for (int id : chipIds) {
    saveHist(c, cOption, {{chipHitMapRawHists.at(id), chipHitMapCfg}});
    saveHist(c, cOption, {{chipHitMapFilteredHists.at(id), chipHitMapCfg}});
  }

  HistPainterConfig dummyCfg;
  TH1D dummy("dummy", "dummy", 100, 0, 100);
  saveHist(c, cOptionLast, {{dummy, dummyCfg}});
}
