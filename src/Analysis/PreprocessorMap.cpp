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
#include "ROI/ButterflyROI.hpp"
#include "ROI/EllipseROI.hpp"
#include "ROI/UnboundedROI.hpp"
#include "TCanvas.h"
#include "detail/Painters.hpp"

std::map<std::string, Preprocessor::PathCollection>
Preprocessor::collectDataPaths(const Options &opt) const {
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

std::map<std::size_t, std::shared_ptr<ROI>> Preprocessor::getROIs(
    const std::string &pathToROIData, const Options &opt) const {
  std::map<std::size_t, std::shared_ptr<ROI>> rois;
  switch (opt.roiType) {
    case Ellipse: {
      throw std::runtime_error("Ellipse ROI reading is not supported yet");
      break;
    }
    case Butterfly: {
      ButterflyROIFactory factory;
      auto butterflies = factory.buildROIForRun(pathToROIData);
      for (std::size_t id : opt.chipIds) {
        rois.insert({id, butterflies.at(id)});
      }
      break;
    }
    case Unbounded: {
      for (std::size_t id : opt.chipIds) {
        rois.insert({id, std::make_shared<UnboundedROI>()});
      }
      break;
    }
    default: {
      throw std::runtime_error("Unknown ROI type");
      break;
    }
  }
  return rois;
}

void Preprocessor::processFiles(const Options &opt) const {
  Stats globalStats;
  globalStats.init("Total", opt.chipIds);

  auto paths = collectDataPaths(opt);

  ULong64_t eventId = 0;
  for (const auto &[runPrefix, collection] : paths) {
    std::cout << "Processing " << runPrefix << "\n";

    std::string roiPath = opt.inROIDir + "/" + runPrefix + ".csv";
    auto rois = getROIs(roiPath, opt);

    // Initialize stats
    Stats localStats;
    localStats.init(runPrefix, opt.chipIds);

    // Get data tree
    auto [inFile, inChain, inTree] =
        getFileHandle(opt.inDataTreeName, collection.inPaths);

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

    // Collect bad events
    std::unordered_set<std::size_t> badEvents =
        getBadEvents(inTree, collection.skip, opt);

    // Initialize tree readout
    detector_event *detEvent = nullptr;

    E320Io::DetectorEvent outDetEvent;
    outDataTree->Branch("event", &outDetEvent);
    outDataTree->Branch("eventId", &eventId);

    inTree->SetBranchAddress("event", &detEvent);
    std::size_t nEntries = std::min(
        static_cast<std::size_t>(inTree->GetEntries()), collection.stop);

    double quad1Grad = 0;
    double quad2Grad = 0;
    double quad3Grad = 0;
    double xCorrectorStrength = 0;
    double dipoleStrength = 0;

    outMagFieldTree->Branch("quad1Grad", &quad1Grad);
    outMagFieldTree->Branch("quad2Grad", &quad2Grad);
    outMagFieldTree->Branch("quad3Grad", &quad3Grad);
    outMagFieldTree->Branch("xCorrectorStrength", &xCorrectorStrength);
    outMagFieldTree->Branch("dipoleStrength", &dipoleStrength);
    outMagFieldTree->Branch("eventId", &eventId);

    // Fill in accelerator events
    std::cout << "Filling in accelerator data\n";
    for (std::size_t i = collection.skip; i < nEntries; i++) {
      inTree->GetEntry(i);

      try {
        // PMTs
        localStats.gr_PMT_S20_3060_Raw.SetPoint(
            localStats.gr_PMT_S20_3060_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, pmt_s20_3060_key));
        localStats.gr_PMT_S20_3070_Raw.SetPoint(
            localStats.gr_PMT_S20_3070_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, pmt_s20_3070_key));
        localStats.gr_PMT_S20_3179_Raw.SetPoint(
            localStats.gr_PMT_S20_3179_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, pmt_s20_3179_key));
        localStats.gr_PMT_S20_3350_Raw.SetPoint(
            localStats.gr_PMT_S20_3350_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, pmt_s20_3350_key));
        localStats.gr_PMT_S20_3360_Raw.SetPoint(
            localStats.gr_PMT_S20_3360_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, pmt_s20_3360_key));
        // Toroids
        localStats.gr_TORO_LI20_1988_TMIT_PC_Raw.SetPoint(
            localStats.gr_TORO_LI20_1988_TMIT_PC_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, toro_li20_1988_tmit_pc_key));
        localStats.gr_TORO_LI20_2040_TMIT_PC_Raw.SetPoint(
            localStats.gr_TORO_LI20_2040_TMIT_PC_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, toro_li20_2040_tmit_pc_key));
        localStats.gr_TORO_LI20_2452_TMIT_PC_Raw.SetPoint(
            localStats.gr_TORO_LI20_2452_TMIT_PC_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, toro_li20_2452_tmit_pc_key));
        localStats.gr_TORO_LI20_3163_TMIT_PC_Raw.SetPoint(
            localStats.gr_TORO_LI20_3163_TMIT_PC_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, toro_li20_3163_tmit_pc_key));
        localStats.gr_TORO_LI20_3255_TMIT_PC_Raw.SetPoint(
            localStats.gr_TORO_LI20_3255_TMIT_PC_Raw.GetN(), detEvent->trg_n,
            getPv(detEvent, toro_li20_3255_tmit_pc_key));
      } catch (std::out_of_range &err) {
        std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
        continue;
      }
      if (badEvents.contains(i)) {
        continue;
      }
      // PMTs
      localStats.gr_PMT_S20_3060_Filtered.SetPoint(
          localStats.gr_PMT_S20_3060_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, pmt_s20_3060_key));
      localStats.gr_PMT_S20_3070_Filtered.SetPoint(
          localStats.gr_PMT_S20_3070_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, pmt_s20_3070_key));
      localStats.gr_PMT_S20_3179_Filtered.SetPoint(
          localStats.gr_PMT_S20_3179_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, pmt_s20_3179_key));
      localStats.gr_PMT_S20_3350_Filtered.SetPoint(
          localStats.gr_PMT_S20_3350_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, pmt_s20_3350_key));
      localStats.gr_PMT_S20_3360_Filtered.SetPoint(
          localStats.gr_PMT_S20_3360_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, pmt_s20_3360_key));
      // Toroids
      localStats.gr_TORO_LI20_1988_TMIT_PC_Filtered.SetPoint(
          localStats.gr_TORO_LI20_1988_TMIT_PC_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, toro_li20_1988_tmit_pc_key));
      localStats.gr_TORO_LI20_2040_TMIT_PC_Filtered.SetPoint(
          localStats.gr_TORO_LI20_2040_TMIT_PC_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, toro_li20_2040_tmit_pc_key));
      localStats.gr_TORO_LI20_2452_TMIT_PC_Filtered.SetPoint(
          localStats.gr_TORO_LI20_2452_TMIT_PC_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, toro_li20_2452_tmit_pc_key));
      localStats.gr_TORO_LI20_3163_TMIT_PC_Filtered.SetPoint(
          localStats.gr_TORO_LI20_3163_TMIT_PC_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, toro_li20_3163_tmit_pc_key));
      localStats.gr_TORO_LI20_3255_TMIT_PC_Filtered.SetPoint(
          localStats.gr_TORO_LI20_3255_TMIT_PC_Filtered.GetN(), detEvent->trg_n,
          getPv(detEvent, toro_li20_3255_tmit_pc_key));
    }

    // Fill in the data
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
      inTree->GetEntry(i);

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

      outDetEvent.trg_n = detEvent->trg_n;
      outDetEvent.ts_begin = detEvent->ts_begin;
      outDetEvent.ts_end = detEvent->ts_end;
      outDetEvent.run_number = detEvent->run_number;
      outDetEvent.epicsParity = isEven ? 0 : 1;
      outDetEvent.epicsPID = std::round((epicsPID - 10) / 36) * 36 + 10;
      outDetEvent.epicsDAQNumber = getPv(detEvent, epics_daq_number_key);

      std::vector<E320Io::StaveEvent> outStaveEvents;
      outStaveEvents.reserve(detEvent->st_ev_buffer.size());

      E320Io::BPMEvent bpm3156Ev;
      bpm3156Ev.meanX = getPv(detEvent, bpm_3156_x_key);
      bpm3156Ev.meanY = getPv(detEvent, bpm_3156_y_key);
      bpm3156Ev.stdDevX = 0.1;
      bpm3156Ev.stdDevY = 0.1;
      bpm3156Ev.id = 3156;

      E320Io::BPMEvent bpm3218Ev;
      bpm3218Ev.meanX = getPv(detEvent, bpm_quad0_3218_x_key);
      bpm3218Ev.meanY = getPv(detEvent, bpm_quad0_3218_y_key);
      bpm3218Ev.stdDevX = 0.1;
      bpm3218Ev.stdDevY = 0.1;
      bpm3218Ev.id = 3218;

      E320Io::BPMEvent bpm3265Ev;
      bpm3265Ev.meanX = getPv(detEvent, bpm_quad1_3265_x_key);
      bpm3265Ev.meanY = getPv(detEvent, bpm_quad1_3265_y_key);
      bpm3265Ev.stdDevX = 0.1;
      bpm3265Ev.stdDevY = 0.1;
      bpm3265Ev.id = 3265;

      E320Io::BPMEvent bpm3315Ev;
      bpm3315Ev.meanX = getPv(detEvent, bpm_quad2_3315_x_key);
      bpm3315Ev.meanY = getPv(detEvent, bpm_quad2_3315_y_key);
      bpm3315Ev.stdDevX = 0.1;
      bpm3315Ev.stdDevY = 0.1;
      bpm3315Ev.id = 3315;

      outDetEvent.bpm_ev_buffer = {bpm3156Ev, bpm3218Ev, bpm3265Ev, bpm3315Ev};

      quad1Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3141_bact));
      quad2Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3261_bact));
      quad3Grad = convertEpicsQuad(getPv(detEvent, li20_lgps_3091_bact));
      dipoleStrength =
          -convertEpicsDipole(getPv(detEvent, li20_lgps_3330_bact) * 1e9,
                              opt.deflectionAngle, opt.dipoleLength);
      xCorrectorStrength = convertEpicsXCorrector(
          getPv(detEvent, li20_xcor_3276_bact), opt.xCorrectorLength);

      for (const auto &staveEv : detEvent->st_ev_buffer) {
        E320Io::StaveEvent outStaveEvent;

        outStaveEvent.stave_id = staveEv.stave_id;
        std::vector<E320Io::ChipEvent> outChipEvents;
        outChipEvents.reserve(staveEv.ch_ev_buffer.size());
        for (const auto &chipEv : staveEv.ch_ev_buffer) {
          int id = 100 * detEvent->run_number + 10 * staveEv.stave_id +
                   chipEv.chip_id;
          E320Io::ChipEvent outChipEvent;

          // ALPIDE readout flags
          outChipEvent.is_busy_violation = chipEv.is_busy_violation;
          outChipEvent.is_flushed_incomplete = chipEv.is_flushed_incomplete;
          outChipEvent.is_strobe_extended = chipEv.is_strobe_extended;
          outChipEvent.is_busy_transition = chipEv.is_busy_transition;

          // MOSAIC readout flags
          outChipEvent.end_of_run = chipEv.end_of_run;
          outChipEvent.overflow = chipEv.overflow;
          outChipEvent.timeout = chipEv.timeout;
          outChipEvent.header_error = chipEv.header_error;
          outChipEvent.decoder_10b8b_error = chipEv.decoder_10b8b_error;
          outChipEvent.event_oversize_error = chipEv.event_oversize_error;

          outChipEvent.chip_id = chipEv.chip_id;
          outChipEvent.channel = chipEv.channel;

          // Collect hits
          std::unordered_set<Hit, PairHash> hits;
          for (const auto &hit : chipEv.hits) {
            hits.insert(hit);
          }

          // Get clusters
          auto clusters = getClusters(hits, opt.maxClusterSize);

          // Continue if the chip is overpopulated
          if (clusters.size() > opt.maxHits) {
            continue;
          }

          // Fill non-ROId hit map
          localStats.chipOccRawGrs.at(outChipEvent.chip_id)
              .SetPoint(
                  localStats.chipOccRawGrs.at(outChipEvent.chip_id).GetN(),
                  eventId, clusters.size());
          for (const auto &hit : clusters) {
            localStats.chipHitMapRawHists.at(chipEv.chip_id)
                .Fill(hit.xCenter, hit.yCenter);
          }

          // Pass hits through the ROI
          outChipEvent.hits.reserve(clusters.size());
          for (const auto &cl : clusters) {
            if (rois.at(outChipEvent.chip_id)
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
          localStats.chipOccFilteredGrs.at(outChipEvent.chip_id)
              .SetPoint(
                  localStats.chipOccFilteredGrs.at(outChipEvent.chip_id).GetN(),
                  eventId, outChipEvent.hits.size());
          for (const auto &hit : outChipEvent.hits) {
            localStats.chipHitMapFilteredHists.at(outChipEvent.chip_id)
                .Fill(hit.xCenter, hit.yCenter);
          }

          outChipEvents.push_back(std::move(outChipEvent));
        }

        // Continue if the stave is not full
        if (outChipEvents.size() != opt.chipIds.size()) {
          continue;
        }
        outStaveEvent.ch_ev_buffer = outChipEvents;
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

      outDetEvent.st_ev_buffer = outStaveEvents;
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

    localStats.print();
    localStats.reset();
    std::cout << "Finished " << runPrefix << "\n";
    std::cout << "\n\n\n\n";
  }

  std::cout << "Finished processing\n";
  globalStats.print();
}

std::pair<AcceleratorState, AcceleratorState>
Preprocessor::getAcceleratorStateMeanStd(TTree *tree, std::size_t skip,
                                         const Options &opt) const {
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
    bool isShuttered =
        static_cast<int>(std::round((epicsPID - 10) / 36)) % 20 == 0;
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
      // TORO:LI20:3163:TMIT_PC
      double toro_li20_3163_tmit_pc_val =
          getPv(detEvent, toro_li20_3163_tmit_pc_key);
      toro_li20_3163_tmit_pc_hist.Fill(toro_li20_3163_tmit_pc_val);
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
  double quantilePos[2] = {opt.noiseQuantileLowAcc, opt.noiseQuantileHighAcc};

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
    double toro_li20_3163_tmit_pc_val = 0;
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
      toro_li20_3163_tmit_pc_val = getPv(detEvent, toro_li20_3163_tmit_pc_key);
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
    if (toro_li20_3163_tmit_pc_val < accStateMin.toro_li20_3163_tmit_pc ||
        toro_li20_3163_tmit_pc_val > accStateMax.toro_li20_3163_tmit_pc) {
      continue;
    }
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
    // TORO:LI20:3163:TMIT_PC
    toro_li20_3163_tmit_pc_hist.Fill(toro_li20_3163_tmit_pc_val);
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

std::unordered_set<std::size_t> Preprocessor::getBadEvents(
    TTree *tree, std::size_t skip, const Options &opt) const {
  std::cout << "Getting bad events with AccStdDev = " << opt.nStdDevsAcc
            << "\n";
  // Get accelerator state parameters
  auto [accStateMean, accStateStdDev] =
      getAcceleratorStateMeanStd(tree, skip, opt);

  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);

  std::unordered_set<std::size_t> badEvents;
  auto isBadCentral = [&opt](double mean, double std, double val) {
    return (std::abs(mean - val) / std > opt.nStdDevsAcc);
  };
  auto isBadUpper = [&opt](double mean, double std, double val) {
    return ((val - mean) / std > opt.nStdDevsAcc);
  };

  for (std::size_t i = skip; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    std::size_t epicsPID = 0;
    try {
      epicsPID = getPv(detEvent, epics_pid_key);
    } catch (std::out_of_range &err) {
      badEvents.insert(i);
      continue;
    }
    std::size_t pid =
        static_cast<std::size_t>(std::round((epicsPID - 10) / 36));

    double pmt_s20_3060_val = 0;
    double pmt_s20_3070_val = 0;
    double pmt_s20_3179_val = 0;
    double pmt_s20_3350_val = 0;
    double pmt_s20_3360_val = 0;

    double toro_li20_1988_tmit_pc_val = 0;
    double toro_li20_2040_tmit_pc_val = 0;
    double toro_li20_2452_tmit_pc_val = 0;
    double toro_li20_3163_tmit_pc_val = 0;
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
      toro_li20_3163_tmit_pc_val = getPv(detEvent, toro_li20_3163_tmit_pc_key);
      toro_li20_3255_tmit_pc_val = getPv(detEvent, toro_li20_3255_tmit_pc_key);
    } catch (std::out_of_range &err) {
      std::cout << "MISSING PVs at " << detEvent->trg_n << "\n";
      badEvents.insert(i);
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
    if (pmt_s20_3070_val > 1000) {
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
    if (isBadCentral(accStateMean.toro_li20_1988_tmit_pc,
                     accStateStdDev.toro_li20_1988_tmit_pc,
                     toro_li20_1988_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2040:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_2040_tmit_pc,
                     accStateStdDev.toro_li20_2040_tmit_pc,
                     toro_li20_2040_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:2452:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_2452_tmit_pc,
                     accStateStdDev.toro_li20_2452_tmit_pc,
                     toro_li20_2452_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:3163:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_3163_tmit_pc,
                     accStateStdDev.toro_li20_3163_tmit_pc,
                     toro_li20_3163_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
    // TORO:LI20:3255:TMIT_PC
    if (isBadCentral(accStateMean.toro_li20_3255_tmit_pc,
                     accStateStdDev.toro_li20_3255_tmit_pc,
                     toro_li20_3255_tmit_pc_val)) {
      badEvents.insert(i);
      continue;
    }
  }
  std::cout << "Found bad events " << badEvents.size() << "\n";
  std::cout << "\n\n\n\n";
  return badEvents;
}

void Preprocessor::bfsClustering(std::vector<Hit> &clusterPixels, Hit &pivot,
                                 HitSet &pixels) const {
  std::deque<Hit> queue{pivot};
  while (!queue.empty()) {
    auto current = queue.front();
    queue.pop_front();
    auto currentPtr = pixels.find(current);

    if (currentPtr != pixels.end()) {
      // Add this pixel to the cluster
      clusterPixels.push_back(current);
      // Remove pixel from live pixels list
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

std::vector<Preprocessor::Cluster> Preprocessor::getClusters(
    HitSet &hits, std::size_t maxCluterSize) const {
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

double Preprocessor::convertEpicsDipole(double val, double theta,
                                        double bz) const {
  const double c_light = 2.99792458E8;  // speed of light in m/s
  const double me = 0.5109989461e6;     // electron mass in eV

  return std::sqrt(val * val - me * me) * std::sin(theta) / (c_light * bz);
}

double Preprocessor::convertEpicsXCorrector(double val, double bz) const {
  return val / bz * 0.1;
}

double Preprocessor::convertEpicsQuad(double val) const {
  return val * 0.1;
}

double Preprocessor::getPv(const detector_event *detEvent,
                           const std::string &pvname) const {
  try {
    return std::stod(detEvent->ev_epics_frame.pv_map.at(pvname).value);
  } catch (std::out_of_range &err) {
    std::cout << "MISSING PV " << pvname << "\n";
    throw err;
  }
}

void Preprocessor::Stats::init(const std::string &prefix,
                               const std::set<std::size_t> &ids) {
  runPrefix = prefix;
  chipIds = ids;

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

  // Raw PMTs
  gr_PMT_S20_3060_Raw = TGraph();
  gr_PMT_S20_3070_Raw = TGraph();
  gr_PMT_S20_3179_Raw = TGraph();
  gr_PMT_S20_3350_Raw = TGraph();
  gr_PMT_S20_3360_Raw = TGraph();
  // Filtered PMTs
  gr_PMT_S20_3060_Filtered = TGraph();
  gr_PMT_S20_3070_Filtered = TGraph();
  gr_PMT_S20_3179_Filtered = TGraph();
  gr_PMT_S20_3350_Filtered = TGraph();
  gr_PMT_S20_3360_Filtered = TGraph();
  // Raw toroids
  gr_TORO_LI20_1988_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_2040_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_2452_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_3163_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_3255_TMIT_PC_Raw = TGraph();
  // Filtered toroids
  gr_TORO_LI20_1988_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_2040_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_2452_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_3163_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_3255_TMIT_PC_Filtered = TGraph();
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

  // PMT_S20_3060
  GraphPainterConfig gr_PMT_S20_3060_Raw_Cfg;
  gr_PMT_S20_3060_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3060_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3060_Raw.SetNameTitle("gr_PMT_S20_3060", "gr_PMT_S20_3060");
  GraphPainterConfig gr_PMT_S20_3060_Filtered_Cfg;
  gr_PMT_S20_3060_Filtered_Cfg.option = "P";
  gr_PMT_S20_3060_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3060_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3060_Raw, gr_PMT_S20_3060_Raw_Cfg},
             {gr_PMT_S20_3060_Filtered, gr_PMT_S20_3060_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3070
  GraphPainterConfig gr_PMT_S20_3070_Raw_Cfg;
  gr_PMT_S20_3070_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3070_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3070_Raw.SetNameTitle("gr_PMT_S20_3070", "gr_PMT_S20_3070");
  GraphPainterConfig gr_PMT_S20_3070_Filtered_Cfg;
  gr_PMT_S20_3070_Filtered_Cfg.option = "P";
  gr_PMT_S20_3070_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3070_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3070_Raw, gr_PMT_S20_3070_Raw_Cfg},
             {gr_PMT_S20_3070_Filtered, gr_PMT_S20_3070_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3179
  GraphPainterConfig gr_PMT_S20_3179_Raw_Cfg;
  gr_PMT_S20_3179_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3179_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3179_Raw.SetNameTitle("gr_PMT_S20_3179", "gr_PMT_S20_3179");
  GraphPainterConfig gr_PMT_S20_3179_Filtered_Cfg;
  gr_PMT_S20_3179_Filtered_Cfg.option = "P";
  gr_PMT_S20_3179_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3179_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3179_Raw, gr_PMT_S20_3179_Raw_Cfg},
             {gr_PMT_S20_3179_Filtered, gr_PMT_S20_3179_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3350
  GraphPainterConfig gr_PMT_S20_3350_Raw_Cfg;
  gr_PMT_S20_3350_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3350_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3350_Raw.SetNameTitle("gr_PMT_S20_3350", "gr_PMT_S20_3350");
  GraphPainterConfig gr_PMT_S20_3350_Filtered_Cfg;
  gr_PMT_S20_3350_Filtered_Cfg.option = "P";
  gr_PMT_S20_3350_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3350_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3350_Raw, gr_PMT_S20_3350_Raw_Cfg},
             {gr_PMT_S20_3350_Filtered, gr_PMT_S20_3350_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3360
  GraphPainterConfig gr_PMT_S20_3360_Raw_Cfg;
  gr_PMT_S20_3360_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3360_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3360_Raw.SetNameTitle("gr_PMT_S20_3360", "gr_PMT_S20_3360");
  GraphPainterConfig gr_PMT_S20_3360_Filtered_Cfg;
  gr_PMT_S20_3360_Filtered_Cfg.option = "P";
  gr_PMT_S20_3360_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3360_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3360_Raw, gr_PMT_S20_3360_Raw_Cfg},
             {gr_PMT_S20_3360_Filtered, gr_PMT_S20_3360_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_1988_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_1988_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_1988_TMIT_PC",
                                             "gr_TORO_LI20_1988_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_1988_TMIT_PC_Raw, gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_1988_TMIT_PC_Filtered,
              gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_2040_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_2040_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_2040_TMIT_PC",
                                             "gr_TORO_LI20_2040_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_2040_TMIT_PC_Raw, gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2040_TMIT_PC_Filtered,
              gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_2452_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_2452_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_2452_TMIT_PC",
                                             "gr_TORO_LI20_2452_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_2452_TMIT_PC_Raw, gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2452_TMIT_PC_Filtered,
              gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_3163_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_3163_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_3163_TMIT_PC",
                                             "gr_TORO_LI20_3163_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_3163_TMIT_PC_Raw, gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_3163_TMIT_PC_Filtered,
              gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_3255_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_3255_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_3255_TMIT_PC",
                                             "gr_TORO_LI20_3255_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_3255_TMIT_PC_Raw, gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_3255_TMIT_PC_Filtered,
              gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg}},
            legCfg);

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
