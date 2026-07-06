#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <csignal>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "Analysis/AnalysisFuncitons.hpp"
#include "DetectorEvent.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "Types/DetectorEvent.hpp"
#include "detail/PairHash.hpp"

struct PathCollection {
  std::string runPrefix = "";
  std::vector<std::string> inPaths{};
  std::string outPath = "";
  std::string magnetPath = "";
};

int main() {
  // Directories with run data
  std::vector<std::string> inDirs = {
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_Feb_2025/raw/data_Run503"};

  // Out path prefix
  std::string outPathPrefix =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_Feb_2025/"
      "processed";

  // Collect paths
  std::unordered_map<std::string, PathCollection> paths;
  for (const auto &inPath : inDirs) {
    std::size_t runIdx = inPath.find_last_of("/") + 1;
    std::string runString = inPath.substr(runIdx);
    paths[runString].runPrefix = runString;

    paths[runString].outPath = outPathPrefix + "/" + runString + ".root";
    paths[runString].magnetPath =
        outPathPrefix + "/" + runString + "_Magnets.root";
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(inPath)) {
      if (!entry.is_regular_file() || entry.path().extension() != ".root") {
        continue;
      }
      paths[runString].inPaths.push_back(entry.path().string());
    }
  }
  for (const auto &[runPrefix, collection] : paths) {
    std::cout << "------------------------------\n";
    std::cout << collection.runPrefix << "\n";
    std::cout << collection.outPath << "\n";
    std::cout << collection.magnetPath << "\n";
    for (const auto &inPath : collection.inPaths) {
      std::cout << inPath << "\n";
    }
  }

  // Readout parameters and cuts
  std::string treeName = "MyTree";

  int nStdDevsAcc = 3;
  int nStdDevsPix = 3;
  std::size_t skip = 0;
  std::size_t stop = std::numeric_limits<std::size_t>::max();

  std::size_t maxClusterSize = 4;

  int requiredStaveSize = 1;
  int requiredChipSize = 5;

  std::unordered_set<int> chipIds = {0, 2, 4, 6, 8};

  // Iterate over the run paths
  for (const auto &[runPrefix, collection] : paths) {
    // Add all the run files to a tree
    const auto &inPaths = collection.inPaths;
    TFile *file = nullptr;
    TTree *tree = nullptr;
    TChain *chainOwner = new TChain(treeName.c_str(), treeName.c_str());
    if (inPaths.size() == 1) {
      file = new TFile(inPaths.at(0).c_str());
      tree = file->Get<TTree>(treeName.c_str());
    } else {
      chainOwner = new TChain(treeName.c_str());
      for (const auto &path : inPaths) {
        chainOwner->Add(path.c_str());
      }
      tree = dynamic_cast<TTree *>(chainOwner);
    }

    // Output files
    TFile *outFile = new TFile(collection.outPath.c_str(), "RECREATE");
    TTree *outTree = new TTree("MyTree", "MyTree");

    // Get accelerator state parameters
    auto [accStateMean, accStateStdDev] =
        getAcceleratorStateMeanStd(tree, skip);

    // Collect bad events
    std::unordered_set<std::size_t> badEvents =
        getBadEvents(tree, accStateMean, accStateStdDev, nStdDevsAcc, skip);

    // Initialize tree readout
    E320Io::DetectorEvent outDetEvent;
    ULong64_t eventId = 0;
    outTree->Branch("event", &outDetEvent);
    outTree->Branch("eventId", &eventId);

    detector_event *detEvent = nullptr;
    tree->SetBranchAddress("event", &detEvent);
    std::size_t nEntries =
        std::min(static_cast<std::size_t>(tree->GetEntries()), stop);

    // Count noisy pixels
    std::unordered_map<int, std::unordered_map<Hit, double, PairHash>>
        totHitCountsPix;
    for (std::size_t i = skip; i < nEntries; i++) {
      if (badEvents.contains(i)) {
        continue;
      }
      tree->GetEntry(i);
      for (const auto &staveEv : detEvent->st_ev_buffer) {
        for (const auto &chipEv : staveEv.ch_ev_buffer) {
          int id = 100 * detEvent->run_number + 10 * staveEv.stave_id +
                   chipEv.chip_id;
          for (const auto &[hitX, hitY] : chipEv.hits) {
            totHitCountsPix[id][{hitX, hitY}]++;
          }
        }
      }
    }
    std::cout << "TOT HIT COUNTS " << totHitCountsPix.size() << "\n";

    std::unordered_map<int, double> meanHitCountsPix;
    std::unordered_map<int, double> mean2HitCountsPix;
    for (const auto &[id, hitCounts] : totHitCountsPix) {
      double norm = hitCounts.size();
      for (const auto &[hit, count] : hitCounts) {
        meanHitCountsPix[id] += count / norm;
        mean2HitCountsPix[id] += count * count / norm;
      }
    }
    std::cout << "MEAN HIT COUNTS " << meanHitCountsPix.size() << "\n";

    std::unordered_map<int, std::unordered_set<Hit, PairHash>> noisyPixels;
    for (const auto &[id, hitCounts] : totHitCountsPix) {
      double mean = meanHitCountsPix[id];
      double mean2 = mean2HitCountsPix[id];
      double std = std::sqrt(mean2 - mean * mean);

      for (const auto &[hit, count] : hitCounts) {
        if (std::abs(count - mean) / std > nStdDevsPix) {
          noisyPixels[id].insert(hit);
        }
      }
    }
    std::cout << "NOISY PIXELS " << noisyPixels.size() << "\n";

    TH1I *rawOcc = new TH1I("rawOcc", "", 1024 * 512, 0, 1024 * 512);
    TH1I *filteredOcc = new TH1I("filteredOcc", "", 1024 * 512, 0, 1024 * 512);

    std::map<int, TH2D *> chipOccRawHists;
    std::map<int, TH2D *> chipOccFilteredHists;
    for (int id : chipIds) {
      std::string nameRaw = "chipOcc" + std::to_string(id) + "Raw";
      std::string nameFiltered = "chipOcc" + std::to_string(id) + "Filtered";
      chipOccRawHists.insert({id, new TH2D(nameRaw.c_str(), nameRaw.c_str(),
                                           200, 0, 1024, 200, 0, 512)});
      chipOccFilteredHists.insert(
          {id, new TH2D(nameFiltered.c_str(), nameFiltered.c_str(), 200, 0,
                        1024, 200, 0, 512)});
    }

    std::ofstream bpmCsvX("bpms_503_X.csv", std::ofstream::out);
    bpmCsvX << "[TRG N],[3156 X],[3218 X],[3265 X],[3315 X],[Q0 BACT],[Q1 "
               "BACT],[Q2 BACT]\n";

    std::ofstream bpmCsvY("bpms_503_Y.csv", std::ofstream::out);
    bpmCsvY << "[TRG N],[3156 Y],[3218 Y],[3265 Y],[3315 X],[Q0 BACT],[Q1 "
               "BACT],[Q2 BACT]\n";

    // Fill in the data
    for (std::size_t i = skip; i < nEntries; i++) {
      // if (badEvents.contains(i)) {
      //   continue;
      // }
      tree->GetEntry(i);

      outDetEvent.trg_n = detEvent->trg_n;
      outDetEvent.ts_begin = detEvent->ts_begin;
      outDetEvent.ts_end = detEvent->ts_end;
      outDetEvent.run_number = detEvent->run_number;

      std::vector<E320Io::StaveEvent> outStaveEvents;
      outStaveEvents.reserve(detEvent->st_ev_buffer.size());

      E320Io::BPMEvent bpm3156Ev;
      bpm3156Ev.meanX = detEvent->epics_frame.bpm_pb_3156_x;
      bpm3156Ev.meanY = detEvent->epics_frame.bpm_pb_3156_y;
      bpm3156Ev.stdDevX = 0.1;
      bpm3156Ev.stdDevY = 0.1;
      bpm3156Ev.id = 3156;

      E320Io::BPMEvent bpm3218Ev;
      bpm3218Ev.meanX = detEvent->epics_frame.bpm_quad0_3218_x;
      bpm3218Ev.meanY = detEvent->epics_frame.bpm_quad0_3218_y;
      bpm3218Ev.stdDevX = 0.1;
      bpm3218Ev.stdDevY = 0.1;
      bpm3218Ev.id = 3218;

      E320Io::BPMEvent bpm3265Ev;
      bpm3265Ev.meanX = detEvent->epics_frame.bpm_quad1_3265_x;
      bpm3265Ev.meanY = detEvent->epics_frame.bpm_quad1_3265_y;
      bpm3265Ev.stdDevX = 0.1;
      bpm3265Ev.stdDevY = 0.1;
      bpm3265Ev.id = 3265;

      E320Io::BPMEvent bpm3315Ev;
      bpm3315Ev.meanX = detEvent->epics_frame.bpm_quad2_3315_x;
      bpm3315Ev.meanY = detEvent->epics_frame.bpm_quad2_3315_y;
      bpm3315Ev.stdDevX = 0.1;
      bpm3315Ev.stdDevY = 0.1;
      bpm3315Ev.id = 3315;

      bpmCsvX << detEvent->trg_n << "," << detEvent->epics_frame.bpm_pb_3156_x
              << "," << detEvent->epics_frame.bpm_quad0_3218_x << ","
              << detEvent->epics_frame.bpm_quad1_3265_x << ","
              << detEvent->epics_frame.bpm_quad2_3315_x << ","
              << detEvent->epics_frame.espec_quad0_bact << ","
              << detEvent->epics_frame.espec_quad1_bact << ","
              << detEvent->epics_frame.espec_quad2_bact << "\n";

      bpmCsvY << detEvent->trg_n << "," << detEvent->epics_frame.bpm_pb_3156_y
              << "," << detEvent->epics_frame.bpm_quad0_3218_y << ","
              << detEvent->epics_frame.bpm_quad1_3265_y << ","
              << detEvent->epics_frame.bpm_quad2_3315_y << ","
              << detEvent->epics_frame.espec_quad0_bact << ","
              << detEvent->epics_frame.espec_quad1_bact << ","
              << detEvent->epics_frame.espec_quad2_bact << "\n";

      outDetEvent.bpm_ev_buffer = {bpm3156Ev, bpm3218Ev, bpm3265Ev, bpm3315Ev};

      for (const auto &staveEv : detEvent->st_ev_buffer) {
        if (staveEv.ch_ev_buffer.size() != requiredChipSize) {
          continue;
        }
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

          std::unordered_set<Hit, PairHash> hits;
          for (const auto &hit : chipEv.hits) {
            rawOcc->Fill(512 * hit.first + hit.second);
            chipOccRawHists.at(chipEv.chip_id)->Fill(hit.first, hit.second);
            if (noisyPixels[id].contains(hit)) {
              continue;
            }
            hits.insert(hit);
            chipOccFilteredHists.at(chipEv.chip_id)
                ->Fill(hit.first, hit.second);
            filteredOcc->Fill(512 * hit.first + hit.second);
          }
          if (hits.empty()) {
            continue;
          }

          outChipEvent.hits = getClusters(hits, maxClusterSize);
          if (outChipEvent.hits.empty()) {
            continue;
          }
          outChipEvents.push_back(std::move(outChipEvent));
        }
        if (outChipEvents.size() != requiredChipSize) {
          continue;
        }
        outStaveEvent.ch_ev_buffer = outChipEvents;
        outStaveEvents.push_back(std::move(outStaveEvent));
      }
      if (outStaveEvents.empty()) {
        continue;
      }
      outDetEvent.st_ev_buffer = outStaveEvents;
      outTree->Fill();

      eventId++;
      if (i % 100 == 0) {
        std::cout << i << "/" << nEntries << "\n";
      }
    }
    bpmCsvX.close();
    bpmCsvY.close();

    saveGraphsHists(tree, skip, stop, badEvents);

    double padLeftMargin = 0.15;
    double padRightMargin = 0.1;
    double padBottomMargin = 0.1;
    double padTopMargin = 0.1;

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetMargin(padLeftMargin, padRightMargin, padBottomMargin, padTopMargin);
    rawOcc->SetStats(0);
    rawOcc->Draw("hist");
    rawOcc->GetXaxis()->SetTitle("Pix. Id");
    rawOcc->GetYaxis()->SetTitle("Count");
    c1->SaveAs("raw.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    c2->SetMargin(padLeftMargin, padRightMargin, padBottomMargin, padTopMargin);
    filteredOcc->SetStats(0);
    filteredOcc->Draw("hist");
    filteredOcc->GetXaxis()->SetTitle("Pix. Id");
    filteredOcc->GetYaxis()->SetTitle("Count");
    c2->SaveAs("filtered.pdf");

    for (const auto &[id, hist] : chipOccRawHists) {
      std::string cName = "cRaw" + std::to_string(id);
      TCanvas *c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 500);
      hist->Draw("colz");
      hist->GetXaxis()->SetTitle("Pix. Id");
      hist->GetYaxis()->SetTitle("Count");

      std::string fileName = "chipOccRaw" + std::to_string(id) + ".pdf";
      c->SaveAs(fileName.c_str());
    }

    for (const auto &[id, hist] : chipOccFilteredHists) {
      std::string cName = "cFiltered" + std::to_string(id);
      TCanvas *c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 500);
      hist->Draw("colz");
      hist->GetXaxis()->SetTitle("Pix. Id");
      hist->GetYaxis()->SetTitle("Count");

      std::string fileName = "chipOccFiltered" + std::to_string(id) + ".pdf";
      c->SaveAs(fileName.c_str());
    }

    outTree->Write();
    outFile->Close();
  }

  return 0;
}
