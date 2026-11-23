#include "Analysis/PreprocessorSim.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <RtypesCore.h>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"

PreprocessorSim::FileHandle PreprocessorSim::getFileHandle(
    const std::string& treeName, const std::vector<std::string>& paths) const {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  TChain* chainOwner = new TChain(treeName.c_str(), treeName.c_str());
  if (paths.size() == 1) {
    file = new TFile(paths.at(0).c_str());
    tree = file->Get<TTree>(treeName.c_str());
  } else {
    chainOwner = new TChain(treeName.c_str());
    for (const auto& path : paths) {
      chainOwner->Add(path.c_str());
    }
    tree = dynamic_cast<TTree*>(chainOwner);
  }

  return {file, chainOwner, tree};
}

void PreprocessorSim::processFiles(const Options& opt) const {
  // --------------------------------------------------
  // Get tracker data tree
  std::vector<std::string> inTrackerPaths;
  for (const auto& entry :
       std::filesystem::directory_iterator(opt.inTrackerDir)) {
    inTrackerPaths.push_back(entry.path().string());
  }
  auto [inTrackerFile, inTrackerChain, inTrackerTree] =
      getFileHandle(opt.inTrackerTreeName, inTrackerPaths);

  // Get beam data tree
  std::vector<std::string> inBeamPaths;
  for (const auto& entry : std::filesystem::directory_iterator(opt.inBeamDir)) {
    inBeamPaths.push_back(entry.path().string());
  }
  auto [inBeamFile, inBeamChain, inBeamTree] =
      getFileHandle(opt.inBeamTreeName, inBeamPaths);

  // --------------------------------------------------
  // Output files
  TFile* outFile = new TFile(opt.outPath.c_str(), "RECREATE");
  TTree* outDataTree = nullptr;

  // --------------------------------------------------
  // Initialize tracker branches

  // Cluster hit in the surface frame
  TVector2* geoCenterLocal = nullptr;

  // Cluster hit in the global frame
  TVector3* geoCenterGlobal = nullptr;

  // Cluster hit covariance in the surface frame
  TMatrixD* clusterCov = nullptr;

  // Momentum direction measurement in the track frame
  TVector3* onSurfaceDirection = nullptr;

  // Angle covariance in the surface frame
  TMatrixD* angleCov = nullptr;

  // Surface geometry ID
  std::size_t geoId = 0;

  // Event ID
  std::size_t trackerEventId = 0;

  // Signal/background flag
  int isSignal = 0;

  // True track hits within the clusters
  std::vector<TVector2>* trackHitsLocal = nullptr;
  std::vector<TVector3>* trackHitsGlobal = nullptr;

  // Track IDs
  std::vector<int>* trackId = nullptr;
  std::vector<int>* parentTrackId = nullptr;
  std::vector<int>* runId = nullptr;

  // Bound origin parameters
  std::vector<TVectorD>* boundTrackParameters = nullptr;
  std::vector<TMatrixD>* boundTrackCov = nullptr;

  // Origin momentum
  std::vector<TLorentzVector>* originMomentum = nullptr;

  // Origin vertex
  std::vector<TVector3>* vertex = nullptr;

  // Momentum at clusters
  std::vector<TLorentzVector>* onSurfaceMomentumTruth = nullptr;

  // Charges of the tracks' particles
  std::vector<int>* charge = nullptr;

  // PDG IDs of the tracks' particles
  std::vector<int>* pdgId = nullptr;

  // Cluster hit in the surface frame
  inTrackerTree->SetBranchAddress("geoCenterLocal", &geoCenterLocal);

  // Cluster hit in the global frame
  inTrackerTree->SetBranchAddress("geoCenterGlobal", &geoCenterGlobal);

  // Cluster hit covariance in the surface frame
  inTrackerTree->SetBranchAddress("clusterCov", &clusterCov);

  // Momentum direction measurement in the track frame
  inTrackerTree->SetBranchAddress("onSurfaceDirection", &onSurfaceDirection);

  // Track angular covariance in the surface frame
  inTrackerTree->SetBranchAddress("angleCov", &angleCov);

  // Surface geometry ID
  inTrackerTree->SetBranchAddress("geoId", &geoId);

  // Signal/background flag
  inTrackerTree->SetBranchAddress("isSignal", &isSignal);

  // True track hits within the clusters
  inTrackerTree->SetBranchAddress("trackHitsGlobal", &trackHitsGlobal);
  inTrackerTree->SetBranchAddress("trackHitsLocal", &trackHitsLocal);

  // Track IDs
  inTrackerTree->SetBranchAddress("trackId", &trackId);
  inTrackerTree->SetBranchAddress("parentTrackId", &parentTrackId);
  inTrackerTree->SetBranchAddress("runId", &runId);

  // Bound origin parameters
  inTrackerTree->SetBranchAddress("boundTrackParameters",
                                  &boundTrackParameters);
  inTrackerTree->SetBranchAddress("boundTrackCov", &boundTrackCov);

  // Origin momentum
  inTrackerTree->SetBranchAddress("originMomentum", &originMomentum);

  // Origin vertex
  inTrackerTree->SetBranchAddress("vertex", &vertex);

  // Momentum at clusters
  inTrackerTree->SetBranchAddress("onSurfaceMomentumTruth",
                                  &onSurfaceMomentumTruth);

  // Charges of the tracks' particles
  inTrackerTree->SetBranchAddress("charge", &charge);

  // PDG IDs of the tracks' particles
  inTrackerTree->SetBranchAddress("pdgId", &pdgId);

  // Event ID
  inTrackerTree->SetBranchAddress("eventId", &trackerEventId);

  // --------------------------------------------------
  // Initialize beam branches

  // Reconstructed beam position at the reference
  TVector3* referencePositionEst = nullptr;

  // Reconstructed beam momentum at the reference
  TLorentzVector* referenceMomentumEst = nullptr;

  // Covariance of the orbit fit at the reference
  TMatrixD* referenceCov = nullptr;

  // Covariance of the beam bunch parameters at the reference
  TMatrixD* referenceBeamBunchCov = nullptr;

  // Covariance at the reference due to the alignment uncertainty
  TMatrixD* referenceAlignmentCov = nullptr;

  // Covariance at the reference of the gradient dependent systematic
  TMatrixD* referenceGradientSystematicCov = nullptr;

  // Shifts at the reference for the alignment samples
  std::vector<TVectorD>* referenceAlignmentSampleShifts = nullptr;

  // Covariance at the reference from the alignment samples
  TMatrixD* referenceAlignmentSampleCov = nullptr;

  // chi2 of the orbit fit
  double referenceFitChi2Est = 0;

  // True beam position at the reference
  TVector3* referencePositionTruth = nullptr;

  // True beam momentum at the reference
  TLorentzVector* referenceMomentumTruth = nullptr;

  // Event ID
  ULong64_t beamEventId = 0;

  inBeamTree->SetBranchAddress("referencePositionEst", &referencePositionEst);
  inBeamTree->SetBranchAddress("referenceMomentumEst", &referenceMomentumEst);
  inBeamTree->SetBranchAddress("referenceCov", &referenceCov);
  inBeamTree->SetBranchAddress("referenceBeamBunchCov", &referenceBeamBunchCov);
  inBeamTree->SetBranchAddress("referenceAlignmentCov", &referenceAlignmentCov);
  inBeamTree->SetBranchAddress("referenceGradientSystematicCov",
                               &referenceGradientSystematicCov);
  inBeamTree->SetBranchAddress("referenceAlignmentSampleShifts",
                               &referenceAlignmentSampleShifts);
  inBeamTree->SetBranchAddress("referenceAlignmentSampleCov",
                               &referenceAlignmentSampleCov);
  inBeamTree->SetBranchAddress("referenceFitChi2Est", &referenceFitChi2Est);
  inBeamTree->SetBranchAddress("referencePositionTruth",
                               &referencePositionTruth);
  inBeamTree->SetBranchAddress("referenceMomentumTruth",
                               &referenceMomentumTruth);
  inBeamTree->SetBranchAddress("eventId", &beamEventId);

  // Map the event IDs to the orbit entries
  std::map<std::size_t, std::size_t> beamEntries;
  for (std::size_t i = 0; i < inBeamTree->GetEntries(); i++) {
    inBeamTree->GetEntry(i);
    auto [it, inserted] = beamEntries.insert({beamEventId, i});
    if (!inserted) {
      throw std::runtime_error("Duplicate orbit reconstruction for event " +
                               std::to_string(beamEventId));
    }
  }

  // --------------------------------------------------
  // Output tree with the tracker branches
  outDataTree = inTrackerTree->CloneTree(0);
  outDataTree->CopyAddresses(inTrackerTree);

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
  TVector3 outOrbitReferencePositionTruth;
  TLorentzVector outOrbitReferenceMomentumTruth;

  std::size_t bufSize = 32000;
  std::size_t splitLvl = 0;
  outDataTree->Branch("orbitReferencePositionEst",
                      &outOrbitReferencePositionEst, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceMomentumEst",
                      &outOrbitReferenceMomentumEst, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceCov", &outOrbitReferenceCov, bufSize,
                      splitLvl);
  outDataTree->Branch("orbitReferenceBeamBunchCov",
                      &outOrbitReferenceBeamBunchCov, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceAlignmentCov",
                      &outOrbitReferenceAlignmentCov, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceGradientSystematicCov",
                      &outOrbitReferenceGradientSystematicCov, bufSize,
                      splitLvl);
  outDataTree->Branch("orbitReferenceAlignmentSampleShifts",
                      &outOrbitReferenceAlignmentSampleShifts, bufSize,
                      splitLvl);
  outDataTree->Branch("orbitReferenceAlignmentSampleCov",
                      &outOrbitReferenceAlignmentSampleCov, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceFitChi2Est", &outOrbitReferenceFitChi2Est,
                      bufSize, splitLvl);
  outDataTree->Branch("orbitReferencePositionTruth",
                      &outOrbitReferencePositionTruth, bufSize, splitLvl);
  outDataTree->Branch("orbitReferenceMomentumTruth",
                      &outOrbitReferenceMomentumTruth, bufSize, splitLvl);

  // --------------------------------------------------
  // Fill in the data
  std::cout << "Filling in data\n";
  std::size_t nEntries = inTrackerTree->GetEntries();
  std::size_t nRowsNoOrbit = 0;
  for (std::size_t i = 0; i < nEntries; i++) {
    inTrackerTree->GetEntry(i);

    // Continue if the event has no orbit reconstruction
    auto beamIt = beamEntries.find(trackerEventId);
    if (beamIt == beamEntries.end()) {
      nRowsNoOrbit++;
      continue;
    }
    inBeamTree->GetEntry(beamIt->second);

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
    outOrbitReferenceAlignmentSampleCov.ResizeTo(*referenceAlignmentSampleCov);
    outOrbitReferenceAlignmentSampleCov = *referenceAlignmentSampleCov;
    outOrbitReferenceFitChi2Est = referenceFitChi2Est;
    outOrbitReferencePositionTruth = *referencePositionTruth;
    outOrbitReferenceMomentumTruth = *referenceMomentumTruth;

    // The orbit at the reference is the IP measurement
    if (geoId == 8) {
      geoCenterLocal->Set(referencePositionEst->Y(),
                          -referencePositionEst->Z());
      geoCenterGlobal->SetXYZ(referencePositionEst->X(),
                              referencePositionEst->Y(),
                              referencePositionEst->Z());
      onSurfaceDirection->SetMagThetaPhi(1.0, referenceMomentumEst->Theta(),
                                         referenceMomentumEst->Phi());

      // Position and angle blocks of the orbit fit covariance,
      // (loc0, loc1, phi, theta) order
      for (std::size_t n = 0; n < clusterCov->GetNrows(); n++) {
        for (std::size_t m = 0; m < clusterCov->GetNcols(); m++) {
          clusterCov->operator()(n, m) = referenceCov->operator()(n, m);
        }
      }
      for (std::size_t n = 0; n < angleCov->GetNrows(); n++) {
        for (std::size_t m = 0; m < angleCov->GetNcols(); m++) {
          angleCov->operator()(n, m) = referenceCov->operator()(n + 2, m + 2);
        }
      }
    }

    outDataTree->Fill();
    if (i % 100 == 0) {
      std::cout << i << "/" << nEntries << "\n";
    }
  }
  std::cout << "Rows of events without orbit reconstruction: " << nRowsNoOrbit
            << "\n";

  outDataTree->Write();
  outFile->Close();
  std::cout << "\n\n\n\n";
  std::cout << "Finished processing\n";
}
