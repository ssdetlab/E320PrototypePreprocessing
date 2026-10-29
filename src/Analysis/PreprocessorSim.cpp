#include "Analysis/PreprocessorSim.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
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

  // Reconstructed beam position
  TVector3* referenceBeamPosition = nullptr;

  // Reconstructed beam momentum
  TLorentzVector* referenceBeamMomentum = nullptr;

  // Reconstructed beam position covariance
  TMatrixD* referenceBeamPositionCov = nullptr;

  // Reconstructed beam momentum covariance
  TMatrixD* referenceBeamAngleCov = nullptr;

  // Event ID
  ULong64_t beamEventId = 0;

  // Reconstructed beam position
  inBeamTree->SetBranchAddress("referencePositionEst", &referenceBeamPosition);

  // Reconstructed beam momentum
  inBeamTree->SetBranchAddress("referenceMomentumEst", &referenceBeamMomentum);

  // Reconstructed beam position covariance
  inBeamTree->SetBranchAddress("referencePositionCov", &referenceBeamPositionCov);

  // Reconstructed beam momentum covariance
  inBeamTree->SetBranchAddress("referenceAngleCov", &referenceBeamAngleCov);

  // Event ID
  inBeamTree->SetBranchAddress("eventId", &beamEventId);

  // Fill in the data
  outDataTree = inTrackerTree->CloneTree(0);
  outDataTree->CopyAddresses(inTrackerTree);

  std::cout << "Filling in data\n";
  std::size_t nEntries = inTrackerTree->GetEntries();
  std::size_t currentBeamIdx = 0;
  for (std::size_t i = 0; i < nEntries; i++) {
    inTrackerTree->GetEntry(i);

    if (geoId == 8) {
      inBeamTree->GetEntry(currentBeamIdx);
      if (trackerEventId < beamEventId) {
        continue;
      }
      while (trackerEventId != beamEventId) {
        currentBeamIdx++;
        inBeamTree->GetEntry(currentBeamIdx);
      }

      geoCenterLocal->Set(referenceBeamPosition->Y(),
                          -referenceBeamPosition->Z());
      geoCenterGlobal->SetXYZ(referenceBeamPosition->X(),
                              referenceBeamPosition->Y(),
                              referenceBeamPosition->Z());
      onSurfaceDirection->SetMagThetaPhi(1.0, referenceBeamMomentum->Theta(),
                                         referenceBeamMomentum->Phi());

      for (std::size_t n = 0; n < clusterCov->GetNrows(); n++) {
        for (std::size_t m = 0; m < clusterCov->GetNcols(); m++) {
          clusterCov->operator()(n, m) = referenceBeamPositionCov->operator()(n, m);
        }
      }
      for (std::size_t n = 0; n < angleCov->GetNrows(); n++) {
        for (std::size_t m = 0; m < angleCov->GetNcols(); m++) {
          angleCov->operator()(n, m) = referenceBeamAngleCov->operator()(n, m);
        }
      }
      currentBeamIdx++;
    }

    outDataTree->Fill();
    if (i % 100 == 0) {
      std::cout << i << "/" << nEntries << "\n";
    }
  }

  outDataTree->Write();
  outFile->Close();
  std::cout << "\n\n\n\n";
  std::cout << "Finished processing\n";
}
