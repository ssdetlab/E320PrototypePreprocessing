#include <filesystem>
#include <memory>
#include <string>

#include "Analysis/ClusteringAlgorithm.hpp"
#include "Analysis/PreprocessorMap.hpp"
#include "ROI/ROIProvider.hpp"

int main() {
  // Initialize clustering algorithm
  ClusteringAlgorithm::Config clusteringAlgoCfg{};
  // Maximum cluster size
  clusteringAlgoCfg.maxClusterSize = 4;
  auto clusteringAlgo =
      std::make_shared<ClusteringAlgorithm>(clusteringAlgoCfg);

  // Initialize ROI provider
  ROIProvider::Config roiProviderCfg{};
  // ROI data directory path
  roiProviderCfg.inROIDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_preprocessing/"
      "out_data/rois/local_alignment_data";
  // ROI type
  roiProviderCfg.roiType = ROIProvider::ROIType::Butterfly;
  // Chip IDs of the dataset
  roiProviderCfg.chipIds = {0, 2, 4, 6, 8};
  // Run prefix to ROI run prefix map,
  // every run with a ROI file uses its own ROIs
  for (const auto &entry :
       std::filesystem::directory_iterator(roiProviderCfg.inROIDir)) {
    if (entry.path().extension() == ".csv") {
      std::string runPrefix = entry.path().stem().string();
      roiProviderCfg.roiRunMap.insert({runPrefix, runPrefix});
    }
  }
  // ROIs of another run, e.g. estimated on signal and applied to background
  // roiProviderCfg.roiRunMap.insert_or_assign("data_RunBkg", "data_RunSig");
  auto roiProvider = std::make_shared<ROIProvider>(roiProviderCfg);

  Preprocessor preprocessor;

  Preprocessor::Options opt;

  // Data directory path
  opt.inDataDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/local_alignment_data/raw/sig";
  // Orbit reconstruction data directory path
  opt.inOrbitDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_bpm_analysis/"
      "out_data/data/global_2026_data";
  // Data tree name
  opt.inDataTreeName = "MyTree";
  // Orbit reconstruction tree name
  opt.inOrbitTreeName = "fitted-beam";
  // Output directory path
  opt.outDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/local_alignment_data/processed";
  // Output data tree name
  opt.outDataTreeName = "data";
  // Output magnet tree name
  opt.outMagnetTreeName = "magnets";
  // Chip IDs of the dataset
  opt.chipIds = {0, 2, 4, 6, 8};
  // Magnets parameters
  opt.dipoleLength = 0.9144;       // in meters
  opt.xCorrectorLength = 0.23368;  // in meters
  opt.deflectionAngle = 0.006;     // in rad
  /// Accelerator data noise low quantile
  opt.noiseQuantileLowAcc = 0.1;
  /// Accelerator data noise high quantile
  opt.noiseQuantileHighAcc = 0.9;
  // Max number of StdDevs for the
  // accelerator data
  opt.nStdDevsAcc = 6;
  // Maximum number of hits
  opt.maxHits = 1e6;
  // Maximum number of hits in the ROI
  opt.maxRoiHits = 1e6;
  // Exclude shuttered events flag
  opt.excludeShuttered = false;
  opt.clusteringAlgo = clusteringAlgo;
  opt.roiProvider = roiProvider;

  preprocessor.processFiles(opt);
};
