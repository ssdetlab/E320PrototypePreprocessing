#include <filesystem>
#include <memory>
#include <string>

#include "Analysis/AcceleratorStateFilter.hpp"
#include "Analysis/ClusteringAlgorithm.hpp"
#include "Analysis/PVList.hpp"
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

  // Initialize accelerator state filter
  AcceleratorStateFilter::Config acceleratorStateFilterCfg{};
  // Accelerator data noise low quantile
  acceleratorStateFilterCfg.noiseQuantileLow = 0.1;
  // Accelerator data noise high quantile
  acceleratorStateFilterCfg.noiseQuantileHigh = 0.9;
  // Max number of StdDevs for the
  // accelerator data
  acceleratorStateFilterCfg.nStdDevs = 6;
  // Accelerator state PVs: PV name, label, cut type,
  // number of bins, lower and upper edge of the distribution
  acceleratorStateFilterCfg.pvs = {
      {pmt_s20_3060_key, "PMT_S20_3060", AcceleratorStateFilter::CutType::Upper,
       10000, 0, 50000},
      {pmt_s20_3070_key, "PMT_S20_3070", AcceleratorStateFilter::CutType::Upper,
       10000, 0, 50000},
      {pmt_s20_3179_key, "PMT_S20_3179", AcceleratorStateFilter::CutType::Upper,
       10000, 0, 50000},
      {pmt_s20_3350_key, "PMT_S20_3350", AcceleratorStateFilter::CutType::Upper,
       10000, 0, 50000},
      {pmt_s20_3360_key, "PMT_S20_3360", AcceleratorStateFilter::CutType::Upper,
       10000, 0, 50000},
      {toro_li20_1988_tmit_pc_key, "TORO_LI20_1988_TMIT_PC",
       AcceleratorStateFilter::CutType::Central, 10000, -5000, 5000},
      {toro_li20_2040_tmit_pc_key, "TORO_LI20_2040_TMIT_PC",
       AcceleratorStateFilter::CutType::Central, 10000, -5000, 5000},
      {toro_li20_2452_tmit_pc_key, "TORO_LI20_2452_TMIT_PC",
       AcceleratorStateFilter::CutType::Central, 10000, -5000, 5000},
      {toro_li20_3163_tmit_pc_key, "TORO_LI20_3163_TMIT_PC",
       AcceleratorStateFilter::CutType::Central, 10000, -5000, 5000},
      {toro_li20_3255_tmit_pc_key, "TORO_LI20_3255_TMIT_PC",
       AcceleratorStateFilter::CutType::Central, 10000, -5000, 5000}};
  auto acceleratorStateFilter =
      std::make_shared<AcceleratorStateFilter>(acceleratorStateFilterCfg);

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
  // Maximum number of hits
  opt.maxHits = 1e6;
  // Maximum number of hits in the ROI
  opt.maxRoiHits = 1e6;
  // Exclude shuttered events flag
  opt.excludeShuttered = false;
  opt.clusteringAlgo = clusteringAlgo;
  opt.roiProvider = roiProvider;
  opt.acceleratorStateFilter = acceleratorStateFilter;

  preprocessor.processFiles(opt);
};
