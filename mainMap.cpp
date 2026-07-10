#include "Analysis/PreprocessorMap.hpp"

int main() {
  Preprocessor preprocessor;

  Preprocessor::Options opt;

  // Data directory path
  opt.inDataDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/local_alignment_data/raw/sig";
  // ROI data directory path
  opt.inROIDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_preprocessing/"
      "out_data/rois/local_alignment_data";
  // Data tree name
  opt.inDataTreeName = "MyTree";
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
  // ROI type
  opt.roiType = Preprocessor::ROIType::Butterfly;
  // Maximum number of hits
  opt.maxHits = 1e6;
  // Maximum number of hits in the ROI
  opt.maxRoiHits = 1e6;
  // Maximum cluster size
  opt.maxClusterSize = 4;
  // Exclude shuttered events flag
  opt.excludeShuttered = false;

  preprocessor.processFiles(opt);
};
