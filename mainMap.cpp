#include "Analysis/PreprocessorMap.hpp"

int main() {
  Preprocessor preprocessor;

  Preprocessor::Options opt;

  // Data directory path
  opt.pathToDataDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_June_2026_positrons/raw/"
      "positron_runs/june_23";
  // ROI data directory path
  opt.pathToROIDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_June_2026_positrons/shift_copy/"
      "misc/june_23/data_Run1207/ev_data";
  // Data tree name
  opt.dataTreeName = "MyTree";
  // ROI data tree name
  opt.roiDataTreeName = "MyTree";
  // Output directory path
  opt.outPathDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/E320Shift_June_2026_positrons/processed/"
      "positron_runs/june_23";
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
  // Max number of StdDevs for the
  // accelerator data
  opt.nStdDevsAcc = 6;
  // ROI type
  opt.roiType = Preprocessor::ROIType::Butterfly;
  // Maximum number of hits
  opt.maxHits = 35;
  // Maximum number of hits in the ROI
  opt.maxRoiHits = 1e6;
  // Maximum cluster size
  opt.maxClusterSize = 4;

  preprocessor.processFiles(opt);
};
