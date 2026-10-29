#include "Analysis/PreprocessorSim.hpp"

int main() {
  PreprocessorSim preprocessor;

  PreprocessorSim::Options opt;

  // Data directory path
  opt.inTrackerDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/sim/test/alignment/clusters_tracker/raw";
  opt.inBeamDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_bpm_analysis/"
      "out_data/sim/test/alignment/orbit_est";

  // Data tree name
  opt.inTrackerTreeName = "clusters";
  opt.inBeamTreeName = "fitted-beam";

  // Output directory path
  opt.outPath =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/sim/clusters.root";
  // Output data tree name
  opt.outTreeName = "clusters";

  preprocessor.processFiles(opt);
};
