#include <cstddef>
#include <filesystem>
#include <iterator>

#include "ROI/ButterflyROI.hpp"

int main() {
  std::string treeName = "MyTree";
  std::string pathToInDir =
      "/home/romanurmanov/work/E320/E320Prototype/"
      "E320Prototype_dataInRootFormat/local_alignment_data/raw/sig";
  std::string pathToOutDir =
      "/home/romanurmanov/work/E320/E320Prototype/E320Prototype_preprocessing/"
      "out_data/rois/local_alignment_data";

  std::map<std::string, ButterflyROIFactory::Options> options{
      {"data_Run856", ButterflyROIFactory::Options{.nSteps = 10,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 2.5,
                                                   .applyNoiseFilter = false,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 1000,
                                                   .chipIds = {0, 2, 4, 6, 8}}},
      {"data_Run857", ButterflyROIFactory::Options{.nSteps = 10,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 2.5,
                                                   .applyNoiseFilter = false,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 20,
                                                   .chipIds = {0, 2, 4, 6, 8}}},
      {"data_Run858", ButterflyROIFactory::Options{.nSteps = 10,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 1.5,
                                                   .applyNoiseFilter = true,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 5,
                                                   .chipIds = {0, 2, 4, 6, 8}}},
      {"data_Run948", ButterflyROIFactory::Options{.nSteps = 5,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 2,
                                                   .applyNoiseFilter = false,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 20,
                                                   .chipIds = {0, 2, 4, 6, 8}}},
      {"data_Run949", ButterflyROIFactory::Options{.nSteps = 10,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 1.5,
                                                   .applyNoiseFilter = false,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 20,
                                                   .chipIds = {0, 2, 4, 6, 8}}},
      {"data_Run953", ButterflyROIFactory::Options{.nSteps = 10,
                                                   .stepWidth = 10,
                                                   .stepSize = 10,
                                                   .nStdDevsWidth = 1.5,
                                                   .applyNoiseFilter = true,
                                                   .noiseQuantileLow = 0.01,
                                                   .noiseQuantileHigh = 0.99,
                                                   .nStdDevsNoise = 8,
                                                   .chipIds = {0, 2, 4, 6, 8}}}

  };

  ButterflyROIFactory::Options defaultOpt{};
  defaultOpt.nSteps = 10;
  defaultOpt.stepWidth = 10;
  defaultOpt.stepSize = 10;
  defaultOpt.nStdDevsWidth = 2;
  defaultOpt.applyNoiseFilter = false;
  defaultOpt.noiseQuantileLow = 0.01;
  defaultOpt.noiseQuantileHigh = 0.99;
  defaultOpt.nStdDevsNoise = 20;
  defaultOpt.chipIds = {0, 2, 4, 6, 8};
  std::size_t idx = 0;
  for (const auto& entry : std::filesystem::directory_iterator(pathToInDir)) {
    auto path = entry.path();
    std::size_t n = std::distance(path.begin(), path.end());
    auto it = path.begin();
    std::advance(it, n - 1);
    std::string prefix = *it;
    if (prefix != "data_Run858") {
      // continue;
    }
    ButterflyROIFactory factory;

    ButterflyROIFactory::Options opt =
        (options.contains(prefix)) ? options.at(prefix) : defaultOpt;
    auto rois =
        factory.buildROIForRun(entry.path().string(), treeName, prefix, opt);
    std::string outPath = pathToOutDir + "/" + prefix + ".csv";
    factory.dumpToCsvForRun(outPath, rois);

    factory.buildROIForRun(prefix + ".csv");

    idx++;
    if (idx > 2) {
      // break;
    }
  }

  return 0;
}
