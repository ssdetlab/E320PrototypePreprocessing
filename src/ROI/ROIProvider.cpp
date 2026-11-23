#include "ROI/ROIProvider.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "ROI/ButterflyROI.hpp"
#include "ROI/UnboundedROI.hpp"

ROIProvider::ROIProvider(const Config &cfg) : m_cfg(cfg) {
  switch (m_cfg.roiType) {
    case Ellipse: {
      throw std::runtime_error("Ellipse ROI reading is not supported yet");
      break;
    }
    case Butterfly: {
      // Read every ROI run once, runs mapped
      // to the same ROI run share its ROIs
      ButterflyROIFactory factory;
      std::map<std::string, ROIContainer> roiRunROIs;
      for (const auto &[runPrefix, roiRunPrefix] : m_cfg.roiRunMap) {
        if (!roiRunROIs.contains(roiRunPrefix)) {
          std::filesystem::path roiPath =
              std::filesystem::path(m_cfg.inROIDir) / (roiRunPrefix + ".csv");
          if (!std::filesystem::is_regular_file(roiPath)) {
            throw std::runtime_error("Missing ROI file " + roiPath.string() +
                                     " for run " + runPrefix);
          }
          auto butterflies = factory.buildROIForRun(roiPath.string());

          ROIContainer rois;
          for (std::size_t id : m_cfg.chipIds) {
            if (!butterflies.contains(id)) {
              throw std::runtime_error("Missing ROI for chip " +
                                       std::to_string(id) + " in " +
                                       roiPath.string());
            }
            rois.insert({id, butterflies.at(id)});
          }
          roiRunROIs.insert({roiRunPrefix, rois});
        }
        m_rois.insert({runPrefix, roiRunROIs.at(roiRunPrefix)});
        std::cout << "Run " << runPrefix << " uses the ROIs of " << roiRunPrefix
                  << "\n";
      }
      break;
    }
    case Unbounded: {
      for (std::size_t id : m_cfg.chipIds) {
        m_unboundedROIs.insert({id, std::make_shared<UnboundedROI>()});
      }
      break;
    }
    default: {
      throw std::runtime_error("Unknown ROI type");
      break;
    }
  }
}

const ROIProvider::ROIContainer &ROIProvider::getROIs(
    const std::string &runPrefix) const {
  if (m_cfg.roiType == Unbounded) {
    return m_unboundedROIs;
  }

  auto it = m_rois.find(runPrefix);
  if (it == m_rois.end()) {
    throw std::runtime_error("No ROI run configured for run " + runPrefix);
  }
  return it->second;
}
