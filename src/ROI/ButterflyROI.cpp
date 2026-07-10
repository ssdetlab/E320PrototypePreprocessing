#include "ROI/ButterflyROI.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "DetectorEvent.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"

double pol1FitFunc(double *x, double *par) {
  return par[0] + par[1] * x[0];
}

double profileFitFunc(double *xx, double *par) {
  double x = xx[0];

  double mean = par[0];
  double sigma = par[1];
  double a = par[2];
  double b = par[3];
  double alpha = par[4];
  return alpha / (sigma * std::sqrt(2 * M_PI)) *
             std::exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma)) +
         (1 - alpha) * (a + b * x);
}

double curvatureFitFunc(double *xx, double *par) {
  double x = xx[0];

  double centerX = par[0];
  double centerY = par[1];
  double width = par[2];
  double curv = par[3];
  return centerY +
         width * std::sqrt(1 + ((x - centerX) * (x - centerX)) / (curv * curv));
}

ButterflyROI::ButterflyROI(double x0, double y0, double w, double theta,
                           double curv)
    : m_x0(x0), m_y0(y0), m_w(w), m_theta(theta), m_curv(curv) {}

bool ButterflyROI::isInside(double x, double y) const {
  // Translate point to ellipse center
  double dx = x - m_x0;
  double dy = y - m_y0;

  // Precompute trig
  double cosT = std::cos(m_theta);
  double sinT = std::sin(m_theta);

  // Rotate into ellipse frame (inverse rotation)
  double xPrime = cosT * dx + sinT * dy;
  double yPrime = -sinT * dx + cosT * dy;

  double boundY = m_w * std::sqrt(1 + (xPrime / m_curv) * (xPrime / m_curv));

  return (std::abs(yPrime) < boundY);
}

ButterflyROIFactory::ROIContainer ButterflyROIFactory::buildROIForRun(
    const std::string &pathToDir, const std::string &treeName,
    const std::string &prefix, const Options &opt) const {
  // Initialize the ROI container
  std::map<std::size_t, std::shared_ptr<ButterflyROI>> rois;

  // Initialize hists, graphs and funcs
  std::map<std::size_t, Stats> roiStats;
  for (std::size_t id : opt.chipIds) {
    Stats stats;
    std::string rawHistName = prefix + "RawHist" + std::to_string(id);
    stats.rawHist = TH2D(rawHistName.c_str(), rawHistName.c_str(), 256, 0, 1024,
                         128, 0, 512);

    std::string slopeFitFuncName = prefix + "SlopeFitFunc" + std::to_string(id);
    stats.slopeGr = TGraph();
    stats.slopeFitFunc = TF1(slopeFitFuncName.c_str(), pol1FitFunc, 0, 1024, 2);

    std::string decorrHistName =
        prefix + "DecorrelatedHist" + std::to_string(id);
    stats.decorrHist = TH2D(decorrHistName.c_str(), decorrHistName.c_str(), 256,
                            0, 1024, 128, 0, 512);

    std::string curvFitFuncName =
        prefix + "CurvatureFitFunc" + std::to_string(id);
    stats.curvGr = TGraph();
    stats.curvFitFunc =
        TF1(curvFitFuncName.c_str(), curvatureFitFunc, 0, 1024, 4);

    std::string testHistName = prefix + "TestHist" + std::to_string(id);
    stats.testHist = TH2D(testHistName.c_str(), testHistName.c_str(), 256, 0,
                          1024, 128, 0, 512);

    roiStats.insert({id, stats});
  };

  // Get ROI data
  std::cout << "Collecting ROI paths from: " << pathToDir << "\n";
  std::vector<std::string> paths;
  for (const auto &entry :
       std::filesystem::recursive_directory_iterator(pathToDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".root") {
      continue;
    }
    paths.push_back(entry.path().string());
  }

  TFile *file = nullptr;
  TTree *tree = nullptr;
  TChain *chainOwner = new TChain("MyTree", "MyTree");
  if (paths.size() == 1) {
    file = new TFile(paths.at(0).c_str());
    tree = file->Get<TTree>("MyTree");
  } else {
    chainOwner = new TChain("MyTree");
    for (const auto &path : paths) {
      chainOwner->Add(path.c_str());
    }
    tree = dynamic_cast<TTree *>(chainOwner);
  }

  detector_event *detEvent = nullptr;
  tree->SetBranchAddress("event", &detEvent);
  for (std::size_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    for (const auto &staveEv : detEvent->st_ev_buffer) {
      for (const auto &chipEv : staveEv.ch_ev_buffer) {
        for (const auto &[hitX, hitY] : chipEv.hits) {
          roiStats.at(chipEv.chip_id).rawHist.Fill(hitX, hitY);
        }
      }
    }
  }

  // Remove noisy pixels
  double quantilePos[2] = {opt.noiseQuantileLow, opt.noiseQuantileHigh};
  for (auto &[id, stats] : roiStats) {
    auto &rawHist = stats.rawHist;
    std::string name = prefix + "BinCounts" + std::to_string(id);

    double maxOcc = rawHist.GetMaximum();
    double minOcc = rawHist.GetMinimum();
    TH1D binCounts(name.c_str(), name.c_str(), maxOcc - minOcc, minOcc, maxOcc);

    std::size_t nBins = rawHist.GetNbinsX() * rawHist.GetNbinsY();
    for (std::size_t i = 1; i <= nBins; i++) {
      binCounts.Fill(rawHist.GetBinContent(i));
    }
    stats.rawBinCounts = binCounts;

    if (opt.applyNoiseFilter) {
      double occQs[2];
      binCounts.GetQuantiles(2, occQs, quantilePos);

      for (std::size_t i = 1; i <= binCounts.FindBin(occQs[0]); i++) {
        binCounts.SetBinContent(i, 0);
      }
      for (std::size_t i = binCounts.FindBin(occQs[1]);
           i <= binCounts.GetNbinsX(); i++) {
        binCounts.SetBinContent(i, 0);
      }
    }
    double mean = binCounts.GetMean();
    double sigma = binCounts.GetStdDev();
    stats.trimmedBinCounts = binCounts;

    if (opt.applyNoiseFilter) {
      for (std::size_t i = 1; i <= nBins; i++) {
        if (std::abs(rawHist.GetBinContent(i) - mean) / sigma >
            opt.nStdDevsNoise) {
          rawHist.SetBinContent(i, 0);
        }
      }
    }
  }

  // Iterate over the chips
  for (auto &[id, stats] : roiStats) {
    auto &rawHist = stats.rawHist;
    auto &slopeGr = stats.slopeGr;
    auto &slopeFitFunc = stats.slopeFitFunc;
    auto &decorrHist = stats.decorrHist;
    auto &curvGr = stats.curvGr;
    auto &curvFitFunc = stats.curvFitFunc;
    auto &testHist = stats.testHist;

    // Get the butterfly centers
    int maxBinX = 0, maxBinY = 0, maxBinZ = 0;
    int maxBin = rawHist.GetMaximumBin();
    rawHist.GetBinXYZ(maxBin, maxBinX, maxBinY, maxBinZ);
    double centerX = rawHist.GetXaxis()->GetBinCenter(maxBinX);
    double centerY = rawHist.GetYaxis()->GetBinCenter(maxBinY);

    // Get the butterfly slope
    for (int i = -opt.nSteps; i <= opt.nSteps; i++) {
      std::size_t currentBin = maxBinX + i * opt.stepSize;
      if (currentBin < 1 || currentBin > rawHist.GetNbinsX()) {
        continue;
      }
      std::string name = "proj" + std::to_string(i);
      auto slice =
          rawHist.ProjectionY(name.c_str(), currentBin - opt.stepWidth / 2,
                              currentBin + opt.stepWidth / 2);

      double x = slice->GetBinCenter(maxBinX + i * opt.stepSize);
      double y = slice->GetBinCenter(slice->GetMaximumBin());
      slopeGr.SetPoint(slopeGr.GetN(), x, y);
    }
    slopeFitFunc.SetParameters(10, 0);
    slopeGr.Fit(&slopeFitFunc, "Q");
    double theta = std::atan(slopeFitFunc.GetParameter(1));

    // Decorrelate the butterfly
    for (std::size_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      for (const auto &staveEv : detEvent->st_ev_buffer) {
        for (const auto &chipEv : staveEv.ch_ev_buffer) {
          for (const auto &[hitX, hitY] : chipEv.hits) {
            if (chipEv.chip_id == 0) {
              double x = hitX - centerX;
              double y = hitY - centerY;

              x = std::cos(theta) * x + std::sin(theta) * y;
              y = -std::sin(theta) * x + std::cos(theta) * y;

              x += centerX;
              y += centerY;

              decorrHist.Fill(x, y);
            }
          }
        }
      }
    }

    // Estimate widths
    std::size_t centerSliceIdx = 0;
    for (int i = -opt.nSteps; i <= opt.nSteps; i++) {
      std::size_t currentBin = maxBinX + i * opt.stepSize;
      if (currentBin < 1 || currentBin > rawHist.GetNbinsX()) {
        continue;
      }
      std::string histName =
          std::string(decorrHist.GetName()) + "Proj" + std::to_string(i);
      auto slice = *decorrHist.ProjectionY(histName.c_str(),
                                           currentBin - opt.stepWidth / 2,
                                           currentBin + opt.stepWidth / 2);
      slice.Scale(1.0 / slice.Integral(), "width");

      std::string funcName = histName + "ProfileFitFunc";
      TF1 fitFunc(funcName.c_str(), profileFitFunc, 0, 1024, 5);

      double mean = slice.GetBinCenter(slice.GetMaximumBin());
      double sigma = 6.0;
      double pedestal = 0.002;
      double slope = -1e-6;
      double alpha = 0.2;

      fitFunc.SetParameters(mean, sigma, pedestal, slope, alpha);
      slice.Fit(&fitFunc, "Q");

      stats.slices.push_back(
          {decorrHist.GetXaxis()->GetBinCenter(currentBin), slice, fitFunc});
      if (i == 0) {
        centerSliceIdx = stats.slices.size() - 1;
      }
    }

    // Estimate curvature
    for (const auto &[x, hist, func] : stats.slices) {
      double y =
          func.GetParameter(0) + func.GetParameter(1) * opt.nStdDevsWidth;
      curvGr.SetPoint(curvGr.GetN(), x, y);
    }
    curvFitFunc.FixParameter(0, centerX);
    curvFitFunc.FixParameter(1, centerY);
    curvFitFunc.SetParameter(2, 30);
    curvFitFunc.SetParameter(3, 300);
    curvGr.Fit(&curvFitFunc, "Q");

    double width = curvFitFunc.GetParameter(2);
    double curv = curvFitFunc.GetParameter(3);
    centerX = std::get<0>(stats.slices.at(centerSliceIdx));
    centerY = std::get<2>(stats.slices.at(centerSliceIdx)).GetParameter(0);

    auto roi =
        std::make_shared<ButterflyROI>(centerX, centerY, width, theta, curv);
    rois.insert({id, roi});

    for (std::size_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      for (const auto &staveEv : detEvent->st_ev_buffer) {
        for (const auto &chipEv : staveEv.ch_ev_buffer) {
          for (const auto &[hitX, hitY] : chipEv.hits) {
            if (roi->isInside(hitX, hitY)) {
              testHist.Fill(hitX, hitY);
            }
          }
        }
      }
    }
  }

  std::string cOption = prefix + ".pdf";
  std::string cOptionFirst = prefix + ".pdf(";
  std::string cOptionLast = prefix + ".pdf)";

  std::string cName = prefix + "Canvas";
  TCanvas c(cName.c_str(), cName.c_str(), 800, 800);

  std::string dummyName = prefix + "Dummy";
  TH1D dummy(dummyName.c_str(), dummyName.c_str(), 100, 0, 100);
  dummy.Draw("hist");
  c.Print(cOptionFirst.c_str());

  for (auto &[id, stats] : roiStats) {
    stats.rawHist.Draw("col");
    stats.slopeGr.Draw("PL same");
    stats.slopeGr.SetMarkerStyle(20);
    stats.slopeGr.SetMarkerColor(kRed);
    stats.slopeGr.SetMarkerSize(1.5);
    stats.slopeGr.SetLineColor(kRed);
    stats.slopeGr.SetLineWidth(2);
    c.Print(cOption.c_str());

    stats.testHist.Draw("col");
    c.Print(cOption.c_str());

    stats.decorrHist.Draw("col");
    stats.curvGr.Draw("P same");
    stats.curvGr.SetMarkerStyle(20);
    stats.curvGr.SetMarkerColor(kRed);
    stats.curvGr.SetMarkerSize(1.5);
    stats.curvGr.SetLineColor(kRed);
    stats.curvGr.SetLineWidth(2);
    stats.curvFitFunc.SetNpx(10000);
    stats.curvFitFunc.Draw("same");
    c.Print(cOption.c_str());
  }
  for (auto &[id, stats] : roiStats) {
    stats.rawBinCounts.Draw("hist");
    c.Print(cOption.c_str());

    stats.trimmedBinCounts.Draw("hist");
    c.Print(cOption.c_str());
  }
  for (auto &[id, stats] : roiStats) {
    for (auto &[x, hist, func] : stats.slices) {
      hist.Draw("hist");
      func.SetNpx(10000);
      func.Draw("same");
      c.Print(cOption.c_str());
    }
  }
  for (auto &[id, stats] : roiStats) {
    for (auto &[x, hist, func] : stats.slices) {
      hist.Draw("hist");
      func.SetNpx(10000);
      func.Draw("same");
      c.Print(cOption.c_str());
    }
  }
  dummy.Draw("hist");
  c.Print(cOptionLast.c_str());

  return rois;
}

ButterflyROIFactory::ROIContainer ButterflyROIFactory::buildROIForRun(
    const std::string &pathToCsv) const {
  ROIContainer rois;

  std::ifstream fin;
  fin.open(pathToCsv);

  for (std::string line; std::getline(fin, line);) {
    std::stringstream lineStream(line);
    std::string token;

    // Chip ID
    std::getline(lineStream, token, ',');
    std::size_t chipId = std::stol(token);

    // Butterfly center x
    std::getline(lineStream, token, ',');
    double x0 = std::stod(token);

    // Butterfly center y
    std::getline(lineStream, token, ',');
    double y0 = std::stod(token);

    // Butterfly minor axis size
    std::getline(lineStream, token, ',');
    double w = std::stod(token);

    // Butterfly rotation angle
    std::getline(lineStream, token, ',');
    double theta = std::stod(token);

    // Butterfly curvature
    std::getline(lineStream, token, ',');
    double curv = std::stod(token);

    rois.insert(
        {chipId, std::make_shared<ButterflyROI>(x0, y0, w, theta, curv)});
  }
  fin.close();
  return rois;
}

void ButterflyROIFactory::dumpToCsvForRun(const std::string &path,
                                          const ROIContainer &rois) const {
  std::ofstream fout;
  fout.open(path);
  for (const auto &[id, roi] : rois) {
    fout << id << "," << roi->centerX() << "," << roi->centerY() << ","
         << roi->width() << "," << roi->theta() << "," << roi->curvature()
         << "\n";
  }
  fout.close();
}
