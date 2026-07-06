#pragma once

#include <csignal>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

struct HistPainterConfig {
  // Canvas
  bool gridX = true;
  bool gridY = true;

  bool logX = false;
  bool logY = false;

  double canvasXSize = 2000;
  double canvasYSize = 2000;

  double padLeftMargin = 0.1;
  double padRightMargin = 0.1;
  double padBottomMargin = 0.1;
  double padTopMargin = 0.1;

  bool stats = false;
  bool printTitle = true;

  // Axes
  std::string xTitle = "";
  double xTitleSize = 0.04;
  double xTitleOffset = 1.0;
  double xLabelSize = 0.03;
  int xNDivisions = 10;
  bool xAdjustDivisions = false;

  std::string yTitle = "";
  double yTitleSize = 0.04;
  double yTitleOffset = 1.0;
  int yNDivisions = 10;
  bool yAdjustDivisions = false;

  double yLabelSize = 0.03;

  // Hist
  EColor color = kBlue;
  double colorAlpha = 0.3;
  double lineWidth = 1;
  std::string option = "hist";

  // Legend
  std::string legTitle = "";

  // Stats
  struct StatBoxConfig {
    enum EntryOpt : int { off = 0, on = 1, onerr = 2 };
    EntryOpt kurtosis = off;
    EntryOpt skewness = off;
    EntryOpt integral = off;
    EntryOpt overflow = off;
    EntryOpt underflow = off;
    EntryOpt rms = on;
    EntryOpt mean = on;
    EntryOpt entries = on;
    EntryOpt name = on;

    double xStart = 0.7;
    double yStart = 0.7;
    double xEnd = 0.9;
    double yEnd = 0.9;

    int getStatBoxOpt() const {
      int opt = 0;

      int k = kurtosis * static_cast<int>(1e8);
      opt += (k == 0 && opt == 0) ? 0 : k;
      int s = skewness * static_cast<int>(1e7);
      opt += (s == 0 && opt == 0) ? 0 : s;
      int i = integral * static_cast<int>(1e6);
      opt += (i == 0 && opt == 0) ? 0 : i;
      int o = overflow * static_cast<int>(1e5);
      opt += (o == 0 && opt == 0) ? 0 : o;
      int u = underflow * static_cast<int>(1e4);
      opt += (u == 0 && opt == 0) ? 0 : u;
      int r = rms * static_cast<int>(1e3);
      opt += (r == 0 && opt == 0) ? 0 : r;
      int m = mean * static_cast<int>(1e2);
      opt += (m == 0 && opt == 0) ? 0 : m;
      int e = entries * static_cast<int>(1e1);
      opt += (e == 0 && opt == 0) ? 0 : e;
      int n = name * static_cast<int>(1e0);
      opt += (n == 0 && opt == 0) ? 0 : n;
      return opt;
    }
  };
  StatBoxConfig statsCfg{};
};

struct GraphPainterConfig {
  // Canvas
  bool gridX = true;
  bool gridY = true;

  bool logX = false;
  bool logY = false;

  double canvasXSize = 2000;
  double canvasYSize = 2000;

  double padLeftMargin = 0.1;
  double padRightMargin = 0.1;
  double padBottomMargin = 0.1;
  double padTopMargin = 0.1;

  bool stats = false;
  bool printTitle = true;

  // Axes
  std::string xTitle = "";
  double xTitleSize = 0.04;
  double xTitleOffset = 1.0;
  double xLabelSize = 0.03;
  int xNDivisions = 10;
  bool xAdjustDivisions = false;

  std::string yTitle = "";
  double yTitleSize = 0.04;
  double yTitleOffset = 1.0;
  int yNDivisions = 10;
  bool yAdjustDivisions = false;

  double yLabelSize = 0.03;

  // Graph
  EColor color = kBlue;
  double lineWidth = 1;
  double markerSize = 2;
  int markerStyle = 20;
  std::string option = "APL";
  std::string legTitle = "";
};

struct LegendPainterConfig {
  bool draw = true;

  double xBottomLeft = 0.7;
  double yBottomLeft = 0.7;
  double xTopRight = 0.9;
  double yTopRight = 0.9;

  std::string option = "p";
};

void applyGraphPainterConfig(TGraph &gr, const GraphPainterConfig &cfg);

void applyHistPainterConfig(TH2 &hist, const HistPainterConfig &cfg);

void applyHistPainterConfig(TH1 &hist, const HistPainterConfig &cfg);

void saveGraph(TCanvas &c, const std::string &cOption,
               const std::vector<std::pair<TGraph &, GraphPainterConfig>> &grs,
               const LegendPainterConfig &legCfg);

void saveHist(TCanvas &c, const std::string &cOption,
              const std::vector<std::pair<TH1 &, HistPainterConfig>> &hists);
