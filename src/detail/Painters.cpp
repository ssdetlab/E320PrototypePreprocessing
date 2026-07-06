#include "detail/Painters.hpp"

#include <csignal>
#include <cstddef>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveStats.h"

void applyGraphPainterConfig(TGraph &gr, const GraphPainterConfig &cfg) {
  gr.Draw(cfg.option.c_str());
  gr.SetStats(cfg.stats);
  if (!cfg.printTitle) {
    gr.SetTitle("");
  }
  gr.GetXaxis()->SetTitle(cfg.xTitle.c_str());
  gr.GetYaxis()->SetTitle(cfg.yTitle.c_str());
  gr.GetXaxis()->SetTitleSize(cfg.xTitleSize);
  gr.GetXaxis()->SetTitleOffset(cfg.xTitleOffset);
  gr.GetXaxis()->SetLabelSize(cfg.xLabelSize);
  gr.GetYaxis()->SetTitleSize(cfg.yTitleSize);
  gr.GetYaxis()->SetTitleOffset(cfg.yTitleOffset);
  gr.GetYaxis()->SetLabelSize(cfg.yLabelSize);
  gr.SetLineColor(cfg.color);
  gr.SetLineWidth(cfg.lineWidth);
  gr.SetMarkerColor(cfg.color);
  gr.SetMarkerStyle(cfg.markerStyle);
  gr.SetMarkerSize(cfg.markerSize);

  if (cfg.xAdjustDivisions) {
    gr.GetXaxis()->SetNdivisions(cfg.xNDivisions);
  }
  if (cfg.yAdjustDivisions) {
    gr.GetXaxis()->SetNdivisions(cfg.yNDivisions);
  }
}

void applyHistPainterConfig(TH1 &hist, const HistPainterConfig &cfg) {
  hist.Draw(cfg.option.c_str());
  hist.SetStats(cfg.stats);
  if (!cfg.printTitle) {
    hist.SetTitle("");
  }
  hist.GetXaxis()->SetTitle(cfg.xTitle.c_str());
  hist.GetYaxis()->SetTitle(cfg.yTitle.c_str());
  hist.GetXaxis()->SetTitleSize(cfg.xTitleSize);
  hist.GetXaxis()->SetTitleOffset(cfg.xTitleOffset);
  hist.GetXaxis()->SetLabelSize(cfg.xLabelSize);
  hist.GetYaxis()->SetTitleSize(cfg.yTitleSize);
  hist.GetYaxis()->SetTitleOffset(cfg.yTitleOffset);
  hist.GetYaxis()->SetLabelSize(cfg.yLabelSize);
  hist.SetLineColor(cfg.color);
  hist.SetFillColorAlpha(cfg.color, cfg.colorAlpha);
  hist.SetLineWidth(cfg.lineWidth);

  if (cfg.stats) {
    gPad->Update();
    const auto &statsCfg = cfg.statsCfg;
    TPaveStats *stats = dynamic_cast<TPaveStats *>(hist.FindObject("stats"));
    stats->SetOptStat(statsCfg.getStatBoxOpt());
    stats->SetX1NDC(statsCfg.xStart);
    stats->SetY1NDC(statsCfg.yStart);
    stats->SetX2NDC(statsCfg.xEnd);
    stats->SetY2NDC(statsCfg.yEnd);
  }

  if (cfg.xAdjustDivisions) {
    hist.GetXaxis()->SetNdivisions(cfg.xNDivisions);
  }
  if (cfg.yAdjustDivisions) {
    hist.GetXaxis()->SetNdivisions(cfg.yNDivisions);
  }
}

void saveGraph(TCanvas &c, const std::string &cOption,
               const std::vector<std::pair<TGraph &, GraphPainterConfig>> &grs,
               const LegendPainterConfig &legCfg) {
  const auto &[topLevelGr, topLevelCfg] = grs.front();

  TLegend *leg = nullptr;
  leg = new TLegend(legCfg.xBottomLeft, legCfg.yBottomLeft, legCfg.xTopRight,
                    legCfg.yTopRight);

  c.GetPad(1);
  c.SetMargin(topLevelCfg.padLeftMargin, topLevelCfg.padRightMargin,
              topLevelCfg.padBottomMargin, topLevelCfg.padTopMargin);

  c.SetGridx(1);
  c.SetGridy(1);

  c.SetLogx(0);
  c.SetLogy(0);

  for (std::size_t i = 0; i < grs.size(); i++) {
    const auto &[gr, cfg] = grs.at(i);
    std::string option = (i == 0) ? cfg.option : cfg.option + " same";
    std::size_t pos = option.find("A");
    if (i != 0 && pos != std::string::npos) {
      option.erase(pos);
    }

    applyGraphPainterConfig(gr, cfg);

    leg->AddEntry(&gr, cfg.legTitle.c_str(), legCfg.option.c_str());
  }

  if (topLevelCfg.xAdjustDivisions) {
    topLevelGr.GetXaxis()->SetNdivisions(topLevelCfg.xNDivisions);
  }
  if (topLevelCfg.yAdjustDivisions) {
    topLevelGr.GetXaxis()->SetNdivisions(topLevelCfg.yNDivisions);
  }

  if (legCfg.draw) {
    leg->Draw("same");
  }
  c.SetCanvasSize(topLevelCfg.canvasXSize, topLevelCfg.canvasYSize);
  c.Print(cOption.c_str(), "pdf");
}

void saveHist(TCanvas &c, const std::string &cOption,
              const std::vector<std::pair<TH1 &, HistPainterConfig>> &hists) {
  const auto &[topLevelHist, topLevelCfg] = hists.front();

  c.GetPad(1);
  c.SetMargin(topLevelCfg.padLeftMargin, topLevelCfg.padRightMargin,
              topLevelCfg.padBottomMargin, topLevelCfg.padTopMargin);

  c.SetGridx(static_cast<Int_t>(topLevelCfg.gridX));
  c.SetGridy(static_cast<Int_t>(topLevelCfg.gridY));

  c.SetLogx(static_cast<Int_t>(topLevelCfg.logX));
  c.SetLogy(static_cast<Int_t>(topLevelCfg.logY));

  for (std::size_t i = 0; i < hists.size(); i++) {
    const auto &[hist, cfg] = hists.at(i);
    std::string option = (i == 0) ? cfg.option : cfg.option + " same";
    std::size_t pos = option.find("A");
    if (i != 0 && pos != std::string::npos) {
      option.erase(pos);
    }

    applyHistPainterConfig(hist, cfg);
  }

  if (topLevelCfg.xAdjustDivisions) {
    topLevelHist.GetXaxis()->SetNdivisions(topLevelCfg.xNDivisions);
  }
  if (topLevelCfg.yAdjustDivisions) {
    topLevelHist.GetXaxis()->SetNdivisions(topLevelCfg.yNDivisions);
  }
  c.SetCanvasSize(topLevelCfg.canvasXSize, topLevelCfg.canvasYSize);
  c.Print(cOption.c_str(), "pdf");
}
