#include "Analysis/Stats.hpp"

#include <string>

#include "detail/Painters.hpp"

Stats::Stats(const std::set<int> &ids) {
  chipIds = ids;
  for (int id : chipIds) {
    // Raw hit map ROI hists
    std::string nameChipHitMapROIRaw =
        "Chip" + std::to_string(id) + "HitMapROIRaw";
    chipHitMapROIRawHists.insert(
        {id, TH2D(nameChipHitMapROIRaw.c_str(), nameChipHitMapROIRaw.c_str(),
                  200, 0, 1024, 200, 0, 512)});
  }
}

void Stats::init(const std::string &prefix) {
  runPrefix = prefix;

  std::string nameEventSplit = runPrefix + "EventSplit";
  eventCategotyHist =
      TH1I(nameEventSplit.c_str(), nameEventSplit.c_str(), 6, 0, 6);
  eventCategotyHist.Fill("nEventsTotal", 0);
  eventCategotyHist.Fill("nEventsAccepted", 0);
  eventCategotyHist.Fill("nEventsUnshuttered", 0);
  eventCategotyHist.Fill("nEventsEvenPID", 0);
  eventCategotyHist.Fill("nEventsOddPID", 0);
  eventCategotyHist.Fill("nEventsFullDetEvenPID", 0);
  eventCategotyHist.Fill("nEventsFullDetOddPID", 0);

  for (int id : chipIds) {
    // Raw hit map hists
    std::string nameChipHitMapRaw =
        runPrefix + "Chip" + std::to_string(id) + "HitMapRaw";
    chipHitMapRawHists.insert(
        {id, TH2D(nameChipHitMapRaw.c_str(), nameChipHitMapRaw.c_str(), 200, 0,
                  1024, 200, 0, 512)});
    // Filtered hit map hists
    std::string nameChipHitMapFiltered =
        runPrefix + "Chip" + std::to_string(id) + "HitMapFiltered";
    chipHitMapFilteredHists.insert(
        {id, TH2D(nameChipHitMapFiltered.c_str(),
                  nameChipHitMapFiltered.c_str(), 200, 0, 1024, 200, 0, 512)});
    // Raw occupancy graphs
    std::string nameChipOccR =
        runPrefix + "Chip" + std::to_string(id) + "OccRaw";
    chipOccRawGrs.insert({id, TGraph()});
    // Filtered occupancy graphs
    std::string nameChipOccFilter =
        runPrefix + "Chip" + std::to_string(id) + "OccFiltered";
    chipOccFilteredGrs.insert({id, TGraph()});
  }

  // Raw PMTs
  gr_PMT_S20_3060_Raw = TGraph();
  gr_PMT_S20_3070_Raw = TGraph();
  gr_PMT_S20_3179_Raw = TGraph();
  gr_PMT_S20_3350_Raw = TGraph();
  gr_PMT_S20_3360_Raw = TGraph();
  // Filtered PMTs
  gr_PMT_S20_3060_Filtered = TGraph();
  gr_PMT_S20_3070_Filtered = TGraph();
  gr_PMT_S20_3179_Filtered = TGraph();
  gr_PMT_S20_3350_Filtered = TGraph();
  gr_PMT_S20_3360_Filtered = TGraph();
  // Raw toroids
  gr_TORO_LI20_1988_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_2040_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_2452_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_3163_TMIT_PC_Raw = TGraph();
  gr_TORO_LI20_3255_TMIT_PC_Raw = TGraph();
  // Filtered toroids
  gr_TORO_LI20_1988_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_2040_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_2452_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_3163_TMIT_PC_Filtered = TGraph();
  gr_TORO_LI20_3255_TMIT_PC_Filtered = TGraph();
}

void Stats::reset() {
  // N total events
  nEventsTotal = 0;
  // N accepted events
  nEventsAccepted = 0;
  // N shuttered events
  nEventsUnshuttered = 0;
  // N even PID events
  nEventsEvenPID = 0;
  // N odd PID events
  nEventsOddPID = 0;
  // N events with full detector occupancy
  // and even PID
  nEventsFullDetEvenPID = 0;
  // N events with full detector occupancy
  // and odd PID
  nEventsFullDetOddPID = 0;

  // Raw hit map hists
  chipHitMapRawHists.clear();
  // Filtered hit map hists
  chipHitMapFilteredHists.clear();
  // Raw occupancy graphs
  chipOccRawGrs.clear();
  // Filtered occupancy graphs
  chipOccFilteredGrs.clear();
}

void Stats::print() {
  TCanvas c = TCanvas(runPrefix.c_str(), "", 800, 800);
  c.GetPad(1);
  c.SetMargin(0.1, 0.1, 0.1, 0.1);

  c.SetGridx(1);
  c.SetGridy(1);

  c.SetLogx(0);
  c.SetLogy(0);

  LegendPainterConfig legCfg;
  legCfg.draw = true;
  legCfg.option = "p";

  std::string cOptionFirst = runPrefix + ".pdf(";
  std::string cOptionLast = runPrefix + ".pdf)";
  std::string cOption = runPrefix + ".pdf";

  HistPainterConfig eventSplitCfg;
  eventSplitCfg.yTitle = "Count";
  eventSplitCfg.option = "hist text";

  eventCategotyHist.Fill("nEventsTotal", nEventsTotal);
  eventCategotyHist.Fill("nEventsAccepted", nEventsAccepted);
  eventCategotyHist.Fill("nEventsUnshuttered", nEventsUnshuttered);
  eventCategotyHist.Fill("nEventsEvenPID", nEventsEvenPID);
  eventCategotyHist.Fill("nEventsOddPID", nEventsOddPID);
  eventCategotyHist.Fill("nEventsFullDetEvenPID", nEventsFullDetEvenPID);
  eventCategotyHist.Fill("nEventsFullDetOddPID", nEventsFullDetOddPID);
  saveHist(c, cOptionFirst, {{eventCategotyHist, eventSplitCfg}});

  // PMT_S20_3060
  GraphPainterConfig gr_PMT_S20_3060_Raw_Cfg;
  gr_PMT_S20_3060_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3060_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3060_Raw.SetNameTitle("gr_PMT_S20_3060", "gr_PMT_S20_3060");
  GraphPainterConfig gr_PMT_S20_3060_Filtered_Cfg;
  gr_PMT_S20_3060_Filtered_Cfg.option = "P";
  gr_PMT_S20_3060_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3060_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3060_Raw, gr_PMT_S20_3060_Raw_Cfg},
             {gr_PMT_S20_3060_Filtered, gr_PMT_S20_3060_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3070
  GraphPainterConfig gr_PMT_S20_3070_Raw_Cfg;
  gr_PMT_S20_3070_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3070_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3070_Raw.SetNameTitle("gr_PMT_S20_3070", "gr_PMT_S20_3070");
  GraphPainterConfig gr_PMT_S20_3070_Filtered_Cfg;
  gr_PMT_S20_3070_Filtered_Cfg.option = "P";
  gr_PMT_S20_3070_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3070_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3070_Raw, gr_PMT_S20_3070_Raw_Cfg},
             {gr_PMT_S20_3070_Filtered, gr_PMT_S20_3070_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3179
  GraphPainterConfig gr_PMT_S20_3179_Raw_Cfg;
  gr_PMT_S20_3179_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3179_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3179_Raw.SetNameTitle("gr_PMT_S20_3179", "gr_PMT_S20_3179");
  GraphPainterConfig gr_PMT_S20_3179_Filtered_Cfg;
  gr_PMT_S20_3179_Filtered_Cfg.option = "P";
  gr_PMT_S20_3179_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3179_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3179_Raw, gr_PMT_S20_3179_Raw_Cfg},
             {gr_PMT_S20_3179_Filtered, gr_PMT_S20_3179_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3350
  GraphPainterConfig gr_PMT_S20_3350_Raw_Cfg;
  gr_PMT_S20_3350_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3350_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3350_Raw.SetNameTitle("gr_PMT_S20_3350", "gr_PMT_S20_3350");
  GraphPainterConfig gr_PMT_S20_3350_Filtered_Cfg;
  gr_PMT_S20_3350_Filtered_Cfg.option = "P";
  gr_PMT_S20_3350_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3350_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3350_Raw, gr_PMT_S20_3350_Raw_Cfg},
             {gr_PMT_S20_3350_Filtered, gr_PMT_S20_3350_Filtered_Cfg}},
            legCfg);

  // PMT_S20_3360
  GraphPainterConfig gr_PMT_S20_3360_Raw_Cfg;
  gr_PMT_S20_3360_Raw_Cfg.legTitle = "Raw";
  gr_PMT_S20_3360_Raw_Cfg.printTitle = true;
  gr_PMT_S20_3360_Raw.SetNameTitle("gr_PMT_S20_3360", "gr_PMT_S20_3360");
  GraphPainterConfig gr_PMT_S20_3360_Filtered_Cfg;
  gr_PMT_S20_3360_Filtered_Cfg.option = "P";
  gr_PMT_S20_3360_Filtered_Cfg.color = kRed;
  gr_PMT_S20_3360_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_PMT_S20_3360_Raw, gr_PMT_S20_3360_Raw_Cfg},
             {gr_PMT_S20_3360_Filtered, gr_PMT_S20_3360_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_1988_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_1988_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_1988_TMIT_PC",
                                             "gr_TORO_LI20_1988_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_1988_TMIT_PC_Raw, gr_TORO_LI20_1988_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_1988_TMIT_PC_Filtered,
              gr_TORO_LI20_1988_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_2040_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_2040_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_2040_TMIT_PC",
                                             "gr_TORO_LI20_2040_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_2040_TMIT_PC_Raw, gr_TORO_LI20_2040_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2040_TMIT_PC_Filtered,
              gr_TORO_LI20_2040_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_2452_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_2452_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_2452_TMIT_PC",
                                             "gr_TORO_LI20_2452_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_2452_TMIT_PC_Raw, gr_TORO_LI20_2452_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_2452_TMIT_PC_Filtered,
              gr_TORO_LI20_2452_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_3163_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_3163_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_3163_TMIT_PC",
                                             "gr_TORO_LI20_3163_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_3163_TMIT_PC_Raw, gr_TORO_LI20_3163_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_3163_TMIT_PC_Filtered,
              gr_TORO_LI20_3163_TMIT_PC_Filtered_Cfg}},
            legCfg);

  // TORO_LI20_3255_TMIT_PC
  GraphPainterConfig gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg;
  gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg.legTitle = "Raw";
  gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg.printTitle = true;
  gr_TORO_LI20_3255_TMIT_PC_Raw.SetNameTitle("gr_TORO_LI20_3255_TMIT_PC",
                                             "gr_TORO_LI20_3255_TMIT_PC");
  GraphPainterConfig gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg;
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.option = "P";
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.color = kRed;
  gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg.legTitle = "Filtered";
  saveGraph(c, cOption,
            {{gr_TORO_LI20_3255_TMIT_PC_Raw, gr_TORO_LI20_3255_TMIT_PC_Raw_Cfg},
             {gr_TORO_LI20_3255_TMIT_PC_Filtered,
              gr_TORO_LI20_3255_TMIT_PC_Filtered_Cfg}},
            legCfg);

  GraphPainterConfig grOccRawCfg;
  grOccRawCfg.xTitle = "EventId";
  grOccRawCfg.yTitle = "Clusters";
  grOccRawCfg.legTitle = "Raw";
  grOccRawCfg.printTitle = true;
  GraphPainterConfig grOccFiltered;
  grOccFiltered.option = "P";
  grOccFiltered.color = kRed;
  grOccFiltered.legTitle = "Filtered";
  for (int id : chipIds) {
    auto &grRaw = chipOccRawGrs.at(id);
    std::string grName = runPrefix + "ChipOccRaw" + std::to_string(id);
    grRaw.SetNameTitle(grName.c_str(), grName.c_str());
    saveGraph(
        c, cOption,
        {{grRaw, grOccRawCfg}, {chipOccFilteredGrs.at(id), grOccFiltered}},
        legCfg);
  }

  HistPainterConfig chipHitMapCfg;
  chipHitMapCfg.xTitle = "Pix X";
  chipHitMapCfg.yTitle = "Pix Y";
  for (int id : chipIds) {
    saveHist(c, cOption, {{chipHitMapRawHists.at(id), chipHitMapCfg}});
    saveHist(c, cOption, {{chipHitMapFilteredHists.at(id), chipHitMapCfg}});
    saveHist(c, cOption, {{chipHitMapROIRawHists.at(id), chipHitMapCfg}});
  }

  HistPainterConfig dummyCfg;
  TH1D dummy("dummy", "dummy", 100, 0, 100);
  saveHist(c, cOptionLast, {{dummy, dummyCfg}});
}
