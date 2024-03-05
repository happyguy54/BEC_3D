// std
#include <iostream>
#include <sstream>
// ROOT
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
// Fit
#include "utils.h"


// Input
InputParser::InputParser (int& argc, char** argv) {
  for (int i = 1; i < argc; ++i)
    this->tokens.push_back(std::string(argv[i]));
}

const std::string& InputParser::getCmdOption(const std::string& option) const {
  std::vector<std::string>::const_iterator itr;
  itr = find(this->tokens.begin(), this->tokens.end(), option);
  if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
    return *itr;
  }
  static const std::string empty_string("");

  return empty_string;
}

bool InputParser::cmdOptionExists(const std::string& option) const {
  return find(this->tokens.begin(),
              this->tokens.end(),
              option) != this->tokens.end();
}

// https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
std::vector<std::string> splitString(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.emplace_back(token);
  }

  return tokens;
}

// https://stackoverflow.com/questions/1120140
std::string const& csvRow::operator[](std::size_t index) const {
  return m_data[index];
}

std::size_t csvRow::size() const {
  return m_data.size();
}

void csvRow::readNextRow(std::istream& str) {
  std::string line;
  std::getline(str, line);

  std::stringstream lineStream(line);
  std::string cell;

  m_data.clear();
  while (std::getline(lineStream, cell, ',')) {
    m_data.push_back(cell);
  }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty()) {
    // If there was a trailing comma then add an empty element.
    m_data.push_back("");
  }
}

std::istream& operator>>(std::istream& str, csvRow& data) {
  data.readNextRow(str);

  return str;
}

int GetVecIndex(const std::vector<double> &binCenterVec,
                double binStart, double binStop) {
  int binCenterIndex = 0;
  for (size_t i = 0; i < binCenterVec.size(); ++i) {
    if (binStart <= binCenterVec.at(i) && binCenterVec.at(i) <= binStop) {
      binCenterIndex = i;
    }
  }

  return binCenterIndex;
}

// Time
tm GetCurrentTime() {
  typedef std::chrono::duration<int,
      std::ratio_multiply<std::chrono::hours::period,
          std::ratio<24> >::type> days;

  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
  std::chrono::system_clock::duration tp = now.time_since_epoch();
  days d = std::chrono::duration_cast<days>(tp);
  tp -= d;
  std::chrono::hours h = std::chrono::duration_cast<std::chrono::hours>(tp);
  tp -= h;
  std::chrono::minutes m = std::chrono::duration_cast<std::chrono::minutes>(tp);
  tp -= m;
  std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(tp);
  tp -= s;

  time_t tt = std::chrono::system_clock::to_time_t(now);
  tm local_tm = *localtime(&tt);

  return local_tm;
}

std::string MonthNumToName(int month) {
  std::vector<std::string> months;
  months.emplace_back("Jan");
  months.emplace_back("Feb");
  months.emplace_back("Mar");
  months.emplace_back("Apr");
  months.emplace_back("May");
  months.emplace_back("Jun");
  months.emplace_back("Jul");
  months.emplace_back("Aug");
  months.emplace_back("Sep");
  months.emplace_back("Oct");
  months.emplace_back("Nov");
  months.emplace_back("Dec");

  return months.at(month);
}

std::string GetCurrentYear() {
  tm now = GetCurrentTime();
  int year = now.tm_year + 1900;

  return std::to_string(year);
}

std::string GetCurrentMonth() {
  tm now = GetCurrentTime();
  int month = now.tm_mon;

  return MonthNumToName(month);
}

std::string GetCurrentDay() {
  tm now = GetCurrentTime();
  int day = now.tm_mday;

  std::string sDay;
  if (day < 10) {
    sDay += "0";
  }
  sDay += std::to_string(day);

  return sDay;
}

std::string GetCurrentHour() {
  tm now = GetCurrentTime();
  int hour = now.tm_hour;

  std::string sHour;
  if (hour < 10) {
    sHour += "0";
  }
  sHour += std::to_string(hour);

  return sHour;
}

std::string GetCurrentMinute() {
  tm now = GetCurrentTime();
  int minute = now.tm_min;

  std::string sMinute;
  if (minute < 10) {
    sMinute += "0";
  }
  sMinute += std::to_string(minute);

  return sMinute;
}

std::string GetCurrentSecond() {
  tm now = GetCurrentTime();
  int second = now.tm_sec;

  std::string sSecond;
  if (second < 10) {
    sSecond += "0";
  }
  sSecond += std::to_string(second);

  return sSecond;
}

std::string Now() {
  std::string now = GetCurrentHour();
  now += ":";
  now += GetCurrentMinute();
  now += ":";
  now += GetCurrentSecond();
  now += " ";
  now += GetCurrentDay();
  now += " ";
  now += GetCurrentMonth();
  now += " ";
  now += GetCurrentYear();

  return now;
}

std::string Today() {
  std::string today = GetCurrentYear();
  today += "_";
  today += GetCurrentMonth();
  today += "_";
  today += GetCurrentDay();

  return today;
}

// Plotter
Plotter::Plotter(const std::string& inFileName) {
  // Setup
  fileName = inFileName;
  yMin = -1e7;
  yMax = -1e7;
  xMin = -1e7;
  xMax = -1e7;
  logX = 0;
  logY = 0;
  legendX1 = .450;
  legendX2 = .98;
  legendY1 = .5;
  legendY2 = .92;
  drawLegend = true;
  atlasLabelX1 = .25;
  atlasLabelX2 = .55;
  atlasLabelY1 = .15;
  atlasLabelY2 = .25;
  drawAtlasLabel = true;

//  color.emplace_back(kBlack);
  color.emplace_back(kBlue);
  color.emplace_back(kRed+1);
  color.emplace_back(kGreen+3);
  color.emplace_back(kViolet-6);
  color.emplace_back(kCyan-0);
  color.emplace_back(kPink-3);
  color.emplace_back(kYellow-4);
  color.emplace_back(kOrange+4);


  marker.emplace_back(20);
  marker.emplace_back(21);
  marker.emplace_back(22);
  marker.emplace_back(23);
  marker.emplace_back(24);
  marker.emplace_back(25);
  marker.emplace_back(26);
}

void Plotter::draw() {
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 350, 350);
  gPad->SetTopMargin(.05);
  gPad->SetLeftMargin(.1);
  gPad->SetBottomMargin(.1);
  gPad->SetRightMargin(.05);

  TLegend *legend;
  if (drawLegend) {
    legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetShadowColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(12);
  }

  TPaveText* atlasLabel = new TPaveText(atlasLabelX1, atlasLabelY1,
                                        atlasLabelX2, atlasLabelY2, "NDC");
  atlasLabel->SetFillStyle(0);
  atlasLabel->SetFillColor(0);
  atlasLabel->SetShadowColor(0);
  atlasLabel->SetBorderSize(0);
  atlasLabel->SetTextFont(43);
  atlasLabel->SetTextSize(14);
  atlasLabel->AddText("#bf{#it{ATLAS} Internal #sqrt{s} = 13 TeV}");

  // Error: gPad->SetLogx() and gPad->SetLogy() erase graph axis titles
  std::vector<std::string> xtitles, ytitles;
  std::string xtitle, ytitle;
  for (size_t i = 0; i < graph.size(); ++i) {
    xtitle = graph.at(i)->GetXaxis()->GetTitle();
    ytitle = graph.at(i)->GetYaxis()->GetTitle();
    xtitles.emplace_back(xtitle);
    ytitles.emplace_back(ytitle);
  }
  gStyle->SetOptStat(0);
  gPad->SetLogx(logX);
  gPad->SetLogy(logY);
  for (size_t i = 0; i < graph.size(); ++i) {
    graph.at(i)->GetXaxis()->SetTitle(xtitles.at(i).c_str());
    graph.at(i)->GetYaxis()->SetTitle(ytitles.at(i).c_str());
  }
  xtitles.clear();
  ytitles.clear();

  int xMinBin = 1;
  int xMaxBin = 1;
  if (xMin > -1e6 && xMax > -1e6) {
    for (size_t i = 0; i < hist.size(); ++i) {
      xMinBin = hist.at(i)->FindBin(xMin) - 1;
      if (xMinBin < 1) {
        xMinBin = 1;
      }
      xMaxBin = hist.at(i)->FindBin(xMax) + 1;
      if (xMaxBin > hist.at(i)->GetNbinsX()) {
        xMaxBin = hist.at(i)->GetNbinsX();
      }
      hist.at(i)->GetXaxis()->SetRange(xMinBin, xMaxBin);
    }
    for (size_t i = 0; i < graph.size(); ++i) {
      graph.at(i)->GetXaxis()->SetRange(xMin, xMax);
    }
  }

  float histMax = 1.;
  float histMin = .5;
  for (size_t i = 0; i < hist.size(); ++i) {
    if (hist.at(i)->GetMaximum() > histMax) {
      histMax = hist.at(i)->GetMaximum();
    }
    if (hist.at(i)->GetMinimum() < histMin) {
      histMin = hist.at(i)->GetMinimum();
    }
  }
  for (size_t i = 0; i < graph.size(); ++i) {
    if (TMath::MaxElement(graph.at(i)->GetN(), graph.at(i)->GetY()) > histMax) {
      histMax = TMath::MaxElement(graph.at(i)->GetN(), graph.at(i)->GetY());
    }
    if (TMath::MinElement(graph.at(i)->GetN(), graph.at(i)->GetY()) < histMin) {
      histMin = TMath::MinElement(graph.at(i)->GetN(), graph.at(i)->GetY());
    }
  }
  if (logY) {
    histMin = histMin * .5;
    histMax = histMax * 2;
  }
  else {
    histMin = histMin * .9;
    histMax = histMax * 1.1;
  }
  if (yMax > -1e6) histMax = yMax;
  if (yMin > -1e6) histMin = yMin;

  // TGaxis::SetMaxDigits(2);

  int nDraw = 0;
  for (size_t i = 0; i < hist.size(); ++i) {
    hist.at(i)->SetLineWidth(histLineWidth.at(i));
    hist.at(i)->SetMarkerSize(.5);

    hist.at(i)->GetXaxis()->SetLabelFont(43);
    hist.at(i)->GetXaxis()->SetLabelSize(12);
    hist.at(i)->GetXaxis()->SetTitleFont(43);
    hist.at(i)->GetXaxis()->SetTitleSize(12);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.3);

    hist.at(i)->GetYaxis()->SetLabelFont(43);
    hist.at(i)->GetYaxis()->SetLabelSize(12);
    hist.at(i)->GetYaxis()->SetTitleFont(43);
    hist.at(i)->GetYaxis()->SetTitleSize(12);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.2);

    hist.at(i)->SetMinimum(histMin);
    hist.at(i)->SetMaximum(histMax);

    if (drawLegend) legend->AddEntry(hist.at(i), hist.at(i)->GetTitle(),
                                     histDrawParams.at(i).c_str());
    hist.at(i)->SetTitle("");

    if (nDraw == 0) hist.at(i)->Draw(histDrawParams.at(i).c_str());
    else hist.at(i)->Draw((histDrawParams.at(i) + "same").c_str());
    ++nDraw;
  }

  for (size_t i = 0; i < graph.size(); ++i) {
    graph.at(i)->SetLineWidth(graphLineWidth.at(i));
    graph.at(i)->SetMarkerSize(.5);

    graph.at(i)->GetXaxis()->SetLabelFont(43);
    graph.at(i)->GetXaxis()->SetLabelSize(12);
    graph.at(i)->GetXaxis()->SetTitleFont(43);
    graph.at(i)->GetXaxis()->SetTitleSize(12);
    graph.at(i)->GetXaxis()->SetTitleOffset(1.3);

    graph.at(i)->GetYaxis()->SetLabelFont(43);
    graph.at(i)->GetYaxis()->SetLabelSize(12);
    graph.at(i)->GetYaxis()->SetTitleFont(43);
    graph.at(i)->GetYaxis()->SetTitleSize(12);
    graph.at(i)->GetYaxis()->SetTitleOffset(1.2);

    graph.at(i)->SetMinimum(histMin);
    graph.at(i)->SetMaximum(histMax);

    if (drawLegend) {
      legend->AddEntry(graph.at(i),
                       graph.at(i)->GetTitle(),
                       graphDrawParams.at(i).c_str());
    }
    graph.at(i)->SetTitle("");

    if (nDraw == 0) {
      graph.at(i)->Draw((graphDrawParams.at(i) + "A").c_str());
    }
    else {
      graph.at(i)->Draw((graphDrawParams.at(i) + "same").c_str());
    }
    ++nDraw;
  }

  for (size_t i = 0; i < func.size(); ++i) {
    func.at(i)->SetLineWidth(funcLineWidth.at(i));
    func.at(i)->SetMarkerSize(.5);

    func.at(i)->GetXaxis()->SetLabelFont(43);
    func.at(i)->GetXaxis()->SetLabelSize(12);
    func.at(i)->GetXaxis()->SetTitleFont(43);
    func.at(i)->GetXaxis()->SetTitleSize(12);
    func.at(i)->GetXaxis()->SetTitleOffset(1.3);

    func.at(i)->GetYaxis()->SetLabelFont(43);
    func.at(i)->GetYaxis()->SetLabelSize(12);
    func.at(i)->GetYaxis()->SetTitleFont(43);
    func.at(i)->GetYaxis()->SetTitleSize(12);
    func.at(i)->GetYaxis()->SetTitleOffset(1.2);

    if (drawLegend) legend->AddEntry(func.at(i), func.at(i)->GetTitle(),
                                     funcDrawParams.at(i).c_str());
    func.at(i)->SetTitle("");

    func.at(i)->Draw(funcDrawParams.at(i).c_str());
    if (nDraw == 0) func.at(i)->Draw(funcDrawParams.at(i).c_str());
    else func.at(i)->Draw((funcDrawParams.at(i) + "same").c_str());
    ++nDraw;
  }

  for (size_t i = 0; i < arrow.size(); ++i) {
    arrow.at(i)->SetLineWidth(arrowLineWidth.at(i));
    arrow.at(i)->Draw();
  }

  for (size_t i = 0; i < note.size(); ++i) {
    legend->AddEntry((TObject*)0, note.at(i).c_str(), "");
  }

  if (drawLegend) legend->Draw();
  if (drawAtlasLabel) atlasLabel->Draw();

  canvas->Update();
  canvas->Print((fileName + ".pdf").c_str());

  delete canvas;
  if (drawLegend) delete legend;
  if (atlasLabel) delete atlasLabel;
  hist.clear();
  graph.clear();
  func.clear();
  arrow.clear();

  return;
}


void Plotter::addHistogram(TH1D* inHist) {
  std::string copy = "_copy";
  TH1D* histHelper = (TH1D*) inHist->Clone((inHist->GetName() + copy).c_str());

  std::string drawParamsHelper = "LE1P";

  int nObj = hist.size() + graph.size();
  histHelper->SetLineColor(color.at(nObj % color.size()));
  histHelper->SetMarkerColor(color.at(nObj % color.size()));
  histHelper->SetMarkerStyle(marker.at(nObj % marker.size()));

  hist.emplace_back(histHelper);
  histDrawParams.emplace_back(drawParamsHelper);
  histLineWidth.emplace_back(2.);

  return;
}

void Plotter::addGraph(TGraphAsymmErrors* inGraph) {
  std::string copy = "_copyafsdfs";
  TGraphAsymmErrors* graphHelper = dynamic_cast<TGraphAsymmErrors*>(
      inGraph->Clone((inGraph->GetName() + copy).c_str()));

  std::string drawParamsHelper = "EP";

  int nObj = hist.size() + graph.size();
  graphHelper->SetLineColor(color.at(nObj % color.size()));
  graphHelper->SetMarkerColor(color.at(nObj % color.size()));
  graphHelper->SetMarkerStyle(marker.at(nObj % marker.size()));

  graph.emplace_back(graphHelper);
  graphDrawParams.emplace_back(drawParamsHelper);
  graphLineWidth.emplace_back(2);

  return;
}

void Plotter::addArrow(TArrow* inArrow) {
  std::string copy = "_copyarrow";
  TArrow* arrowHelper = dynamic_cast<TArrow*>(
      inArrow->Clone((inArrow->GetName() + copy).c_str()));
  int nObj = hist.size() + arrow.size();
  arrowHelper->SetLineColor(color.at(nObj % color.size()));
  arrow.emplace_back(arrowHelper);
  arrowLineWidth.emplace_back(2);

  return;
}

void Plotter::addFunc(TF1* inFunc) {
  std::string copy = "_copy";
  TF1* funcHelper = (TF1*) inFunc->Clone((inFunc->GetName() + copy).c_str());


  std::string drawParamsHelper = "L";

  int nObj = hist.size() + func.size();
  funcHelper->SetLineColor(color.at(nObj % color.size()));
  funcHelper->SetMarkerColor(color.at(nObj % color.size()));
  funcHelper->SetMarkerStyle(marker.at(nObj % marker.size()));

  func.emplace_back(funcHelper);
  funcDrawParams.emplace_back(drawParamsHelper);
  funcLineWidth.emplace_back(2);

  return;
}

void Plotter::addNote(const std::string& inNote) {
  note.emplace_back(inNote);

  return;
}

void Plotter::addNotes(const std::vector<std::string>& inNotes) {
  for (size_t i = 0; i < inNotes.size(); ++i) {
    note.emplace_back(inNotes.at(i));
  }

  return;
}

TH1D* Plotter::getHist(int index) {

  return hist.at(index);
}


TGraphAsymmErrors* Plotter::getGraph(int index) {

  return graph.at(index);
}


TF1* Plotter::getFunc(int index) {

  return func.at(index);
}


void Plotter::setHistDrawParam(int index, const std::string& param) {
  histDrawParams.at(index) = param;

  return;
}


void Plotter::setGraphDrawParam(int index, const std::string& param) {
  graphDrawParams.at(index) = param;

  return;
}


void Plotter::setFuncDrawParam(int index, const std::string& param) {
  funcDrawParams.at(index) = param;

  return;
}


void Plotter::setHistLineWidth(int index, double width) {
  histLineWidth.at(index) = width;

  return;
}


void Plotter::setGraphLineWidth(int index, double width) {
  graphLineWidth.at(index) = width;

  return;
}


void Plotter::setFuncLineWidth(int index, double width) {
  funcLineWidth.at(index) = width;

  return;
}
