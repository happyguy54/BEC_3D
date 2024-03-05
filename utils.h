#ifndef UTILS_H
#define UTILS_H

// std
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
// ROOT
#include <TGraph.h>
#include <TH3.h>
#include <TF3.h>
#include <TPaveText.h>
#include "Fit/BinData.h"


extern int projMin;
extern int projMax;
extern int projRange;
extern double zMin;
extern double zMax;
extern bool showC0;
extern bool fixC0;
extern bool useEps;

extern TString c2IndexText;
extern TString qOutText;
extern TString qSideText;
extern TString qLongText;


using std::vector;
using std::string;
using std::find;
using std::to_string;


// Projections
TGraph* makeProjection(TH3D*,TF3*, const string&, double, double, int, int);
TH1D* makeProjectionHist(ROOT::Fit::BinData&, const string&, double, int, int);
void plot2Dprojection(TH3D*, TF3*, const std::string&,
                      std::vector<std::string>,
                      const TString&,
                      int, int);

// Histogram manipulation
void enlargeUncertainties(TH3D*);
void scaleAxisMeVtoGeV(TH2D*, std::string);
void scaleAxesMeVtoGeV(TH2D*);

// Fit results
void addFitResults(TPaveText*, TF3*);
//void outputCSV(TH3D*, const std::string&);
void outputCSV(TF3*, const std::string&);

// Input
class InputParser {
  public:
    InputParser (int&, char**);
    const string& getCmdOption(const string&) const;
    bool cmdOptionExists(const string&) const;
  private:
    vector<string> tokens;
};
// https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
std::vector<std::string> splitString(const std::string& s, char delimiter);

// Time
tm GetCurrentTime();
string MonthNumToName(int);
string GetCurrentYear();
string GetCurrentMonth();
string GetCurrentDay();
string GetCurrentHour();
string GetCurrentMinute();
string GetCurrentSecond();
string Now();
string Today();

#endif /* UTILS_H */
