#ifndef UTILS_H
#define UTILS_H

// std
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
// ROOT
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TArrow.h>
#include <TF1.h>
#include <TLatex.h>


// Input
class InputParser {
  public:
    InputParser (int&, char**);
    const std::string& getCmdOption(const std::string&) const;
    bool cmdOptionExists(const std::string&) const;
  private:
    std::vector<std::string> tokens;
};
// https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
std::vector<std::string> splitString(const std::string& s, char delimiter);
// https://stackoverflow.com/questions/1120140
class csvRow {
  public:
    std::string const& operator[](std::size_t) const;
    std::size_t size() const;
    void readNextRow(std::istream&);

  private:
    std::vector<std::string> m_data;
};
std::istream& operator>>(std::istream& str, csvRow& data);

int GetVecIndex(const std::vector<double> &binCenterVec, double binStart, double binStop);

// Time
tm GetCurrentTime();
std::string MonthNumToName(int);
std::string GetCurrentYear();
std::string GetCurrentMonth();
std::string GetCurrentDay();
std::string GetCurrentHour();
std::string GetCurrentMinute();
std::string GetCurrentSecond();
std::string Now();
std::string Today();
std::string Yesterday();


// Plotter
class Plotter {
  private:
    std::vector<TH1D*> hist;
    std::vector<TGraphAsymmErrors*> graph;
    std::vector<TArrow*> arrow;
    std::vector<TF1*> func;
    std::vector<TLatex*> textVec;

    std::vector<std::string> histDrawParams;
    std::vector<std::string> graphDrawParams;
    std::vector<std::string> funcDrawParams;

    std::vector<double> histLineWidth;
    std::vector<double> graphLineWidth;
    std::vector<double> arrowLineWidth;
    std::vector<double> funcLineWidth;

    std::vector<std::string> note;
    std::vector<int> color;

    std::vector<int> marker;

    std::string fileName;
    float tickLength;

  public:
    float yMin;
    float yMax;
    float xMin;
    float xMax;
    int logX;
    int logY;
    float legendX1;
    float legendX2;
    float legendY1;
    float legendY2;
    bool drawLegend;
    float atlasLabelX1;
    float atlasLabelX2;
    float atlasLabelY1;
    float atlasLabelY2;
    bool drawAtlasLabel;

    Plotter(const std::string&);
    void draw();
    void addHistogram(TH1D*);
    void addGraph(TGraphAsymmErrors*);
    void addGraph(TGraphAsymmErrors*, int);
    void addArrow(TArrow*);
    void addArrow(TArrow*, int);
    void addFunc(TF1*);
    void addText(TLatex* text);
    void addNote(const std::string&);
    void addNotes(const std::vector<std::string>&);
    TH1D* getHist(int);
    TGraphAsymmErrors* getGraph(int);
    TF1* getFunc(int);
    void setHistDrawParam(int, const std::string&);
    void setGraphDrawParam(int, const std::string&);
    void setFuncDrawParam(int, const std::string&);
    void setHistLineWidth(int, double);
    void setGraphLineWidth(int, double);
    void setFuncLineWidth(int, double);
};



#endif /* UTILS_H */
