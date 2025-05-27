#ifndef PLOT_H
#define PLOT_H

//std
#include <vector>
#include <string>
// ROOT
#include <TGraphAsymmErrors.h>
#include <TArrow.h>
#include <TLatex.h>

bool help;
bool verbose;
std::vector<std::string> paramNameVec;
std::vector<std::string> paramNameSlugVec;
std::string plotName;
std::string c2_type;
std::string c2_number;
std::string ref_type;
std::vector<std::string> ref_type_now;
std::vector<std::string> fileVec;
std::vector<TGraphAsymmErrors*> graphVec;
std::vector<TArrow*> arrowVec;
std::vector<std::string> titleVec;
std::vector<std::string> commentVec;
std::vector<TLatex*> textVec;
int ref_num;
double max_binStop;
double yMin;
double yMax;
double legX[2];
double legY[2];
std::string xLabel;
std::string yLabel;
bool logY;
bool plotMult;
bool plotKt;
bool plottt;


void LoadHistograms();
void Plot();

#endif /* PLOT_H */
