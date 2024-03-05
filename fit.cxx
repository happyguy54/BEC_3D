// std
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <ios>
#include <memory>
#include <set>
#include <algorithm>
// ROOT
#include <TFile.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TH2.h>
#include <TArrow.h>
#include <TLatex.h>
// Fit
#include "fit.h"
#include "c2.h"
#include "utils.h"
#include "TError.h"
#include "Fit/BinData.h"
#include <Fit/DataOptions.h>
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include <Math/WrappedParamFunction.h>
#include <TMatrixDSym.h>


// std
using std::cout;
using std::endl;
using std::string;
using std::vector;

bool help = false;
bool verbose = false;//true;//
double cubes = true;

string fileData =
    "input/bec_alln_W_trPt_100_0_Qcut_0_Qbin_20_100_0_3D_2_GeV.root";
string file2Data = "";
string fileMC =
    "input/bec_mc_alln_W_sc_trPt_100_0_Qcut_0_Qbin_20_100_0_3D_2_GeV.root";
string file2MC = "";

string histData = "ppmm_qosl_g";
string hist2Data = "";
string histMC = "";
string hist2MC = "";
bool r2 = false;

int c2Index = 1;
double qMin = -1.;
double qMax = -1.;
double rejFrom = 1.;
double rejTo = -1.;
vector<double> rejFromOsl = {1., 1., 1.};
vector<double> rejToOsl = {-1., -1., -1.};
// rej[0,1,2] is equivalent to rej[out,side,long]
double rej2From = 1.;
double rej2To = -1.;
bool useEps = false;
bool showC0 = false;
bool fixC0 = false;
bool enlargeUncert = false;
int nParams = 0;

string fitParams = "RE";

string plotTitle = "FIT 3D";
string plotName = "bec_3d_fit";
int sign_chi = 1;
double yMin = .8;
double yMax = 1.3;
double zMin = 0.89;
double zMax = 1.25;
double ar_y;
TLatex lexc;
int projMin = 1;
int projMax = 12;
int projRange = 1;
size_t nSamples = 50;
double qStep = 30;
vector<double> qBoud;
double projQmin;
double projQmax;
double projQrange;
TString c2IndexText;
TString qOutText;
TString qSideText;
TString qLongText;
string comment = "";
vector<string> commentVec;
ROOT::Fit::BinData DataStruc(100000, 3);
ROOT::Fit::BinData DataStruc_0(100000, 3);
int npoints_DataStruc;
std::vector<std::unique_ptr<TH3D>> InHistRatio;
std::vector<double> qMin3D_0 = {.0, .28, .6, 1.2, 4.2};
std::vector<double> qMin3D_02 = {0.02, .28, .6, 1.2, 4.2};
// std::vector<TH3D*> InHistRatio;

std::vector<double> varQbins = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
                                0.22, 0.24, 0.26, 0.28, 0.30, 0.35, 0.4, 0.45, 0.5, 0.6,
                                0.7, 0.8, 1.0, 1.2, 1.5, 2.0};
/*{0., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200.,
                     220., 240., 260., 280., 300., 360., 420., 480., 540., 720.,
900., 1080., 1620., 2160., 2700., 3240., 3780.};*/
int varQbinsSize = varQbins.size() - 1;

std::string outFilePathBase = "output/";

int main(int argc, char** argv) {
  InputParser input(argc, argv);
  if (input.cmdOptionExists("-h") ||
      input.cmdOptionExists("--help")) {
    help = true;
  }
  if (input.cmdOptionExists("-v") ||
      input.cmdOptionExists("--verbose")) {
    verbose = true;
  }

  if (input.cmdOptionExists("--file-data")) {
    fileData = input.getCmdOption("--file-data");
  }
  if (input.cmdOptionExists("--file2-data")) {
    file2Data = input.getCmdOption("--file2-data");
  }
  if (input.cmdOptionExists("--file-mc")) {
    fileMC = input.getCmdOption("--file-mc");
  }
  if (input.cmdOptionExists("--file2-mc")) {
    file2MC = input.getCmdOption("--file2-mc");
  }

  if (input.cmdOptionExists("--hist-data")) {
    histData = input.getCmdOption("--hist-data");
  }
  if (input.cmdOptionExists("--hist2-data")) {
    hist2Data = input.getCmdOption("--hist2-data");
  }
  if (input.cmdOptionExists("--hist-mc")) {
    histMC = input.getCmdOption("--hist-mc");
  }
  if (input.cmdOptionExists("--hist2-mc")) {
    hist2MC = input.getCmdOption("--hist2-mc");
  }

  if (input.cmdOptionExists("-i")) {
    string param = input.getCmdOption("-i");
    c2Index = std::stoi(param);
  }
  if (input.cmdOptionExists("--c2-index")) {
    string param = input.getCmdOption("--c2-index");
    c2Index = std::stoi(param);
  }
  if (input.cmdOptionExists("--q-min")) {
    string param = input.getCmdOption("--q-min");
    qMin = std::stod(param);
  }
  if (input.cmdOptionExists("--q-max")) {
    string param = input.getCmdOption("--q-max");
    qMax = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-from")) {
    string param = input.getCmdOption("--rej-from");
    rejFrom = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-from-out")) {
    string param = input.getCmdOption("--rej-from-out");
    rejFromOsl[0] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-from-side")) {
    string param = input.getCmdOption("--rej-from-side");
    rejFromOsl[1] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-from-long")) {
    string param = input.getCmdOption("--rej-from-long");
    rejFromOsl[2] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-to")) {
    string param = input.getCmdOption("--rej-to");
    rejTo = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-to-out")) {
    string param = input.getCmdOption("--rej-to-out");
    rejToOsl[0] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-to-side")) {
    string param = input.getCmdOption("--rej-to-side");
    rejToOsl[1] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej-to-long")) {
    string param = input.getCmdOption("--rej-to-long");
    rejToOsl[2] = std::stod(param);
  }
  if (input.cmdOptionExists("--rej2-from")) {
    string param = input.getCmdOption("--rej2-from");
    rej2From = std::stod(param);
  }
  if (input.cmdOptionExists("--rej2-to")) {
    string param = input.getCmdOption("--rej2-to");
    rej2To = std::stod(param);
  }
  if (input.cmdOptionExists("--use-eps")) {
    useEps = true;
  }
  if (input.cmdOptionExists("--show-c0")) {
    showC0 = true;
  }
  if (input.cmdOptionExists("--fix-c0")) {
    fixC0 = true;
  }
  if (input.cmdOptionExists("--enlarge-uncert")) {
    enlargeUncert = true;
  }

  if (input.cmdOptionExists("--fit-params")) {
    fitParams = input.getCmdOption("--fit-params");
  }

  if (input.cmdOptionExists("--plot-title")) {
    plotTitle = input.getCmdOption("--plot-title");
  }
  if (input.cmdOptionExists("--plot-name")) {
    plotName = input.getCmdOption("--plot-name");
  }
  if (input.cmdOptionExists("--y-min")) {
    string param = input.getCmdOption("--y-min");
    yMin = std::stod(param);
  }
  if (input.cmdOptionExists("--y-max")) {
    string param = input.getCmdOption("--y-max");
    yMax = std::stod(param);
  }
  if (input.cmdOptionExists("--proj-min")) {
    string param = input.getCmdOption("--proj-min");
    projMin = std::stoi(param);
  }
  if (input.cmdOptionExists("--proj-max")) {
    string param = input.getCmdOption("--proj-max");
    projMax = std::stoi(param);
  }
  if (input.cmdOptionExists("--num-samples")) {
    string param = input.getCmdOption("--num-samples");
    nSamples = std::stoi(param);
  }
  if (input.cmdOptionExists("--comment")) {
    comment = input.getCmdOption("--comment");
  }

  if (help) {
    cout << "INFO: Usage: " << endl;
  }

  c2IndexText.Form("C_{2}^{fit}: %i", c2Index);
  commentVec = splitString(comment, ';');

  // qMin = -1;
  if (qMin > 0) {
    for(auto& value : varQbins) {
      value += qMin;
    }
    varQbins.insert(varQbins.begin(), 0.0);
    cout << "Bins for 1d histograms changed" << endl;
    varQbinsSize = varQbins.size() - 1;
  }

  try {
    MakeHistogram();
    Fit();
    Output();
    Plot();
  } catch (const char* msg) {
    cout << msg << endl;

    return 1;
  }

  return 0;
}

bool rej_OSL() {
  std::vector<int> diff(rejToOsl.size());
  std::transform(rejToOsl.begin(), rejToOsl.end(), rejFromOsl.begin(), diff.begin(), std::minus<int>());
  return std::all_of(diff.begin(), diff.end(), [](int i){ return i > 0; });
}

void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, ROOT::Fit::BinData& data, const std::vector<double>& qMin3D = {0., .28, .6, 1.2, 4.2});
// void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, vector<std::unique_ptr<TH3D>>& InHistRatiotest, ROOT::Fit::BinData& data) {
void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, ROOT::Fit::BinData& data, const std::vector<double>& qMin3D) {
  const int hist_number = qMin3D.size()-1;
  InHistRatio.clear();
  double ix_down_limit = 0, iy_down_limit = 0, iz_down_limit = 0;
  for (int i=0;i<4-hist_number;i++) {
    // push empty histograms to InHistRatio
    InHistRatio.push_back(std::make_unique<TH3D>());
  }
  for(int i=(4-hist_number);i<4;i++) {
    // auto inHistRatio = std::make_unique<TH3D>(*HistData[i]);
    std::unique_ptr<TH3D> inHistRatio;
    //both vectors cant be empty, the if in which this function is called checks for that
    if (!HistData.empty()) {
      inHistRatio = std::make_unique<TH3D>(*HistData[i]);
      // cout<<"Size of HistData: "<<inHistRatio->GetSize()<<endl;
      if (!HistMC.empty()) {
        inHistRatio->Divide(HistData[i].get(), HistMC[i].get());
           //                    1./HistData[i]->Integral(), 1./HistMC[i]->Integral());
        cout<<"Fitting R2 function"<<endl;
      } else 
        cout<<"Fitting C2 function for data"<<endl;
    } else if (!HistMC.empty()) {
      inHistRatio = std::make_unique<TH3D>(*HistMC[i]);
      cout<<"Fitting C2 function for MC"<<endl;
    } else {
      throw "ERROR: Can't make histogram!";
    }
    inHistRatio->SetDirectory(0);
    // std::vector<double> qMin3D = {0., .28, .6, 1.2, 4.2};
    if (qMin3D[0] > 0.06)
      throw "ERROR: lower limit must be less than 0.06!"; //otherwise we do not fit the BEC peak
    // std::vector<double> qMin3D = {lower_limit, .28, .6, 1.2, 4.2};
    if (i==(4-hist_number)) {
      ix_down_limit = inHistRatio->GetXaxis()->FindBin(qMin3D[0]);
      cout<<"ix_down_limit "<<ix_down_limit<<endl;
      iy_down_limit = inHistRatio->GetYaxis()->FindBin(qMin3D[0]);
      iz_down_limit = inHistRatio->GetZaxis()->FindBin(qMin3D[0]);
    }
    int ixStart = inHistRatio->GetXaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int iyStart = inHistRatio->GetYaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int izStart = inHistRatio->GetZaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    for (int ix=ix_down_limit; ix<=inHistRatio->GetNbinsX(); ++ix) {
      for (int iy=iy_down_limit; iy<=inHistRatio->GetNbinsY(); ++iy) {
        for (int iz=iz_down_limit; iz<=inHistRatio->GetNbinsZ(); ++iz) {
          if ((ix<ixStart)&&(iy<iyStart)&&(iz<izStart))
            continue;
          double binContent = inHistRatio->GetBinContent(ix, iy, iz);
          double binError = inHistRatio->GetBinError(ix, iy, iz);
          double binCoord[3] ={inHistRatio->GetXaxis()->GetBinCenter(ix),
                                inHistRatio->GetYaxis()->GetBinCenter(iy),
                                inHistRatio->GetZaxis()->GetBinCenter(iz)};
          data.Add(binCoord, binContent, binError); //binCoordErr,
          // cout<<inHistRatio->GetXaxis()->GetBinCenter(ix)<<"  "<<inHistRatio->GetYaxis()->GetBinCenter(iy)<<"  "<<inHistRatio->GetZaxis()->GetBinCenter(iz)<<"\t"<<binError<<endl;
        }
      }
    }
    InHistRatio.push_back(std::move(inHistRatio));
    // InHistRatio.push_back(inHistRatio.get());
    // check InHistRatio content
    //cout size of InHistRatio
    // cout<<"InHistRatio size: "<<InHistRatio[i]->GetSize()<<endl;
    // for (unsigned int n = 0; n < InHistRatio[i]->GetSize(); n++) {
    //   std::cout << "InHistRatio content: " << InHistRatio[i]->GetBinContent(n) << std::endl;
    // }
    //if (inHistRatio) delete inHistRatio;
  }
   std::cout << "Number of points in data object: " << data.Size() << std::endl;
}

void draw_hist_chi(std::vector<std::unique_ptr<TH3D>>& Hist, const std::string& suffix) {
    int iter = 1;
    std::string plotNewNameAll = plotName+"_chi_all_"+suffix;
    auto c_combined = std::make_unique<TCanvas>(plotNewNameAll.c_str(), plotNewNameAll.c_str(), 800, 600);

    for (auto rit = Hist.rbegin(); rit != Hist.rend(); ++rit, ++iter) {
        auto& histtemp = *rit;
        // histtemp->SetDirectory(0);
        std::string plotNewName = plotName+"_chi_"+suffix+"_"+std::to_string(5-iter);
        std::make_unique<TCanvas>(plotNewName.c_str(), plotNewName.c_str(), 800, 600);
        auto c_individual = std::make_unique<TCanvas>(plotNewName.c_str(), plotNewName.c_str(), 800, 600);
        
        histtemp->SetMarkerStyle(20);
        histtemp->SetMarkerSize(1);
        histtemp->SetMarkerColor(kRed);
        histtemp->SetLineColor(kRed);
        histtemp->SetLineWidth(1);
        histtemp->SetFillColorAlpha(kRed, 0.5);
        histtemp->SetFillStyle(3002);
        histtemp->SetMinimum(0.01);
        histtemp->SetMaximum(0.04);
        // histtemp->SetMinimum(avgChi);
        // histtemp->SetMaximum(avgChi*6);
        // Draw on individual canvas and write it
        c_individual->cd();
        histtemp->Draw("BOX2Z");
        c_individual->Write(plotNewName.c_str());

        // Draw on combined canvas
        c_combined->cd();
        if (rit == Hist.rbegin()) {
            histtemp->Draw("BOX2Z");
        } else {
            histtemp->Draw("BOX2Z SAME");
        }
    }

    // Write the combined canvas
    c_combined->Write(plotNewNameAll.c_str());
}

// void DataToHist(ROOT::Fit::BinData& data, vector<TH3D*>& HistData, vector<TH3D*>& InHistRatio) {
//     std::vector<double> qMin3D = {0., .28, .6, 1.2, 4.2};

//     for(int i=0; i<4; i++) {
//         inHistRatio = new TH3D(*HistData[i]); // Assuming HistData has the correct binning

//         int ixStart = HistData[i]->GetXaxis()->FindBin(qMin3D[i]);
//         int iyStart = HistData[i]->GetYaxis()->FindBin(qMin3D[i]);
//         int izStart = HistData[i]->GetZaxis()->FindBin(qMin3D[i]);

//         for (int ix=qMin3D[0]+1; ix<=HistData[i]->GetNbinsX(); ++ix) {
//             for (int iy=qMin3D[0]+1; iy<=HistData[i]->GetNbinsY(); ++iy) {
//                 for (int iz=qMin3D[0]+1; iz<=HistData[i]->GetNbinsZ(); ++iz) {
//                     if ((ix<ixStart)&&(iy<iyStart)&&(iz<izStart))
//                         continue;

//                     double binCoord[3] ={inHistRatio->GetXaxis()->GetBinCenter(ix), 
//                                           inHistRatio->GetYaxis()->GetBinCenter(iy),
//                                           inHistRatio->GetZaxis()->GetBinCenter(iz)};

//                     unsigned int bin = data.FindBin(binCoord);
//                     double value, inv_error, error;
//                     const double* coords = data.GetPoint(bin, value, inv_error);
//                     error = 1./inv_error;

//                     inHistRatio->SetBinContent(ix, iy, iz, value);
//                     inHistRatio->SetBinError(ix, iy, iz, error);
//                 }
//             }
//         }

//         InHistRatio->push_back(inHistRatio); // Assuming InHistRatio is the output
//     }
// }


void MakeHistogram() {
    cout << "INFO: Creating histogram..." << endl;
    std::unique_ptr<TH3D> inHistData;
    std::unique_ptr<TH3D> inHistMC;
    std::vector<std::unique_ptr<TH3D>> inHistRatio;
    string histtemp;
    std::vector<std::unique_ptr<TH3D>> inHistVec;
    std::vector<std::unique_ptr<TH3D>> inHistDataVec;
    std::vector<std::unique_ptr<TH3D>> inHistMCVec;
    string sign = "ppmm";
    string sign2 = "pm";
    std::vector<double> qMin3D_4_0 = {.0, .28, .6, 1.2, 4.2};
    std::vector<double> qMin3D_4_02 = {0.02, .28, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_0 = {.0, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_04 = {0.04, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_02 = {0.02, .62, 1.22, 4.22};

    qMin3D_0 = qMin3D_3_0;
    qMin3D_02 = qMin3D_3_02;

    if (!fileData.empty() && !histData.empty()) {
      std::unique_ptr<TFile> inFile = std::make_unique<TFile>(fileData.c_str(), "READ");
      if (!inFile->IsOpen()) {
        throw "ERROR: Data file not found!";
      }
      std::unique_ptr<TH3D> inHist(dynamic_cast<TH3D*>(inFile->Get(histData.c_str())));
      if (!inHist) {
        throw "ERROR: Data histogram not found!";
      }
      inHist->SetDirectory(0);
      for (int i = 1; i <= 4; i++) {
        histtemp = histData;
        size_t pos = histtemp.find(sign);
        histtemp.insert(pos, to_string(i) + "_");
        std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
        inHisttemp->SetDirectory(0);
        inHistVec.push_back(std::move(inHisttemp));
      }

      if (!hist2Data.empty()) {
        std::unique_ptr<TH3D> inHist2;
        std::vector<std::unique_ptr<TH3D>> inHistVec2;
        if (file2Data.empty()) {
          inHist2 = std::make_unique<TH3D>(*dynamic_cast<TH3D*>(inFile->Get(hist2Data.c_str())));
          if (!inHist2) {
            throw "ERROR: Second data histogram not found!";
          }
          for (int i = 1; i <= 4; i++) {
            histtemp = hist2Data;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
        } else {
          std::unique_ptr<TFile> inFile2 = std::make_unique<TFile>(file2Data.c_str(), "READ");
          if (!inFile2->IsOpen()) {
            throw "ERROR: Second data file not found!";
          }
          inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile2->Get(hist2Data.c_str())));
          if (!inHist2) {
            throw "ERROR: Second data histogram not found!";
          }
          inHist2->SetDirectory(0);
          for (int i = 1; i <= 4; i++) {
            histtemp = hist2Data;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile2->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          inFile2->Close();
        }
        for (int i = 0; i < 4; i++) {
          if (!inHistVec2[i]) {
            throw "ERROR: Second data histogram not found!";
          }
          inHistVec[i]->Divide(inHistVec[i].get(), inHistVec2[i].get(),
                              1./inHist->Integral(), 1./inHist2->Integral());

          string histName = inHistVec[i]->GetName();
          histName += "_over_";
          histName += inHistVec2[i]->GetName();
          inHistVec[i]->SetName(histName.c_str());
        }
        
        inHist->Divide(inHist.get(), inHist2.get(),
                       1./inHist->Integral(), 1./inHist2->Integral());

        string histName = inHist->GetName();
        histName += "_over_";
        histName += inHist2->GetName();
        inHist->SetName(histName.c_str());
        inHist2.reset();
        inHistVec2.clear();
      }

      // inHistData = std::make_unique<TH3D>(*inHist);
      inHistData = std::move(inHist);
      if (verbose) inHistData->Print();
      for (auto& inHisttemp : inHistVec) {
        inHistDataVec.push_back(std::move(inHisttemp));
      }
      inHistVec.clear();
      inFile->Close();
    }

    if (!fileMC.empty() && !histMC.empty()) {
      std::unique_ptr<TFile> inFile = std::make_unique<TFile>(fileMC.c_str(), "READ");
      if (!inFile->IsOpen()) {
        throw "ERROR: MC file not found!";
      }
      std::unique_ptr<TH3D> inHist(dynamic_cast<TH3D*>(inFile->Get(histMC.c_str())));
      if (!inHist) {
        throw "ERROR: MC histogram not found!";
      }
      inHist->SetDirectory(0);
      for (int i = 1; i <= 4; i++) {
        histtemp = histMC;
        size_t pos = histtemp.find(sign);
        histtemp.insert(pos, to_string(i) + "_");
        std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
        inHisttemp->SetDirectory(0);
        inHistVec.push_back(std::move(inHisttemp));
      }
      if (!hist2MC.empty()) {
        std::unique_ptr<TH3D> inHist2;
        std::vector<std::unique_ptr<TH3D>> inHistVec2;
        if (file2MC.empty()) {
          inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile->Get(hist2MC.c_str())));
          if (!inHist2) {
            throw "ERROR: Second data histogram not found!";
          }
          for (int i = 1; i <= 4; i++) {
            histtemp = hist2MC;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
        } else {
          std::unique_ptr<TFile> inFile2 = std::make_unique<TFile>(file2MC.c_str(), "READ");
          if (!inFile2->IsOpen()) {
            throw "ERROR: Second data file not found!";
          }
          inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile2->Get(hist2MC.c_str())));
          if (!inHist2) {
            throw "ERROR: Second data histogram not found!";
          }
          inHist2->SetDirectory(0);
          for (int i = 1; i <= 4; i++) {
            histtemp = hist2MC;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile2->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          inFile2->Close();
        }
        for (int i = 0; i < 4; i++) {
          inHistVec[i]->Divide(inHistVec[i].get(), inHistVec2[i].get(),
                              1./inHist->Integral(), 1./inHist2->Integral());

          string histName = inHistVec[i]->GetName();
          histName += "_over_";
          histName += inHistVec2[i]->GetName();
          inHistVec[i]->SetName(histName.c_str());
        }
        inHistVec2.clear();
        inHist->Divide(inHist.get(), inHist2.get(),
                       1./inHist->Integral(), 1./inHist2->Integral());

        string histName = inHist->GetName();
        histName += "_over_";
        histName += inHist2->GetName();
        inHist->SetName(histName.c_str());
        inHist2.reset();
      }

      // inHistMC = std::make_unique<TH3D>(*inHist);
      inHistMC = std::move(inHist);
      if (verbose) inHistMC->Print();
      for (auto& inHisttemp : inHistVec) {
        inHistMCVec.push_back(std::move(inHisttemp));
      }
      inHistVec.clear();
      inFile->Close();
    }

    if (inHistData && inHistMC) {
      string histName = inHistData->GetName();
      histName += "_over_";
      histName += inHistMC->GetName();
      inHistRatio.push_back(std::make_unique<TH3D>(*inHistData));
      if (!hist2Data.empty() || !hist2MC.empty()) {
        inHistRatio[0]->Divide(inHistData.get(), inHistMC.get());
      } else {
        inHistRatio[0]->Divide(inHistData.get(), inHistMC.get(),
                            1./inHistData->Integral(), 1./inHistMC->Integral());
      }
      // inHistData.reset();
      // inHistMC.reset();
      // hist = inHistRatio[0].get();
      hist = new TH3D(*(inHistRatio[0].get()));
      //clear InHistRatio
      // inHistRatio.clear();
      // HisttoData(inHistDataVec, inHistMCVec, inHistRatio, DataStruc);
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      std::cout << "DataStruc size: " << DataStruc.Size() << std::endl;
      
      r2 = true;
    } else if (!inHistData && inHistMC) {
      hist = new TH3D(*(inHistMC.get()));
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      r2 = false;
    } else if (inHistData && !inHistMC) {
      hist = new TH3D(*(inHistData.get()));
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      r2 = false;
    } else {
      throw "ERROR: Can't make histogram!";
    }

    if (verbose) {
      cout << "INFO: Histogram:" << endl;
      hist->Print();
    }

    if (enlargeUncert) {
      enlargeUncertainties(hist);
    }

    if (qMin < 0.) {
      qMin = hist->GetXaxis()->GetXmin();
    }
    if (qMax < 0.) {
      qMax = hist->GetXaxis()->GetXmax();
    }
    qStep = (qMax - qMin) / nSamples;
    std::cout<<"qmax "<<qMax<<" qmin "<<qMin<<" qstep "<<qStep<<std::endl;
    for (size_t i = 0; i < nSamples; ++i) {
      qBoud.emplace_back(qMin + i * qStep);
    }
    qBoud.emplace_back(qMax);
    projRange = projMax - projMin + 1;
    projQmin = 1000*hist->GetXaxis()->GetBinLowEdge(projMin);
    projQmax = 1000*(hist->GetXaxis()->GetBinLowEdge(projMax) +
                hist->GetXaxis()->GetBinWidth(projMax));
    projQrange = projQmax - projQmin;
    std::cout<<"projqmax "<<projQmax/1000<<" projqmin "<<projQmin/1000<<std::endl;

    qOutText.Form("%.0f #leq Q_{out} < %.0f MeV", projQmin, projQmax);
    qSideText.Form("%.0f #leq Q_{side} < %.0f MeV", projQmin, projQmax);
    qLongText.Form("%.0f #leq Q_{long} < %.0f MeV", projQmin, projQmax);
}
  /*
  void Fit() {
  cout << "INFO: Fitting..." << endl;

  C2* c2;
  switch (c2Index) {
    case 1: c2 = new C2_1();
            break;
    case 2: c2 = new C2_2();
            break;
    case 3: c2 = new C2_3();
            break;
    case 4: c2 = new C2_4();
            break;    
    default: throw "ERROR: C2 function index not found!";
             break;
  }

  c2->rejFrom = rejFrom;
  c2->rejTo = rejTo;
  c2->rej2From = rej2From;
  c2->rej2To = rej2To;
  c2->fixC0 = fixC0;
  c2->useEps = useEps;

  func = new TF3(("C2_" + std::to_string(c2Index)).c_str(), c2,
                 qMin, qMax,
                 qMin, qMax,
                 qMin, qMax,
                 c2->nParams);
  c2->setup(func);

  if (verbose) {
    cout << "INFO: Fit function:" << endl;
    func->Print();
    cout << "INFO: Q range" << endl;
    cout << "      from: " << qMin << endl;
    cout << "      to: " << qMax << endl;
    cout << "INFO: Rejection region" << endl;
    cout << "      from: " << c2->rejFrom << endl;
    cout << "      to: " << c2->rejTo << endl;
    cout << "INFO: Second rejection region " << endl;
    cout << "      from: " << c2->rej2From << endl;
    cout << "      to: " << c2->rej2To << endl;
    cout << "INFO: Use epsilon: " << c2->useEps << endl;

    fitParams += "V";
  }

  hist->Fit(func, fitParams.c_str());
  gMinuit->mnmatu(1);
}
*/


void Fit() {
  cout << "INFO: Fitting..." << endl;
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
  // minuit2->SetPrintLevel(3);
  // ROOT::Minuit2::Minuit2Minimizer *minuit2 = new ROOT::Minuit2::Minuit2Minimizer();
  C2* c2;
  switch (c2Index) {
    case 1: c2 = new C2_1();
            break;
    case 2: c2 = new C2_2();
            break;
    case 3: c2 = new C2_3();
            break;
    case 4: c2 = new C2_4();
            break;
    case 5: c2 = new C2_5();
            break;
    case 6: c2 = new C2_6();
            break;
    default: throw "ERROR: C2 function index not found!";
             break;
  }

  c2->rejFrom = rejFrom;
  c2->rejFromOsl = rejFromOsl;
  c2->rejTo = rejTo;
  c2->rejToOsl = rejToOsl;
  c2->rej2From = rej2From;
  c2->rej2To = rej2To;
  c2->fixC0 = fixC0;
  c2->useEps = useEps;
  nParams = c2->nParams;
  cout<<"nParams "<<nParams<<endl;
  
  func = new TF3(("C2_" + std::to_string(c2Index)).c_str(), c2,
                 qMin, qMax,
                 qMin, qMax,
                 qMin, qMax,
                 c2->nParams);
  c2->setup(func);

  ROOT::Fit::Fitter fitterMigrad;
  ROOT::Math::WrappedMultiTF1 wf(*func, 3);
  fitterMigrad.SetFunction(wf);
  fitterMigrad.Config().MinimizerOptions().SetMinimizerAlgorithm("Migrad");

  ROOT::Fit::Fitter fitterSimulatedAnnealing;
  ROOT::Math::WrappedMultiTF1 wf2(*func, 3);
  fitterSimulatedAnnealing.SetFunction(wf2);

  fitterSimulatedAnnealing.Config().MinimizerOptions().SetMinimizerAlgorithm("SimulatedAnnealing");
  // fitterSimulatedAnnealing.Config().MinimizerOptions().SetMaxIterations(10000);
  // fitterSimulatedAnnealing.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
  fitterMigrad.Config().MinimizerOptions().SetMaxIterations(10000);
  fitterMigrad.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
  // fitterMigrad.Config().MinimizerOptions().SetTolerance(3000000.0);
  // fitterSimulatedAnnealing.Config().MinimizerOptions().SetTolerance(3000000.0);
  fitterSimulatedAnnealing.Config().MinimizerOptions().SetPrintLevel(2);
  for (int i = 0; i < nParams; i++) {
    double min, max;
    func->GetParLimits(i, min, max);
    fitterMigrad.Config().ParSettings(i).SetLimits(min, max);
    fitterSimulatedAnnealing.Config().ParSettings(i).SetLimits(min, max);
  }

  bool retSimulatedAnnealing = fitterSimulatedAnnealing.Fit(DataStruc);
  const ROOT::Fit::FitResult & resSimulatedAnnealing = fitterSimulatedAnnealing.Result();
  bool retMigrad = fitterMigrad.Fit(DataStruc);
  const ROOT::Fit::FitResult & resMigrad = fitterMigrad.Result();

  bool ret = false;
  const ROOT::Fit::FitResult *res;
  // Compare the results
  if (retMigrad || retSimulatedAnnealing) {
    std::cout<<"Migrad chi2 "<<resMigrad.Chi2()<<" Simulated Annealing chi2 "<<resSimulatedAnnealing.Chi2()<<std::endl;
    if (resMigrad.Chi2() == resSimulatedAnnealing.Chi2()) {
      std::cout << "Both fits are equally good" << std::endl;
      ret = retMigrad;
      res = &resMigrad;
    } else {
      std::cout<<"THERE IS DIFFERENCE IN CHI2: "<<resMigrad.Chi2()-resSimulatedAnnealing.Chi2()<<std::endl;
      if (resMigrad.Chi2() < resSimulatedAnnealing.Chi2()) {
        std::cout << "Migrad fit is better" << std::endl;
        ret = retMigrad;
        res = &resMigrad;
      } else {
        std::cout << "Simulated Annealing fit is better" << std::endl;
        ret = retSimulatedAnnealing;
        res = &resSimulatedAnnealing;
      }
    }
  // if (retMigrad) {
  //   std::cout << "Migrad fit is successfil" << std::endl;
  //   ret = retMigrad;
  //   res = &resMigrad;
  } else {
    std::cout << "One or both fits failed" << std::endl;
  }

  if (ret) {
    //const ROOT::Fit::FitResult & res = fitter.Result();  THIS WAS WORKING BEFORE!!!
    // print result
    //res.Print(std::cout); THIS WAS WORKING BEFORE!!!
    res->Print(std::cout);
    // copy all fit result info (values, chi2, etc..) in TF3
    // func->SetFitResult(res); THIS WAS WORKING BEFORE!!!
    func->SetFitResult(*res);
    // test fit p-value (chi2 probability)
    // double prob = res.Prob(); THIS WAS WORKING BEFORE!!!
    double prob = res->Prob();
    if (prob < 1.E-2)
      Error("Fit", "Bad data fit - fit p-value is %f", prob);
    else
      std::cout << "Good fit : p-value  = " << prob << std::endl;

    gMinuit->mnmatu(1);
    // TFitResult t_res (res);



    //matrix print
    // TMatrixDSym correlationMatrix(res->NPar());
    // res->GetCorrelationMatrix(correlationMatrix);
    // // //skip parameters with fixed values
    // std::vector<int> fixedIndices;
    // std::vector<std::string> parameterNames;
    // for (int i = 0; i < res->NPar(); ++i) {
    //   if (fitterMigrad.Config().ParSettings(i).IsFixed()) {
    //     fixedIndices.push_back(i);
    //   } else 
    //     parameterNames.push_back(fitterMigrad.Config().ParSettings(i).Name());
    // }

    // TMatrixDSym tempMatrix(res->NPar() - fixedIndices.size());
    // int ii = 0;
    // for (int row = 0; row < correlationMatrix.GetNrows(); ++row) {
    //   if (std::find(fixedIndices.begin(), fixedIndices.end(), row) != fixedIndices.end()) continue;
    //   int jj = 0;
    //   for (int col = 0; col < correlationMatrix.GetNcols(); ++col) {
    //     if (std::find(fixedIndices.begin(), fixedIndices.end(), col) != fixedIndices.end()) continue;
    //     tempMatrix(ii, jj) = correlationMatrix(row, col);
    //     ++jj;
    //   }
    //   ++ii;
    // }
    // correlationMatrix.ResizeTo(tempMatrix.GetNrows(), tempMatrix.GetNcols());
    // correlationMatrix = tempMatrix;
    // for (int i = 0; i < correlationMatrix.GetNcols(); ++i) {
    //   std::cout << "\t" << parameterNames[i];
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < correlationMatrix.GetNrows(); ++i) {
    //   std::cout << parameterNames[i] << "\t";
    //   for (int j = 0; j < correlationMatrix.GetNcols(); ++j) {
    //       std::cout << correlationMatrix(i, j) << "\t";
    //   }
    //   std::cout << std::endl;
    // }
    // correlationMatrix.Print();
  }
  else
    Error("Fit", "3D fit failed");
}


void Output() {
  cout << "INFO: Saving..." << endl;

  outFilePathBase = "output/" + Today() + "/c2_" + std::to_string(c2Index);
  system(("mkdir -p " + outFilePathBase).c_str());

  /**
   * ROOT files
   */
  string outFilePathRoot = outFilePathBase + "/root";
  system(("mkdir -p " + outFilePathRoot).c_str());
  outFilePathRoot += "/";
  outFilePathRoot += plotName;
  outFilePathRoot += ".root";

  if (verbose) {
    cout << "INFO: Saving ROOT file here:" << endl;
    cout << "      " << outFilePathRoot << endl;
  }

  TFile* outFile = new TFile(outFilePathRoot.c_str(), "RECREATE");
  cout<<"outFilePathRoot "<<outFilePathRoot<<"plotName "<<plotName<<endl;
  hist->Write(plotName.c_str());
  func->Write(("fitFunc_" + plotName).c_str());

  std::map<std::tuple<double, double, double>, int> coordinates;

  for (unsigned int i = 0; i < DataStruc.Size(); ++i) {
      double value;
      const double* coords = DataStruc.GetPoint(i, value);

      coordinates[std::make_tuple(coords[0], coords[1], coords[2])] = i;
      // cout<<coords[0]<<"  "<<coords[1]<<"  "<<coords[2]<<"bin "<<i<<endl;
  }
  // copy hist to new hist
  std::vector<double> qMin3D = qMin3D_02;
  std::vector<int> qNBins3D = {14, 10, 7, 0};

  // int nBins = 0;
  // // int projRan = 0;
  // vector<double> binVal;
  // vector<double> binErr;
  // vector<double> newEdges;
  // binVal.push_back(0);
  // binErr.push_back(0);
  // newEdges.push_back(0);
  // for (int i = 0; i < qMin3D.size()-1; i++) {
  //   int bins = qNBins3D[i]*((qMin3D[i+1]-qMin3D[i])/qMin3D[i+1]);
  //   nBins += bins;
  //   double binWidth = qMin3D[i+1] / qNBins3D[i];
  //   for (int j = 0; j < bins; j++) {
  //     newEdges.push_back(qMin3D[i] + j*binWidth);
  //     //std::cout<<qMin3D[i] + j*binWidth;
  //     binVal.push_back(0);
  //     binErr.push_back(0);
  //   }
  // }
  int Q_values = 100;
  // std::unique_ptr<TH1D> hist_Q_1D(new TH1D("hist_Q_1D", "hist_Q_1D", nBins, newEdges.data()));
  std::unique_ptr<TH1D> hist_Q_1D(new TH1D("hist_Q_1D", "hist_Q_1D", Q_values, 0, 4.2));
  hist_Q_1D->SetDirectory(0);

  std::vector<double> chi_values(Q_values);
  std::vector<double> top_Q(Q_values);
  double peak_value_all=0;
  double peak_fit_all=0;
  double bins_peak_all=0;
  // find chi quadrat between data and fit for every bin, create new datastructure
  ROOT::Fit::BinData DataStruc_chi_p(100000, 3);
  ROOT::Fit::BinData DataStruc_chi_m(100000, 3);
  double avgChi = 0;
  double maxChi = 0;
  double chi = 0;
  sign_chi = 1;
  for (unsigned int i = 0; i < DataStruc.Size(); ++i) {
      double value;
      double error, inv_error;

      // Get the coordinates and value for the bin
      // const double* coords = DataStruc.GetPoint(i, value, inv_error);
      const double* coords = DataStruc.GetPoint(i, value);
      error = DataStruc.Error(i);
      //const double* coordsErr = DataStruc.GetPointError(i, error);
      // error = 1./inv_error;

      // Calculate the expected value
      double expected = func->Eval(coords[0], coords[1], coords[2]);
      // double chi = sqrt(pow(value - expected, 2));

      double min_chi = *std::min_element(chi_values.begin(), chi_values.end());
      int min_index = std::distance(chi_values.begin(), std::find(chi_values.begin(), chi_values.end(), min_chi));

      chi = sqrt(pow(value - expected, 2));

      if (chi > min_chi) {
          chi_values[min_index] = chi;
          top_Q[min_index] = 0;
          for (int j = 0; j < 3; ++j) {
              top_Q[min_index] += pow(coords[j], 2);
          }
          top_Q[min_index] = sqrt(top_Q[min_index]);
      }
      // cout<<"chi "<<chi<<" value "<<value<<" expected "<<expected<<" coords "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<endl;
      if ((value - expected) > 0) {
        DataStruc_chi_p.Add(coords, chi);
        DataStruc_chi_m.Add(coords, 0);
      } else if ((value - expected) < 0) {
        DataStruc_chi_m.Add(coords, chi);
        DataStruc_chi_p.Add(coords, 0);
      }
      // Add the chi-squared value to the new BinData object
      avgChi += chi;
      if (chi > maxChi) {
        maxChi = chi;
      }
  }
  // fill hist_Q_1D with chi values
  for (int i = 0; i < Q_values; ++i) {
      hist_Q_1D->Fill(top_Q[i]);
  }
  // Draw the histogram
  // hist_Q_1D->Draw("BOX2Z");
  hist_Q_1D->Write(("hist_Q_1D_" + plotName).c_str());
  // cout<<"peak_value_all "<<peak_value_all<<" peak_fit_all "<<peak_fit_all<<" bins_peak_all "<<bins_peak_all<<" ratio "<<(peak_value_all-peak_fit_all)/bins_peak_all<<endl;
  avgChi /= DataStruc.Size();
  cout<<"avgChi "<<avgChi<<" maxChi "<<maxChi<<endl;
  // scale average chi
  // avgChi *= 3;
  avgChi *= 1;
  std::vector<std::unique_ptr<TH3D>> Hist_p;
  std::vector<std::unique_ptr<TH3D>> Hist_m;

  for (auto& histtemp : InHistRatio) {
      std::unique_ptr<TH3D> newHist;
      std::unique_ptr<TH3D> newHist2;
  //     try {
  //       histtemp->GetName();  // Or any other TH3D method
  //   } catch (...) {
  //       std::cerr << "Error: Invalid TH3D object in InHistRatio\n";
  //       continue;
  //   }
  // // for (auto& histtemp : InHistRatio) {
  //   if (histtemp == nullptr) {
  //       std::cerr << "Error: histtemp is a null pointer\n";
  //       continue;
  //   }
    newHist = std::make_unique<TH3D>(*histtemp);
    newHist2 = std::make_unique<TH3D>(*histtemp);
    newHist->SetDirectory(0);
    newHist2->SetDirectory(0);
    // histtemp.reset();
    // newHist = std::move(histtemp);
    // check newHist size
    cout<<"newHist size "<<newHist->GetSize()<<endl;
    newHist->Reset();
    newHist2->Reset();
    Hist_p.push_back(std::move(newHist));
    Hist_m.push_back(std::move(newHist2));
    cout<<"Hist_p size "<<Hist_p.size()<<" Hist_m size "<<Hist_m.size()<<endl;
  }
  // check InHistRatio contents size
  // for (auto& histtemp : InHistRatio) {
  //   cout<<"InHistRatio size "<<histtemp->GetSize()<<endl;
  // }
  // check InHistRatio[i] size
  // for(int i=0;i<4;i++) {
  //   cout<<"InHistRatio[i] size "<<InHistRatio[i]->GetSize()<<endl;
  // }
  // Get the number of points in the BinData object
  // Create a 3D vector with the same dimensions


  double integral_peak = 0;
  double bins_peak = 0;
  // std::vector<int> qNBins3D = {14, 15, 10, 7, 0};
  // std::vector<double> qWidth3D = {0.02, 0.04, 0.12, 0.6};
  const int hist_number = qMin3D.size()-1;
  cout<<"hist_number "<<hist_number<<endl;
  for(int i=4-hist_number;i<4;i++) {
    //
    int ixStart = InHistRatio[i]->GetXaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int iyStart = InHistRatio[i]->GetYaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int izStart = InHistRatio[i]->GetZaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    // cout<<"ix starting from "<<InHistRatio[4-hist_number]->GetXaxis()->FindBin(qMin3D[0])<<" total bins "<<InHistRatio[i]->GetNbinsX()<<endl;
    for (int ix=InHistRatio[4-hist_number]->GetXaxis()->FindBin(qMin3D[0]); ix<=InHistRatio[i]->GetNbinsX(); ++ix) {
      for (int iy=InHistRatio[4-hist_number]->GetYaxis()->FindBin(qMin3D[0]); iy<=InHistRatio[i]->GetNbinsY(); ++iy) {
        for (int iz=InHistRatio[4-hist_number]->GetZaxis()->FindBin(qMin3D[0]); iz<=InHistRatio[i]->GetNbinsZ(); ++iz) {
          if ((ix<ixStart)&&(iy<iyStart)&&(iz<izStart)) {
            // InHistRatio[i]->SetBinContent(ix, iy, iz, 0);
            // InHistRatio[i]->SetBinError(ix, iy, iz, 0);
            Hist_p[i]->SetBinContent(ix, iy, iz, 0);
            Hist_p[i]->SetBinError(ix, iy, iz, 0);
            Hist_m[i]->SetBinContent(ix, iy, iz, 0);
            Hist_m[i]->SetBinError(ix, iy, iz, 0);
            continue;
          }
            
          // double binContent = inHistRatio->GetBinContent(ix, iy, iz);
          // double binError = inHistRatio->GetBinError(ix, iy, iz); 
          double binCoord[3] ={InHistRatio[i]->GetXaxis()->GetBinCenter(ix), 
                                InHistRatio[i]->GetYaxis()->GetBinCenter(iy),
                                InHistRatio[i]->GetZaxis()->GetBinCenter(iz)};
          //unsigned int bin = DataStruc_chi.FindBin(binCoord);
          // int bin = InHistRatio[i]->FindBin(binCoord[0], binCoord[1], binCoord[2]);
          int bin = coordinates[std::make_tuple(binCoord[0], binCoord[1], binCoord[2])];
          double value_p, value_m, inv_error, error;
          // get bin number for coordinates
          // int bin_test = coordinates[std::make_tuple(-2, -2, -2)];
          // const double* coords_test = DataStruc_chi.GetPoint(bin_test, value, inv_error);
          // cout<<"value "<<value<<"bin "<<bin_test<<"<- coords "<<coords_test[0]<<"  "<<coords_test[1]<<"  "<<coords_test[2]<<endl;
          const double* coords_p = DataStruc_chi_p.GetPoint(bin, value_p, inv_error);
          const double* coords_m = DataStruc_chi_m.GetPoint(bin, value_m, inv_error);
          //error = 1./inv_error;
          //data.Add(binCoord, binContent, binError); //binCoordErr,
          // Hist[i]->SetBinContent(ix, iy, iz, value);
          // if only bins with chi 3x bigger than average chi
          float bin_scale=1;//pow(qWidth3D[i]/qWidth3D[0], 1.5);
          // InHistRatio[i]->SetBinContent(ix, iy, iz, value/bin_scale);
          if(value_p > (avgChi))
            Hist_p[i]->SetBinContent(ix, iy, iz, value_p/bin_scale);
          else
            Hist_p[i]->SetBinContent(ix, iy, iz, 0);
          if(value_m > (avgChi))
            Hist_m[i]->SetBinContent(ix, iy, iz, value_m/bin_scale);
          else
            Hist_m[i]->SetBinContent(ix, iy, iz, 0);
            //cout<<InHistRatio[i]->GetXaxis()->GetBinCenter(ix)<<"  "<<InHistRatio[i]->GetYaxis()->GetBinCenter(iy)<<"  "<<InHistRatio[i]->GetZaxis()->GetBinCenter(iz)<<"\t"<<value<<endl;
            // Hist[i]->SetBinError(ix, iy, iz, 0);
            
            //cout<<inHistRatio->GetXaxis()->GetBinCenter(ix)<<"  "<<inHistRatio->GetYaxis()->GetBinCenter(iy)<<"  "<<inHistRatio->GetZaxis()->GetBinCenter(iz)<<"\t"<<binError<<endl;
            
            // InHistRatio[i]->SetBinContent(ix, iy, iz, 0);
          // InHistRatio[i]->SetBinError(ix, iy, iz, 0);
          Hist_p[i]->SetBinError(ix, iy, iz, 0);
          Hist_m[i]->SetBinError(ix, iy, iz, 0);
        }
      }
    }
    // cout<<"integral_peak "<<integral_peak/bins_peak<<endl;
  }

  draw_hist_chi(Hist_p, "p");
  draw_hist_chi(Hist_m, "m");
  //define canvas

  //draw Hist wich each bin as cube with the cube size as chi value
  // for (auto& histtemp : Hist) {
  // loop in descending order: for (auto& histtemp : InHistRatio) {
  // std::string plotNewNameAll = plotName+"_chi_all";
  // auto c_combined = std::make_unique<TCanvas>(plotNewNameAll.c_str(), plotNewNameAll.c_str(), 800, 600);

  // int iter = 1;
  // // for (auto rit = InHistRatio.rbegin(); rit != InHistRatio.rend(); ++rit, ++iter) {
  // for (auto rit = Hist.rbegin(); rit != Hist.rend(); ++rit, ++iter) {
  //     auto& histtemp = *rit;
  //     histtemp->SetDirectory(0);
  //     std::string plotNewName = plotName+"_chi_"+std::to_string(5-iter);
  //     std::make_unique<TCanvas>(plotNewName.c_str(), plotNewName.c_str(), 800, 600);
  //     auto c_individual = std::make_unique<TCanvas>(plotNewName.c_str(), plotNewName.c_str(), 800, 600);
      
  //     histtemp->SetMarkerStyle(20);
  //     histtemp->SetMarkerSize(1);
  //     histtemp->SetMarkerColor(kRed);
  //     histtemp->SetLineColor(kRed);
  //     histtemp->SetLineWidth(1);
  //     histtemp->SetFillColorAlpha(kRed, 0.5);
  //     histtemp->SetFillStyle(3002);
  //     histtemp->SetMinimum(0.01);
  //     histtemp->SetMaximum(0.04);
  //     // histtemp->SetMinimum(avgChi);
  //     // histtemp->SetMaximum(avgChi*6);
  //     // Draw on individual canvas and write it
  //     c_individual->cd();
  //     histtemp->Draw("BOX2Z");
  //     c_individual->Write(plotNewName.c_str());

  //     // Draw on combined canvas
  //     c_combined->cd();
  //     // if (rit == InHistRatio.rbegin()) {
  //     if (rit == Hist.rbegin()) {
  //         histtemp->Draw("BOX2Z");
  //     } else {
  //         histtemp->Draw("BOX2Z SAME");
  //     }
  // }

  // // Write the combined canvas
  // c_combined->Write(plotNewNameAll.c_str());

  // InHistRatio.clear();
  //   int iter = 1;
  // for (auto rit = InHistRatio.rbegin(); rit != InHistRatio.rend(); ++rit) {
  //   auto& histtemp = *rit;
  //   std::string plotNewName = plotName+"_chi_"+std::to_string(iter);
  //   auto c1 = std::make_unique<TCanvas>(plotNewName.c_str(), plotNewName.c_str(), 800, 600);
  //   histtemp->SetMarkerStyle(20);
  //   histtemp->SetMarkerSize(1);
  //   histtemp->SetMarkerColor(kRed);
  //   histtemp->SetLineColor(kRed);
  //   histtemp->SetLineWidth(1);
  //   histtemp->SetFillColorAlpha(kRed, 0.5);
  //   histtemp->SetFillStyle(3002);
  //   histtemp->Draw("BOX2");
  //   // histtemp->SetMinimum(avgChi);
  //   // histtemp->SetMaximum(maxChi);
  //   c1->Write(plotNewName.c_str());
  //   iter++;
  // }
  // delete Hist vector



  // Hist.clear();
  std::size_t pos = plotName.find("c2");
  std::string firstPart = "";
  if (pos != std::string::npos) {
      firstPart = plotName.substr(0, pos + 4); // +4 to include "c2" and the next 2 character
  } else {
      std::cout << "Substring not found." << std::endl;
  }
  string outFilePathTxt = outFilePathBase + "/txt";
  system(("mkdir -p " + outFilePathTxt).c_str());
  outFilePathTxt += "/";
  outFilePathTxt += firstPart;
  outFilePathTxt += ".txt";
  std::ofstream outTxt;
  // open txt file, clear it and write the plot name

  outTxt.open(outFilePathTxt);
  outTxt << plotName << endl << endl;
  // outTxt.open(outFilePathTxt, std::ios::app);

  // cout only bins and their chi if its bigger than average chi
  for (unsigned int i = 0; i < DataStruc_chi_p.Size(); ++i) {
      double value_chi_p, value_chi_m, value;
      double error_chi, error, inv_error;

      const double* coords_chi_p = DataStruc_chi_p.GetPoint(i, value_chi_p, error_chi);
      const double* coords_chi_m = DataStruc_chi_m.GetPoint(i, value_chi_m, error_chi);
      const double* coords = DataStruc.GetPoint(i, value, inv_error);
      error = 1./inv_error;
      double expected = func->Eval(coords[0], coords[1], coords[2]);

      // Add the chi-squared value and its error to the new BinData object
      if (value_chi_p > (3*avgChi) or value_chi_m > (3*avgChi)) {
        // cout << "Bin " << i << " has chi = " << value_chi << " and avgChi = " << avgChi << endl;
        // cout << "x = " << coords[0] << ", y = " << coords[1] << ", z = " << coords[2] << endl;
        // cout << "expected = " << expected << ", value = " << value << endl << endl;
        //output the same cout in txt file
        
        outTxt << "Bin " << i << ": chi_above_fit = " << value_chi_p << "; chi_below_fit = " << value_chi_m << "; 3*avgChi = " << 3*avgChi;
        outTxt << "; x = " << coords[0] << ", y = " << coords[1] << ", z = " << coords[2];
        outTxt << "; expected = " << expected << ", value = " << value << ", error = " << error << endl << endl;
      }
  }
  outTxt << endl;
  outTxt.close();

  // for (int ix=1; ix<=hist->GetNbinsX(); ++ix) {
  //   for (int iy=1; iy<=hist->GetNbinsY(); ++iy) {
  //     for (int iz=1; iz<=hist->GetNbinsZ(); ++iz) {
  //       double binContent = hist->GetBinContent(ix, iy, iz);
  //       double binError = hist->GetBinError(ix, iy, iz);
  //       double binCoord[3] ={hist->GetXaxis()->GetBinCenter(ix), 
  //                             hist->GetYaxis()->GetBinCenter(iy),
  //                             hist->GetZaxis()->GetBinCenter(iz)};
  //       double binContent_fit = func->Eval(binCoord[0], binCoord[1], binCoord[2]);
  //       double binError_fit = 0;
  //       double binCoord_fit[3] ={hist->GetXaxis()->GetBinCenter(ix), 
  //                             hist->GetYaxis()->GetBinCenter(iy),
  //                             hist->GetZaxis()->GetBinCenter(iz)};
  //       double binContent_chi = (binContent - binContent_fit)*(binContent - binContent_fit)/(binError*binError)/(DataStruc.Size()-nParams);
  //       double binError_chi = 0;
  //       DataStruc_chi.Add(binCoord, binContent_chi, binError_chi); //binCoordErr,
  //       cout<<hist->GetXaxis()->GetBinCenter(ix)<<"  "<<hist->GetYaxis()->GetBinCenter(iy)<<"  "<<hist->GetZaxis()->GetBinCenter(iz)<<endl;
  //     }
  //   }
  // }
  

  outFile->Write();
  outFile->Close();

  /**
   * LaTex table line output
   */
  string outFilePathLatex = outFilePathBase + "/tex";
  system(("mkdir -p " + outFilePathLatex).c_str());
  outFilePathLatex += "/";
  outFilePathLatex += plotName;
  outFilePathLatex += ".tex";

  if (verbose) {
    cout << "INFO: Saving LaTex table line file here:" << endl;
    cout << "      " << outFilePathLatex << endl;
  }

  std::ofstream tableOutput;
  tableOutput.open(outFilePathLatex);
  tableOutput << std::fixed;
  tableOutput << std::setprecision(3);

  string parName = "";
  tableOutput << "& ";

  for (int i = 0; i < func->GetNpar(); ++i) {
    parName = func->GetParName(i);
    if (!showC0 && parName.compare("C_{0}") == 0) {
      continue;
    }
    if (!useEps && parName.compare("#epsilon") == 0) {
      continue;
    }
    tableOutput << func->GetParameter(i) << " $\\pm$ " << func->GetParError(i);
    tableOutput << " & ";
  }
  tableOutput << func->GetChisquare() / func->GetNDF();
  tableOutput << " \\\\";
  tableOutput << endl;
  tableOutput.close();


  /**
   * CSV output
   */
  string outFilePathCsv = outFilePathBase + "/csv";
  system(("mkdir -p " + outFilePathCsv).c_str());
  outFilePathCsv += "/";
  outFilePathCsv += plotName;
  outFilePathCsv += ".csv";

  if (verbose) {
    cout << "INFO: Saving CSV file here:" << endl;
    cout << "      " << outFilePathCsv << endl;
  }
//outputCSV(hist, outFilePathCsv);
  outputCSV(func, outFilePathCsv);
}


void Plot() {
  cout << "INFO: Plotting..." << endl;

  string outFilePathBase = "output";
  string today = Today();
  string c2_func = "c2_" + std::to_string(c2Index);
  TString outFilePath;
  vector<string> formatVec;
  formatVec.emplace_back("pdf");
  formatVec.emplace_back("png");

  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s", outFilePathBase.c_str(),
                                 today.c_str(),
                                 c2_func.c_str(),
                                 formatVec.at(i).c_str());
    system(("mkdir -p " + outFilePath).Data());
    if (verbose) {
      cout << "INFO: Saving " << formatVec.at(i) << " files here:" << endl;
      cout << "      " << outFilePath << endl;
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 50, 50, 600, 600);

  /**
   * 1D Plots
   */
  canvas->SetTopMargin(.05);
  canvas->SetRightMargin(.05);
  canvas->SetBottomMargin(.11);
  canvas->SetLeftMargin(.15);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  std::string histName = hist->GetName();

  /**
   * Q_out
   */
  TH1D* histX;
  if (histData.find("_cube_") != std::string::npos) {
    TH1D* histXpart1 = dynamic_cast<TH1D*>(hist->ProjectionX(
                       (histName + "_out_part1").c_str(),
                       1, 9, 1, 9));
    TH1D* histXpart2 = dynamic_cast<TH1D*>(hist->ProjectionX(
                       (histName + "_out_part2").c_str(),
                       1, 9, 1, 9));
    TH1D* histXpart3 = dynamic_cast<TH1D*>(hist->ProjectionX(
                       (histName + "_out_part3").c_str(),
                       1, 9, 1, 9));
    TH1D* histXpart4 = dynamic_cast<TH1D*>(hist->ProjectionX(
                       (histName + "_out_part4").c_str(),
                       1, 27, 1, 27));

    histX = new TH1D((histName + "_out").c_str(),
                     hist->GetTitle(),
                     varQbinsSize, &varQbins[0]);

    for (int i = 1; i <= 15; ++i) {
      histX->SetBinContent(i, histXpart1->GetBinContent(i) / (9 * 9));
      histX->SetBinError(i, histXpart1->GetBinError(i) / (9 * 9));
    }

    histX->SetBinContent(16, histXpart2->GetBinContent(17) / (3 * 3));
    histX->SetBinError(16, histXpart2->GetBinError(17) / (3 * 3));
    histX->SetBinContent(17, histXpart2->GetBinContent(20) / (3 * 3));
    histX->SetBinError(17, histXpart2->GetBinError(20) / (3 * 3));
    histX->SetBinContent(18, histXpart2->GetBinContent(23) / (3 * 3));
    histX->SetBinError(18, histXpart2->GetBinError(23) / (3 * 3));
    histX->SetBinContent(19, histXpart2->GetBinContent(26) / (3 * 3));
    histX->SetBinError(19, histXpart2->GetBinError(26) / (3 * 3));

    histX->SetBinContent(20, histXpart3->GetBinContent(32));
    histX->SetBinError(20, histXpart3->GetBinError(32));
    histX->SetBinContent(21, histXpart3->GetBinContent(41));
    histX->SetBinError(21, histXpart3->GetBinError(41));
    histX->SetBinContent(22, histXpart3->GetBinContent(50));
    histX->SetBinError(22, histXpart3->GetBinError(50));

    histX->SetBinContent(23, histXpart4->GetBinContent(68));
    histX->SetBinError(23, histXpart4->GetBinError(68));
    histX->SetBinContent(24, histXpart4->GetBinContent(95));
    histX->SetBinError(24, histXpart4->GetBinError(95));
    histX->SetBinContent(25, histXpart4->GetBinContent(122));
    histX->SetBinError(25, histXpart4->GetBinError(122));
    histX->SetBinContent(26, histXpart4->GetBinContent(149));
    histX->SetBinError(26, histXpart4->GetBinError(149));
    histX->SetBinContent(27, histXpart4->GetBinContent(176));
    histX->SetBinError(27, histXpart4->GetBinError(176));

    delete histXpart1;
    delete histXpart2;
    delete histXpart3;
    delete histXpart4;
  } else {
    histX = dynamic_cast<TH1D*>(hist->ProjectionX(
                  (histName + "_out").c_str(),
                  projMin, projMax,
                  projMin, projMax));
    histX->Scale(1. / (projRange * projRange));
    /*histDataStruc = dynamic_cast<TH1D*>(hist->ProjectionX(
                  (histName + "_out").c_str(),
                  projMin, projMax,
                  projMin, projMax));
    histDataStruc->Reset("ICESM");*/
  }
  /*
  unsigned int n=0;
  int maxn = 100000;//DataStruc.GetFitData().GetNpoints();
  double val = 1;
  double invErr = 1;
  double value = 1;

  // define edges and number of bins for qMin
  std::vector<double> qMin3D = {0.0, 0.28, 0.6, 1.2, 4.2};
  std::vector<int> qNBins3D = {14, 15, 10, 7, 0};

  int nBins = 0;
  int projRan = 0;
  vector<double> binVal;
  vector<double> binErr;
  vector<double> newEdges;
  binVal.push_back(0);
  binErr.push_back(0);
  newEdges.push_back(0);
  for (int i = 0; i < qMin3D.size()-1; i++) {
    int bins = qNBins3D[i]*((qMin3D[i+1]-qMin3D[i])/qMin3D[i+1]);
    nBins += bins;
    double binWidth = qMin3D[i+1] / qNBins3D[i];
    for (int j = 0; j < bins; j++) {
      newEdges.push_back(qMin3D[i] + j*binWidth);
      cout<<qMin3D[i] + j*binWidth;
      binVal.push_back(0);
      binErr.push_back(0);
    }
  }
  TH1D* histDataStruc = new TH1D("histDataStruc", "histDataStruc", nBins, newEdges.data());
  histDataStruc->Sumw2();
    while (n<=DataStruc.Size()) {
      bool skip = true;
      const double *x = DataStruc.GetPoint(n, val, invErr);
      n++;
      for (int i = 0; i < 3; i++) {
        if (x[i] < 0) {
            skip = false;
            break;
        }
      } 
      if (!skip) continue;
      if ((((x[1]>projQmax/1000)||(x[2]>projQmax/1000))&&(x[0]<0.6))||(val==0)) {
        cout<<"skipping point N= "<<n<<"value= "<<val;
        for (int i = 0; i < 3; i++) {
          cout << ", x[" << i << "] = " << x[i];
        }
        cout<<endl;
        continue;
      }
      if (((x[1]>0.6)||(x[2]>0.6))&&(x[0]>=0.6))
        continue;
      //if ((x[1]>projQmax)&&(x[2]>projQmax)) break;
      cout<<" N= "<<n<<"value= "<<val;
      for (int i = 0; i < 3; i++) {
        cout << ", x[" << i << "] = " << x[i];
      }
      cout<<endl;
      projRan = 0;
      for (int i = 0; i<qMin3D.size(); i++) {
        if ((qMin3D[i]<=x[0])&&(qMin3D[i+1]>x[0]))
          projRan = i;
      }
      switch(projRan){
        case 0: case 1: projRan = projRange/(projRan+1); //projRange/(100*qMin3D[i+1]/qNBins3D[i])*2;
          break;
        case 2: projRan = 5;
          break;
        case 3: projRan = 1;
          break;
      }
      binVal.at(histDataStruc->FindBin(x[0])-1) += val;
      binErr.at(histDataStruc->FindBin(x[0])-1) += 1.0/(invErr*invErr);
      cout<<"Binerror"<<binErr.at(histDataStruc->FindBin(x[0])-1)<<endl;
    //histDataStruc->Fill(x[0],val/(projRan * projRan));
    }
    
    //histDataStruc->Scale(1. / (projRange * projRange));
    cout<<"N = "<<n<<endl;
    for (int i = 1; i<histDataStruc->GetNbinsX();i++) {
      for (int j = 0; j<qMin3D.size(); j++) {
        if ((qMin3D[j]<=histDataStruc->GetBinCenter(i))&&(qMin3D[j+1]>histDataStruc->GetBinCenter(i)))
          projRan = j;
      }
      switch(projRan){
        case 0: case 1: projRan = projRange/(projRan+1); //projRange/(100*qMin3D[i+1]/qNBins3D[i])*2;
          break;
        case 2: projRan = 5;
          break;
        case 3: projRan = 1;
          break;
      }
      histDataStruc->SetBinContent(i, binVal.at(i-1)/(projRan * projRan));
      histDataStruc->SetBinError(i, sqrt(binErr.at(i-1))/(projRan * projRan));
      cout<<"Bin "<<i<<", xaxis "<<histDataStruc->GetXaxis()->GetBinCenter(i)<<", value: "<<histDataStruc->GetBinContent(i)<<endl;
    }*/

  histX->Print();

  // TH1D* histDataStruc = makeProjectionHist(DataStruc_0, "out", projQmax, projRange, qMin3D_0.size()-1);
  TH1D* histDataStruc = makeProjectionHist(DataStruc, "out", projQmax, projRange, qMin3D_02.size()-1);
  histDataStruc->SetLineColor(6);
  histDataStruc->SetMarkerColor(4);
  histX->SetLineColor(4);
  histX->SetMarkerColor(4);

  histX->SetMinimum(yMin);
  histX->SetMaximum(yMax);

  histX->SetTitle("");
  histX->GetXaxis()->SetTitle("Q_{out} [GeV]");
  if (r2) {
    histX->GetYaxis()->SetTitle("R_{2} (Q_{out})");
  } else {
    histX->GetYaxis()->SetTitle("C_{2} (Q_{out})");
  }

  histX->GetXaxis()->SetTitleSize(0.05);
  histX->GetYaxis()->SetTitleSize(0.05);

  TGraph* graphX;
  if (histData.find("_cube_") != std::string::npos) {
    graphX = makeProjection(hist, func, "out", qMin, qMax, 1, 9);
  } else {
    graphX = makeProjection(hist, func, "out", qMin, qMax, projMin, projMax);
  }
  graphX->SetLineColor(2);
  graphX->SetLineWidth(2);

  TPaveText *textX = new TPaveText(.53, .53, .93, .93, "brNDC");
  textX->SetFillStyle(0);
  textX->SetFillColor(0);
  textX->SetBorderSize(0);
  textX->SetTextColor(1);
  textX->SetTextFont(42);
  textX->SetTextAlign(11);
  textX->SetTextSize(0.045);
  addFitResults(textX, func);

  TPaveText *text2X = new TPaveText(.53, .14, .93, .38, "brNDC");
  text2X->SetFillStyle(0);
  text2X->SetFillColor(0);
  text2X->SetBorderSize(0);
  text2X->SetTextColor(1);
  text2X->SetTextFont(42);
  text2X->SetTextAlign(11);
  text2X->SetTextSize(0.045);
  for (size_t i = 0; i < commentVec.size(); ++i) {
    text2X->AddText(commentVec.at(i).c_str());
  }
  text2X->AddText(c2IndexText);
  if (histData.find("_cube_") != std::string::npos) {
    text2X->AddText("0 #leq Q_{side} < 180 MeV");
    text2X->AddText("0 #leq Q_{long} < 180 MeV");
  } else {
    text2X->AddText(qSideText.Data());
    text2X->AddText(qLongText.Data());
  }

  TLine* lineX = new TLine(histX->GetXaxis()->GetXmin(), 1.,
                           histX->GetXaxis()->GetXmax(), 1.);
  lineX->SetLineStyle(7);

  histX->Draw();
  histDataStruc->Draw("same");
  graphX->Draw("Csame");
  lineX->Draw();
  textX->Draw();
  text2X->Draw();

  if (rejTo > rejFrom) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar1 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
    ar1 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  }
  if (rej_OSL) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar2 = new TArrow(rejFromOsl[0], ar_y, rejToOsl[0], ar_y, 0.015, "<>");
    ar2 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFromOsl[0] + rejToOsl[0]) - 0.1, ar_y, "excluded");
  }

  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_out.%s", outFilePathBase.c_str(),
                                           today.c_str(),
                                           c2_func.c_str(),
                                           formatVec.at(i).c_str(),
                                           plotName.c_str(),
                                           formatVec.at(i).c_str());
    canvas->Print(outFilePath);
  }

  delete histX;
  delete graphX;
  delete textX;
  delete text2X;
  delete lineX;

  /**
   * Q_side
   */
  TH1D* histY;
  if (histData.find("_cube_") != std::string::npos) {
    TH1D* histYpart1 = dynamic_cast<TH1D*>(hist->ProjectionY(
                       (histName + "_side_part1").c_str(),
                       1, 9, 1, 9));
    TH1D* histYpart2 = dynamic_cast<TH1D*>(hist->ProjectionY(
                       (histName + "_side_part2").c_str(),
                       1, 9, 1, 9));
    TH1D* histYpart3 = dynamic_cast<TH1D*>(hist->ProjectionY(
                       (histName + "_side_part3").c_str(),
                       1, 9, 1, 9));
    TH1D* histYpart4 = dynamic_cast<TH1D*>(hist->ProjectionY(
                       (histName + "_side_part4").c_str(),
                       1, 27, 1, 27));

    histY = new TH1D((histName + "_side").c_str(),
                     hist->GetTitle(),
                     varQbinsSize, &varQbins[0]);

    for (int i = 1; i <= 15; ++i) {
      histY->SetBinContent(i, histYpart1->GetBinContent(i) / (9 * 9));
      histY->SetBinError(i, histYpart1->GetBinError(i) / (9 * 9));
    }

    histY->SetBinContent(16, histYpart2->GetBinContent(17) / (3 * 3));
    histY->SetBinError(16, histYpart2->GetBinError(17) / (3 * 3));
    histY->SetBinContent(17, histYpart2->GetBinContent(20) / (3 * 3));
    histY->SetBinError(17, histYpart2->GetBinError(20) / (3 * 3));
    histY->SetBinContent(18, histYpart2->GetBinContent(23) / (3 * 3));
    histY->SetBinError(18, histYpart2->GetBinError(23) / (3 * 3));
    histY->SetBinContent(19, histYpart2->GetBinContent(26) / (3 * 3));
    histY->SetBinError(19, histYpart2->GetBinError(26) / (3 * 3));

    histY->SetBinContent(20, histYpart3->GetBinContent(32));
    histY->SetBinError(20, histYpart3->GetBinError(32));
    histY->SetBinContent(21, histYpart3->GetBinContent(41));
    histY->SetBinError(21, histYpart3->GetBinError(41));
    histY->SetBinContent(22, histYpart3->GetBinContent(50));
    histY->SetBinError(22, histYpart3->GetBinError(50));

    histY->SetBinContent(23, histYpart4->GetBinContent(68));
    histY->SetBinError(23, histYpart4->GetBinError(68));
    histY->SetBinContent(24, histYpart4->GetBinContent(95));
    histY->SetBinError(24, histYpart4->GetBinError(95));
    histY->SetBinContent(25, histYpart4->GetBinContent(122));
    histY->SetBinError(25, histYpart4->GetBinError(122));
    histY->SetBinContent(26, histYpart4->GetBinContent(149));
    histY->SetBinError(26, histYpart4->GetBinError(149));
    histY->SetBinContent(27, histYpart4->GetBinContent(176));
    histY->SetBinError(27, histYpart4->GetBinError(176));

    delete histYpart1;
    delete histYpart2;
    delete histYpart3;
    delete histYpart4;
  } else {
    histY = dynamic_cast<TH1D*>(hist->ProjectionY(
                  (histName + "_side").c_str(),
                  projMin, projMax,
                  projMin, projMax));
    histY->Scale(1. / (projRange * projRange));
  }

  // histDataStruc = makeProjectionHist(DataStruc_0, "side", projQmax, projRange, qMin3D_0.size()-1);
  histDataStruc = makeProjectionHist(DataStruc, "side", projQmax, projRange, qMin3D_02.size()-1);
  histDataStruc->SetLineColor(6);
  histDataStruc->SetMarkerColor(4);
  histY->SetLineColor(4);
  histY->SetMarkerColor(4);

  histY->SetMinimum(yMin);
  histY->SetMaximum(yMax);

  histY->SetTitle("");
  histY->GetXaxis()->SetTitle("Q_{side} [GeV]");
  if (r2) {
    histY->GetYaxis()->SetTitle("R_{2} (Q_{side})");
  } else {
    histY->GetYaxis()->SetTitle("C_{2} (Q_{side})");
  }
  histY->GetXaxis()->SetTitleOffset(1.0);

  histY->GetXaxis()->SetTitleSize(0.05);
  histY->GetYaxis()->SetTitleSize(0.05);

  TGraph* graphY;
  if (histData.find("_cube_") != std::string::npos) {
    graphY = makeProjection(hist, func, "side", qMin, qMax, 1, 9);
  } else {
    graphY = makeProjection(hist, func, "side", qMin, qMax, projMin, projMax);
  }
  graphY->SetLineColor(2);
  graphY->SetLineWidth(2);

  TPaveText *textY = new TPaveText(.53, .53, .93, .93, "brNDC");
  textY->SetFillStyle(0);
  textY->SetFillColor(0);
  textY->SetBorderSize(0);
  textY->SetTextColor(1);
  textY->SetTextFont(42);
  textY->SetTextAlign(11);
  textY->SetTextSize(0.045);
  addFitResults(textY, func);

  TPaveText *text2Y = new TPaveText(.53, .14, .93, .38, "brNDC");
  text2Y->SetFillStyle(0);
  text2Y->SetFillColor(0);
  text2Y->SetBorderSize(0);
  text2Y->SetTextColor(1);
  text2Y->SetTextFont(42);
  text2Y->SetTextAlign(11);
  text2Y->SetTextSize(0.045);
  for (size_t i = 0; i < commentVec.size(); ++i) {
    text2Y->AddText(commentVec.at(i).c_str());
  }
  text2Y->AddText(c2IndexText);
  if (histData.find("_cube_") != std::string::npos) {
    text2Y->AddText("0 #leq Q_{out} < 180 MeV");
    text2Y->AddText("0 #leq Q_{long} < 180 MeV");
  } else {
    text2Y->AddText(qOutText.Data());
    text2Y->AddText(qLongText.Data());
  }

  TLine* lineY = new TLine(histY->GetXaxis()->GetXmin(), 1.,
                           histY->GetXaxis()->GetXmax(), 1.);
  lineY->SetLineStyle(7);

  histY->Draw();
  histDataStruc->Draw("same");
  graphY->Draw("Csame");
  lineY->Draw();
  textY->Draw();
  text2Y->Draw();

  if (rejTo > rejFrom) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar2 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
    ar2 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  }
  if (rej_OSL) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar2 = new TArrow(rejFromOsl[1], ar_y, rejToOsl[1], ar_y, 0.015, "<>");
    ar2 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFromOsl[1] + rejToOsl[1]) - 0.1, ar_y, "excluded");
  }

  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_side.%s", outFilePathBase.c_str(),
                                            today.c_str(),
                                            c2_func.c_str(),
                                            formatVec.at(i).c_str(),
                                            plotName.c_str(),
                                            formatVec.at(i).c_str());
    canvas->Print(outFilePath);
  }

  delete histY;
  delete graphY;
  delete textY;
  delete text2Y;
  delete lineY;

  /**
   * Q_long
   */
  TH1D* histZ;
  if (histData.find("_cube_") != std::string::npos) {
    TH1D* histZpart1 = dynamic_cast<TH1D*>(hist->ProjectionZ(
                       (histName + "_long_part1").c_str(),
                       1, 9, 1, 9));
    TH1D* histZpart2 = dynamic_cast<TH1D*>(hist->ProjectionZ(
                       (histName + "_long_part2").c_str(),
                       1, 9, 1, 9));
    TH1D* histZpart3 = dynamic_cast<TH1D*>(hist->ProjectionZ(
                       (histName + "_long_part3").c_str(),
                       1, 9, 1, 9));
    TH1D* histZpart4 = dynamic_cast<TH1D*>(hist->ProjectionZ(
                       (histName + "_long_part4").c_str(),
                       1, 27, 1, 27));

    histZ = new TH1D((histName + "_long").c_str(),
                     hist->GetTitle(),
                     varQbinsSize, &varQbins[0]);

    for (int i = 1; i <= 15; ++i) {
      histZ->SetBinContent(i, histZpart1->GetBinContent(i) / (9 * 9));
      histZ->SetBinError(i, histZpart1->GetBinError(i) / (9 * 9));
    }

    histZ->SetBinContent(16, histZpart2->GetBinContent(17) / (3 * 3));
    histZ->SetBinError(16, histZpart2->GetBinError(17) / (3 * 3));
    histZ->SetBinContent(17, histZpart2->GetBinContent(20) / (3 * 3));
    histZ->SetBinError(17, histZpart2->GetBinError(20) / (3 * 3));
    histZ->SetBinContent(18, histZpart2->GetBinContent(23) / (3 * 3));
    histZ->SetBinError(18, histZpart2->GetBinError(23) / (3 * 3));
    histZ->SetBinContent(19, histZpart2->GetBinContent(26) / (3 * 3));
    histZ->SetBinError(19, histZpart2->GetBinError(26) / (3 * 3));

    histZ->SetBinContent(20, histZpart3->GetBinContent(32));
    histZ->SetBinError(20, histZpart3->GetBinError(32));
    histZ->SetBinContent(21, histZpart3->GetBinContent(41));
    histZ->SetBinError(21, histZpart3->GetBinError(41));
    histZ->SetBinContent(22, histZpart3->GetBinContent(50));
    histZ->SetBinError(22, histZpart3->GetBinError(50));

    histZ->SetBinContent(23, histZpart4->GetBinContent(68));
    histZ->SetBinError(23, histZpart4->GetBinError(68));
    histZ->SetBinContent(24, histZpart4->GetBinContent(95));
    histZ->SetBinError(24, histZpart4->GetBinError(95));
    histZ->SetBinContent(25, histZpart4->GetBinContent(122));
    histZ->SetBinError(25, histZpart4->GetBinError(122));
    histZ->SetBinContent(26, histZpart4->GetBinContent(149));
    histZ->SetBinError(26, histZpart4->GetBinError(149));
    histZ->SetBinContent(27, histZpart4->GetBinContent(176));
    histZ->SetBinError(27, histZpart4->GetBinError(176));

    delete histZpart1;
    delete histZpart2;
    delete histZpart3;
    delete histZpart4;
  } else {
    histZ = dynamic_cast<TH1D*>(hist->ProjectionZ(
                  (histName + "_long").c_str(),
                  projMin, projMax,
                  projMin, projMax));
    histZ->Scale(1. / (projRange * projRange));
  }

  // histDataStruc = makeProjectionHist(DataStruc_0, "long", projQmax, projRange, qMin3D_0.size()-1);
  histDataStruc = makeProjectionHist(DataStruc, "long", projQmax, projRange, qMin3D_02.size()-1);
  histDataStruc->SetLineColor(6);
  histDataStruc->SetMarkerColor(4);
  histZ->SetLineColor(4);
  histZ->SetMarkerColor(4);

  histZ->SetMinimum(yMin);
  histZ->SetMaximum(yMax);

  histZ->SetTitle("");
  histZ->GetXaxis()->SetTitle("Q_{long} [GeV]");
  if (r2) {
    histZ->GetYaxis()->SetTitle("R_{2} (Q_{long})");
  } else {
    histZ->GetYaxis()->SetTitle("C_{2} (Q_{long})");
  }

  histZ->GetXaxis()->SetTitleSize(0.05);
  histZ->GetYaxis()->SetTitleSize(0.05);

  TGraph* graphZ;
  if (histData.find("_cube_") != std::string::npos) {
    graphZ = makeProjection(hist, func, "long", qMin, qMax, projMin, projMax);
  } else {
    graphZ = makeProjection(hist, func, "long", qMin, qMax, projMin, projMax);
  }
  graphZ->SetLineColor(2);
  graphZ->SetLineWidth(2);

  TPaveText *textZ = new TPaveText(.53, .53, .93, .93, "brNDC");
  textZ->SetFillStyle(0);
  textZ->SetFillColor(0);
  textZ->SetBorderSize(0);
  textZ->SetTextColor(1);
  textZ->SetTextFont(42);
  textZ->SetTextAlign(11);
  textZ->SetTextSize(0.045);
  addFitResults(textZ, func);

  TPaveText *text2Z = new TPaveText(.53, .14, .93, .38, "brNDC");
  text2Z->SetFillStyle(0);
  text2Z->SetFillColor(0);
  text2Z->SetBorderSize(0);
  text2Z->SetTextColor(1);
  text2Z->SetTextFont(42);
  text2Z->SetTextAlign(11);
  text2Z->SetTextSize(0.045);
  for (size_t i = 0; i < commentVec.size(); ++i) {
    text2Z->AddText(commentVec.at(i).c_str());
  }
  text2Z->AddText(c2IndexText);
  if (histData.find("_cube_") != std::string::npos) {
    text2Z->AddText("0 #leq Q_{out} < 180 MeV");
    text2Z->AddText("0 #leq Q_{side} < 180 MeV");
  } else {
    text2Z->AddText(qOutText.Data());
    text2Z->AddText(qSideText.Data());
  }

  TLine* lineZ = new TLine(histZ->GetXaxis()->GetXmin(), 1.,
                           histZ->GetXaxis()->GetXmax(), 1.);
  lineZ->SetLineStyle(7);

  histZ->Draw();
  histDataStruc->Draw("same");
  graphZ->Draw("Csame");
  lineZ->Draw();
  textZ->Draw();
  text2Z->Draw();

  if (rejTo > rejFrom) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar3 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
    ar3 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  }
  if (rej_OSL) {
    ar_y = yMin + 0.14 * (yMax - yMin);
    TArrow *ar3 = new TArrow(rejFromOsl[2], ar_y, rejToOsl[2], ar_y, 0.015, "<>");
    ar3 -> Draw();

    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    lexc.DrawLatex(0.5 * (rejFromOsl[2] + rejToOsl[2]) - 0.1, ar_y, "excluded");
  }

  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_long.%s", outFilePathBase.c_str(),
                                            today.c_str(),
                                            c2_func.c_str(),
                                            formatVec.at(i).c_str(),
                                            plotName.c_str(),
                                            formatVec.at(i).c_str());
    canvas->Print(outFilePath);
  }

  delete histZ;
  delete graphZ;
  delete textZ;
  delete text2Z;
  delete lineZ;
  delete canvas;

  /**
   * 2D Plots
   */

  // Q_out_side
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          "out_side",
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "out_side", commentVec, outFilePath,
                     projMin, projMax);
  }

  // Q_side_long
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          "side_long",
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "side_long", commentVec, outFilePath,
                     projMin, projMax);
  }

  // Q_out_long
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          "out_long",
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "out_long", commentVec, outFilePath,
                     projMin, projMax);
  }
}
