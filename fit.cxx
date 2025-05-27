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
#include "THStack.h"
#include <TLegend.h>
// Fit
#include "fit.h"
#include "c2.h"
#include "utils.h"
#include "TError.h"
#include "Fit/BinData.h"
#include "Fit/DataRange.h"
#include <Fit/DataOptions.h>
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"
#include <Math/WrappedParamFunction.h>
#include <TMatrixDSym.h>
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <Fit/Chi2FCN.h>
#include <Math/IntegratorMultiDim.h>

// std
using std::cout;
using std::endl;
using std::string;
using std::vector;

bool help = false;
bool verbose = false;//true;//
double cubes = true;

string fileData = "";
    // "input/bec_alln_W_trPt_100_0_Qcut_0_Qbin_20_100_0_3D_2_GeV.root";
string file2Data = "";
string fileMC ="";
    // "input/bec_mc_alln_W_sc_trPt_100_0_Qcut_0_Qbin_20_100_0_3D_2_GeV.root";
string file2MC = "";

string histData = "ppmm_qosl_g";
string hist2Data = "";
string histMC = "";
string hist2MC = "";
bool r2 = false;

bool reject_final = false;
bool mean_fit = false;

int c2Index = 1;
string c2_func ="";
double qMin = -1.;
double qMax = -1.;
double rejFrom = 1.;
double rejTo = -1.;
// vector<double> rejFromOsl = {1., 1., 1.};
// vector<double> rejToOsl = {-1., -1., -1.};
std::vector<std::string> sets = {"", "2", "3", "4", "5"};
std::vector<std::string> types = {"out", "side", "long"};
vector<vector<double>> rejFromOsl;//= {{1., 1., 1.}, {1., 1., 1.}, {1., 1., 1.}};
vector<double> rejOsl_temp;
vector<vector<double>> rejToOsl;//= {{-1., -1., -1.}, {-1., -1., -1.}, {1., 1., 1.}};
vector<double> alpha = {-1., -1., -1.};
// rej[0,1,2] is equivalent to rej[out,side,long]
double rej2From = 1.;
double rej2To = -1.;
bool useEps = false;
bool showC0 = false;
bool fixC0 = false;
bool enlargeUncert = false;
int nParams = 0;
bool with_errors = false;
string fitParams = "RE";
std::pair<double, std::string> chi_firstBins = {0.0, "PEAK"};

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
FitResultData fitData;

ROOT::Fit::BinData DataStruc;
ROOT::Fit::BinData DataStruc_0;
int npoints_DataStruc;
vector<std::unique_ptr<TH3D>> InHistRatio;
vector<double> qBins3D_width = {0.02, 0.04, 0.12, 0.6};
vector<double> qMin3D_0 = {.0, .28, .6, 1.2, 4.2};
vector<double> qMin3D_02 = {0.02, .28, .6, 1.2, 4.2};
// std::vector<TH3D*> InHistRatio;

// std::vector<double> varQbins = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
//                                 0.22, 0.24, 0.26, 0.28, 0.30, 0.35, 0.4, 0.45, 0.5, 0.6,
//                                 0.7, 0.8, 1.0, 1.2, 1.5, 2.0};
vector<double> varQbins = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 
                                0.32, 0.36, 0.4, 0.45, 0.48, 0.52, 0.56, 0.6,
                                0.72, 0.84, 0.96, 1.08, 1.2, 1.8};
/*{0., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200.,
                     220., 240., 260., 280., 300., 360., 420., 480., 540., 720.,
900., 1080., 1620., 2160., 2700., 3240., 3780.};*/
int varQbinsSize = varQbins.size() - 1;

string outFilePathBase = "output/";

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
  if (input.cmdOptionExists("--rej-fin")) {
    reject_final = true;
  }
  if (input.cmdOptionExists("--rej-from")) {
    string param = input.getCmdOption("--rej-from");
    rejFrom = std::stod(param);
  }
  
  for (const auto& set : sets) {
    for (const auto& type : types) {
      if (input.cmdOptionExists("--rej" + set + "-from-" + type)) {
        string param = input.getCmdOption("--rej" + set + "-from-" + type);
        rejOsl_temp.push_back(param.empty() ? 1.0 : std::stod(param)); // Default value 1.0 if not defined
      }
    }
    if (rejOsl_temp != std::vector<double>{1.0, 1.0, 1.0}) {
      rejFromOsl.push_back(rejOsl_temp); 
    }
    rejOsl_temp.clear();

    for (const auto& type : types) {
      if (input.cmdOptionExists("--rej" + set + "-to-" + type)) {
        string param = input.getCmdOption("--rej" + set + "-to-" + type);
        rejOsl_temp.push_back(param.empty() ? -1.0 : std::stod(param)); // Default value -1.0 if not defined
      }
    }
    if (rejOsl_temp != std::vector<double>{-1.0, -1.0, -1.0}) {
      rejToOsl.push_back(rejOsl_temp);
    }
    rejOsl_temp.clear();
  }

  if (input.cmdOptionExists("--rej-to")) {
    string param = input.getCmdOption("--rej-to");
    rejTo = std::stod(param);
  }
  if (input.cmdOptionExists("--rej2-from")) {
    string param = input.getCmdOption("--rej2-from");
    rej2From = std::stod(param);
  }
  if (input.cmdOptionExists("--rej2-to")) {
    string param = input.getCmdOption("--rej2-to");
    rej2To = std::stod(param);
  }
  if (input.cmdOptionExists("--alpha-out")) {
    string param = input.getCmdOption("--alpha-out");
    alpha[0] = std::stod(param);
  }
  if (input.cmdOptionExists("--alpha-side")) {
    string param = input.getCmdOption("--alpha-side");
    alpha[1] = std::stod(param);
  }
  if (input.cmdOptionExists("--alpha-long")) {
    string param = input.getCmdOption("--alpha-long");
    alpha[2] = std::stod(param);
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

  //define range of bindata object, min and max defined with qMin and qMax
  ROOT::Fit::DataOptions DataStruc_options;
  ROOT::Fit::DataRange DataStruc_range(3);

  DataStruc_options.fUseRange = true;
  double DataStruc_qmin = qMin > 0 ? qMin : 0;
  double DataStruc_qmax = qMax > 0 ? qMax : 4.2;
  // DataStruc_range.SetRange(DataStruc_qmin, DataStruc_qmax, DataStruc_qmin, DataStruc_qmax, DataStruc_qmin, DataStruc_qmax);
  DataStruc_range.SetRange(0, DataStruc_qmin, DataStruc_qmax); // X range
  DataStruc_range.SetRange(1, DataStruc_qmin, DataStruc_qmax); // Y range
  DataStruc_range.SetRange(2, DataStruc_qmin, DataStruc_qmax); // Z range
  cout << "DataStruc_qmin: " << DataStruc_qmin << " DataStruc_qmax: " << DataStruc_qmax << endl;


  DataStruc = ROOT::Fit::BinData (DataStruc_options, DataStruc_range, 100000, 3);
  // DataStruc_0 = ROOT::Fit::BinData (DataStruc_options, DataStruc_range, 100000, 3);
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

void configure3DAxis(TH3* histos, bool stats = false, double offset1 = 1.8, double offset2 = 2.1, double offset3 = 1.2) {
  if (!histos) {
      std::cerr << "Error: Histogram is a null pointer!" << std::endl;
      return;
  }
  histos->SetStats(stats);
  // Center the axis titles
  histos->GetXaxis()->CenterTitle();
  histos->GetYaxis()->CenterTitle();
  histos->GetZaxis()->CenterTitle();

  // Set the title offsets
  histos->GetXaxis()->SetTitleOffset(offset1);
  histos->GetYaxis()->SetTitleOffset(offset2);
  histos->GetZaxis()->SetTitleOffset(offset3);
}
// bool rej_OSL() {
//   std::vector<double> diff(rejToOsl.size());
//   std::transform(rejToOsl.begin(), rejToOsl.end(), rejFromOsl.begin(), diff.begin(), std::minus<double>());
//   return std::all_of(diff.begin(), diff.end(), [](double i){ return i > 0; });
// }
std::vector<bool> rej_OSL() {
  // Lambda function to check if all elements in the difference vector are greater than 0
  if (rejToOsl.empty() || rejFromOsl.empty()) {
    return std::vector<bool>();
  }
  auto check_vector = [](const std::vector<double>& toOslVec, const std::vector<double>& fromOslVec) {
    std::vector<double> diff(toOslVec.size());
    std::transform(toOslVec.begin(), toOslVec.end(), fromOslVec.begin(), diff.begin(), std::minus<double>());
    return std::all_of(diff.begin(), diff.end(), [](double i){ return i > 0; });
  };

  // Vector to store the results
  std::vector<bool> results(rejToOsl.size());

  // Ensure rejToOsl and rejFromOsl have the same size
  if (rejToOsl.size() != rejFromOsl.size()) {
    throw std::runtime_error("rejToOsl and rejFromOsl must have the same size");
  }

  // Apply the lambda function to each pair of vectors from rejToOsl and rejFromOsl
  std::transform(rejToOsl.begin(), rejToOsl.end(), rejFromOsl.begin(), results.begin(), check_vector);

  // Print the results
  // for (bool res : results) {
  //   std::cout << "Result: " << res << std::endl;
  // }
  // Resize the results vector to only include the relevant rejection regions
  results.erase(std::remove(results.begin(), results.end(), false), results.end());

  // Return the results vector
  return results;
}

double ComputeChi2(const std::vector<std::tuple<std::vector<double>, double, double>> fDataPoints,
                   TF3& func, 
                   const std::vector<double>& fQMin3D,
                   const std::vector<int>& fQNBins3D) {
    ROOT::Math::WrappedMultiTF1 wrappedFunc(func, 3);
    ROOT::Math::IntegratorMultiDim integrator(ROOT::Math::IntegrationMultiDim::kADAPTIVE);

    double chi2 = 0.0;
    for (const auto& point : fDataPoints) {
        const auto& coords = std::get<0>(point);
        double value = std::get<1>(point);
        double error = std::get<2>(point);

        // std::cout << "Coords: ";
        // for (const auto& coord : coords) {
        //     std::cout << coord << " ";
        // }
        // std::cout << "Value: " << value << " Error: " << error << std::endl;
        double maxCoord = *std::max_element(coords.begin(), coords.end());
        double binWidth = 0.0;
        for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
            if (fQMin3D[k] <= maxCoord && fQMin3D[k + 1] > maxCoord) {
                binWidth = (fQMin3D[k + 1] - fQMin3D[k]) / fQNBins3D[k];
                break;
            }
        }

        // std::cout << "Max coordinate: " << maxCoord << ", Bin width: " << binWidth << std::endl;
        double binMin[3], binMax[3];
        for (int j = 0; j < 3; ++j) {
            binMin[j] = coords[j] - binWidth / 2.0;
            binMax[j] = coords[j] + binWidth / 2.0;
        }

        double integral = integrator.Integral(wrappedFunc, binMin, binMax);
        double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
        double meanValue = integral / binVolume;

        double diff = value - meanValue;
        chi2 += (diff * diff) / (error * error);
    }

    return chi2;
}

double ComputeChi2ForFirstBins(const std::vector<std::tuple<std::vector<double>, double, double>>& fDataPoints,
                               TF3& func, 
                               const std::vector<double>& fQMin3D,
                               const std::vector<int>& fQNBins3D,
                               const int nBinsToConsider) {
    ROOT::Math::WrappedMultiTF1 wrappedFunc(func, 3);
    ROOT::Math::IntegratorMultiDim integrator(ROOT::Math::IntegrationMultiDim::kADAPTIVE);

    double chi2 = 0.0;
    int chis = 0;
    int Nfd = pow(nBinsToConsider,3) - func.GetNumberFreeParameters() + 2; //2 comes from C0 and eps not contributing to the first bins
    // cout<<"Nfd: "<<Nfd<<endl;
    for (const auto& point : fDataPoints) {
        const auto& coords = std::get<0>(point);
        double value = std::get<1>(point);
        double error = std::get<2>(point);

        // Restrict to the first nBinsToConsider bins in all 3 dimensions
        if (coords[0] > fQMin3D[0] + nBinsToConsider * (fQMin3D[1] - fQMin3D[0]) / fQNBins3D[0] ||
            coords[1] > fQMin3D[0] + nBinsToConsider * (fQMin3D[1] - fQMin3D[0]) / fQNBins3D[0] ||
            coords[2] > fQMin3D[0] + nBinsToConsider * (fQMin3D[1] - fQMin3D[0]) / fQNBins3D[0]) {
            continue; // Skip points outside the first 5 bins
        }

        double maxCoord = *std::max_element(coords.begin(), coords.end());
        double binWidth = 0.0;
        for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
            if (fQMin3D[k] <= maxCoord && fQMin3D[k + 1] > maxCoord) {
                binWidth = (fQMin3D[k + 1] - fQMin3D[k]) / fQNBins3D[k];
                break;
            }
        }

        double binMin[3], binMax[3];
        for (int j = 0; j < 3; ++j) {
            binMin[j] = coords[j] - binWidth / 2.0;
            binMax[j] = coords[j] + binWidth / 2.0;
        }

        double integral = integrator.Integral(wrappedFunc, binMin, binMax);
        double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
        double meanValue = integral / binVolume;

        double diff = value - meanValue;

        chi2 += (diff * diff) / (error * error)/ Nfd;
        // print current chi2 and total chi2
        // cout<<"coordinates: "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<" current chi2: "<<
        // (diff * diff) / (error * error)<<" total chi2: "<< chi2<< " chis: "<<chis<<endl;
        // chis++;
    }

    return chi2;
}

void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, ROOT::Fit::BinData& data, const std::vector<double>& qMin3D = {0., .28, .6, 1.2, 4.2});
// void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, vector<std::unique_ptr<TH3D>>& InHistRatiotest, ROOT::Fit::BinData& data) {
void HisttoData(vector<std::unique_ptr<TH3D>>& HistData, vector<std::unique_ptr<TH3D>>& HistMC, ROOT::Fit::BinData& data, const std::vector<double>& qMin3D) {
  double DataStruc_qmin = qMin > 0 ? qMin : 0;
  double DataStruc_qmax = qMax > 0 ? qMax : 4.2;
  cout<<"DataStruc_qmin: "<<DataStruc_qmin<<" DataStruc_qmax: "<<DataStruc_qmax<<endl;
  const int hist_number = qMin3D.size()-1;
  InHistRatio.clear();
  double ix_down_limit = 0, iy_down_limit = 0, iz_down_limit = 0;
  // bool reject_region = rej_OSL();
  std::vector<bool> reject_region = rej_OSL();
  // for (int i=0;i<4-hist_number;i++) {
  //   // push empty histograms to InHistRatio
  //   InHistRatio.push_back(std::make_unique<TH3D>());
  // }
  for(int i=0;i<hist_number;i++) {
    // auto inHistRatio = std::make_unique<TH3D>(*HistData[i]);
    std::unique_ptr<TH3D> inHistRatio;
    //both vectors cant be empty, the if in which this function is called checks for that
    if (!HistData.empty()) {
      cout<<"i before making histogram: "<<i<<endl;
      inHistRatio = std::make_unique<TH3D>(*HistData[i]);
      cout<<"Size of HistData: "<<inHistRatio->GetSize()<<endl;
      if (!HistMC.empty()) {
        inHistRatio->Divide(HistData[i].get(), HistMC[i].get());//,
                              // 1./HistData[i]->Integral(), 1./HistMC[i]->Integral());
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
    double inHistRatioMin = inHistRatio->GetMinimum();
    double inHistRatioMax = inHistRatio->GetMaximum();
    // std::vector<double> qMin3D = {0., .28, .6, 1.2, 4.2};
    if (qMin3D[0] > 0.06)
      throw "ERROR: lower limit must be less than 0.06!"; //otherwise we do not fit the BEC peak
    // std::vector<double> qMin3D = {lower_limit, .28, .6, 1.2, 4.2};
    if (i==0) {
      ix_down_limit = inHistRatio->GetXaxis()->FindBin(qMin3D[0]);
      cout<<"ix_down_limit "<<ix_down_limit<<endl;
      iy_down_limit = inHistRatio->GetYaxis()->FindBin(qMin3D[0]);
      iz_down_limit = inHistRatio->GetZaxis()->FindBin(qMin3D[0]);
    }
    // int ixStart = inHistRatio->GetXaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    // int iyStart = inHistRatio->GetYaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    // int izStart = inHistRatio->GetZaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int ixStart = inHistRatio->GetXaxis()->FindBin(qMin3D[i]);
    int iyStart = inHistRatio->GetYaxis()->FindBin(qMin3D[i]);
    int izStart = inHistRatio->GetZaxis()->FindBin(qMin3D[i]);
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
          bool reject = false;
          for (int k = 0; k < reject_region.size(); k++) {
            if (reject_region[k]) {
              if ((binCoord[0] >= rejFromOsl[k][0]) && (binCoord[0] <= rejToOsl[k][0]) &&
                  (binCoord[1] >= rejFromOsl[k][1]) && (binCoord[1] <= rejToOsl[k][1]) &&
                  (binCoord[2] >= rejFromOsl[k][2]) && (binCoord[2] <= rejToOsl[k][2])) {
                reject = true;
                break;
              }
            }
          }
          //check if the coords are in the range of qMin and qMax
          if ((binCoord[0] < DataStruc_qmin) || (binCoord[0] > DataStruc_qmax) ||
              (binCoord[1] < DataStruc_qmin) || (binCoord[1] > DataStruc_qmax) ||
              (binCoord[2] < DataStruc_qmin) || (binCoord[2] > DataStruc_qmax)) {
            continue;
          }
          // redefine reject, should be false if bool reject_final is false, otherwise equal to reject
          reject = reject_final ? reject : false;
          // if (reject) continue;
          if (reject) {
            binError *= 100;
            // make bincontent of inHistRatio zero
            inHistRatio->SetBinContent(ix, iy, iz, inHistRatioMax);
            inHistRatio->SetBinError(ix, iy, iz, 0);
          }
          data.Add(binCoord, binContent, binError); //binCoordErr,
          // if (i==3)
          //   cout<<inHistRatio->GetXaxis()->GetBinCenter(ix)<<"  "<<inHistRatio->GetYaxis()->GetBinCenter(iy)<<"  "<<inHistRatio->GetZaxis()->GetBinCenter(iz)<<"\t"<<binContent<<"  "<<binError<<endl;
          // cout lower part of the bins and also the higher part of the bins
          // cout<<inHistRatio->GetXaxis()->GetBinLowEdge(ix)<<"  "<<inHistRatio->GetXaxis()->GetBinUpEdge(ix)<<"\t"
          //     <<inHistRatio->GetYaxis()->GetBinLowEdge(iy)<<"  "<<inHistRatio->GetYaxis()->GetBinUpEdge(iy)<<"\t"
          //     <<inHistRatio->GetZaxis()->GetBinLowEdge(iz)<<"  "<<inHistRatio->GetZaxis()->GetBinUpEdge(iz)<<endl;
        }
      }
    }
    inHistRatio->SetMinimum(inHistRatioMin);
    //set max range of inHistRatio to qMax
    inHistRatio->GetXaxis()->SetRangeUser(qMin3D[0], qMin3D[hist_number] - qMin3D[0]);
    inHistRatio->GetYaxis()->SetRangeUser(qMin3D[0], qMin3D[hist_number] - qMin3D[0]);
    inHistRatio->GetZaxis()->SetRangeUser(qMin3D[0], qMin3D[hist_number] - qMin3D[0]);
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
        c_individual->SetRightMargin(0.12); // Increase right margin to make space for the palette
        histtemp->SetMarkerStyle(20);
        histtemp->SetMarkerSize(1);
        histtemp->SetMarkerColor(kRed);
        histtemp->SetLineColor(kRed);
        histtemp->SetLineWidth(1);
        histtemp->SetFillColorAlpha(kRed, 0.5);
        histtemp->SetFillStyle(3002);
        histtemp->SetMinimum(0.005);
        histtemp->SetMaximum(0.04);
        configure3DAxis(histtemp.get());
        // histtemp->SetMinimum(avgChi/2);
        // histtemp->SetMaximum(avgChi*4);
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



void MakeHistogram() {
    cout << "INFO: Creating histogram..." << endl;
    string c2_type = "data";
    std::unique_ptr<TH3D> inHistData;
    std::unique_ptr<TH3D> inHistMC;
    std::vector<std::unique_ptr<TH3D>> inHistRatio;
    string histtemp;
    std::vector<std::unique_ptr<TH3D>> inHistVec;
    std::vector<std::unique_ptr<TH3D>> inHistDataVec;
    std::vector<std::unique_ptr<TH3D>> inHistMCVec;
    string sign = "ppmm";
    string sign2 = "pm";
    // std::vector<double> qBins3D_width = {0.04, 0.12, 0.6};
    std::vector<double> qMin3D_4_0 = {.0, .28, .6, 1.2, 4.2};
    std::vector<double> qMin3D_4_02 = {0.02, .28, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_0 = {.0, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_04 = {0.04, .6, 1.2, 4.2};
    std::vector<double> qMin3D_3_02 = {0.02, .62, 1.22, qMax};
    double integral, integral2;

    qMin3D_0 = qMin3D_3_0;
    qMin3D_02 = qMin3D_3_02;
    const int hist_number = qMin3D_02.size()-1;
    // qMin3D_02 = qMin3D_3_0;
    //redefine varqbins according to qMin3D_02 and qbins3D_width
    // varQbins.clear();
    // for (int i = 0; i < qMin3D_02.size(); i++) {
    //   for 

    if (!fileData.empty() && !histData.empty()) {
      std::unique_ptr<TFile> inFile = std::make_unique<TFile>(fileData.c_str(), "READ");
      if (!inFile->IsOpen()) {
        throw "ERROR: Data file not found!";
      }
      // std::unique_ptr<TH3D> inHist(dynamic_cast<TH3D*>(inFile->Get(histData.c_str())));
      // if (!inHist) {
      //   throw "ERROR: Data histogram not found!";
      // }
      //check for the 4th histogram, which is always needed for fitting
      // histtemp = histData;
      // size_t pos = histtemp.find(sign);
      // histtemp.insert(pos, to_string(4) + "_");
      // std::unique_ptr<TH3D> inHist(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
      // inHist->SetDirectory(0);
      for (int i = 4-hist_number+1; i <= 4; i++) {
        histtemp = histData;
        size_t pos = histtemp.find(sign);
        histtemp.insert(pos, to_string(i) + "_");
        cout<<"histtemp: "<<histtemp<<endl;
        std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
        inHisttemp->SetDirectory(0);
        inHistVec.push_back(std::move(inHisttemp));
      }
      if (!inHistVec[0]) {
        throw "ERROR: Data histogram not found!";
      }
      if (!hist2Data.empty()) {
        std::unique_ptr<TH3D> inHist2;
        std::vector<std::unique_ptr<TH3D>> inHistVec2;
        if (file2Data.empty()) {
          // inHist2 = std::make_unique<TH3D>(*dynamic_cast<TH3D*>(inFile->Get(hist2Data.c_str())));
          // if (!inHist2) {
          //   throw "ERROR: Second data histogram not found!";
          // }
          // check for the 4th histogram, which is always needed for fitting
          // histtemp = hist2Data;
          // size_t pos = histtemp.find(sign2);
          // histtemp.insert(pos, to_string(4) + "_");
          // inHist2 = std::make_unique<TH3D>(*dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
          // inHist2->SetDirectory(0);
          for (int i = 4-hist_number+1; i <= 4; i++) {
            histtemp = hist2Data;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          if (!inHistVec2[0]) {
            throw "ERROR: Second data histogram not found!";
          }
        } else {
          std::unique_ptr<TFile> inFile2 = std::make_unique<TFile>(file2Data.c_str(), "READ");
          if (!inFile2->IsOpen()) {
            throw "ERROR: Second data file not found!";
          }
          // inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile2->Get(hist2Data.c_str())));
          // if (!inHist2) {
          //   throw "ERROR: Second data histogram not found!";
          // }
          // inHist2->SetDirectory(0);
          for (int i = 4-hist_number+1; i <= 4; i++) {
            histtemp = hist2Data;
            size_t pos = histtemp.find(sign); //works for ohp, mixing, we use ppmm (or sign)
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile2->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          if (!inHistVec2[0]) {
            throw "ERROR: Second data histogram not found!";
          }
          inFile2->Close();
        }
        for (int i = 0; i < hist_number; i++) {
          if (!inHistVec2[i]) {
            throw "ERROR: Second data histogram not found!";
          }
          cout<<"integral data: "<<inHistVec[i]->Integral()<<" integral data ref: "<<inHistVec2[i]->Integral()<<endl;
          integral = inHistVec[inHistVec.size()-1]->Integral(1, 5, 1, 3, 1, inHistVec[inHistVec.size()-1]->GetNbinsZ());
          std::cout<<"integral data: "<<integral<<"integral all data: "<<inHistVec[inHistVec.size()-1]->Integral()<<std::endl;
          integral2 = inHistVec2[inHistVec2.size()-1]->Integral(1, 5, 1, 3, 1, inHistVec2[inHistVec2.size()-1]->GetNbinsZ());
          inHistVec[i]->Divide(inHistVec[i].get(), inHistVec2[i].get(),
                              1./integral, 1./integral2);
                              // 1./inHistVec[inHistVec.size()-1]->Integral(), 1./inHistVec2[inHistVec2.size()-1]->Integral());
                              // 1./inHist->Integral(), 1./inHist2->Integral());
          cout<<"integral data_ratio: "<<inHistVec[i]->Integral()<<endl;
          string histName = inHistVec[i]->GetName();
          histName += "_over_";
          histName += inHistVec2[i]->GetName();
          inHistVec[i]->SetName(histName.c_str());
        }
        // cout<<"before inHist divide"<<endl;
        // inHist->Divide(inHist.get(), inHist2.get(),
        //                1./inHist->Integral(), 1./inHist2->Integral());
        // cout<<"after inHist divide"<<endl;
        // string histName = inHist->GetName();
        // histName += "_over_";
        // histName += inHist2->GetName();
        // inHist->SetName(histName.c_str());
        // inHist2.reset();
        // inHistVec2.clear();
      }

      // inHistData = std::make_unique<TH3D>(*inHist);
      // inHistData = std::move(inHist);
      // if (verbose) inHistData->Print();
      if (verbose) inHistVec[0]->Print();
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
      // std::unique_ptr<TH3D> inHist(dynamic_cast<TH3D*>(inFile->Get(histMC.c_str())));
      // if (!inHist) {
      //   throw "ERROR: MC histogram not found!";
      // }
      // inHist->SetDirectory(0);
      for (int i = 4-hist_number+1; i <= 4; i++) {
        histtemp = histMC;
        size_t pos = histtemp.find(sign);
        histtemp.insert(pos, to_string(i) + "_");
        cout<<"histtemp: "<<histtemp<<endl;
        std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
        inHisttemp->SetDirectory(0);
        inHistVec.push_back(std::move(inHisttemp));
      }
      if (!inHistVec[0]) {
        throw "ERROR: MC histogram not found!";
      }
      if (!hist2MC.empty()) {
        // std::unique_ptr<TH3D> inHist2;
        std::vector<std::unique_ptr<TH3D>> inHistVec2;
        if (file2MC.empty()) {
          // inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile->Get(hist2MC.c_str())));
          // if (!inHist2) {
          //   throw "ERROR: Second data histogram not found!";
          // }
          for (int i = 4-hist_number+1; i <= 4; i++) {
            histtemp = hist2MC;
            size_t pos = histtemp.find(sign2);
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          if (!inHistVec2[0]) {
            throw "ERROR: Second data histogram not found!";
          }
        } else {
          std::unique_ptr<TFile> inFile2 = std::make_unique<TFile>(file2MC.c_str(), "READ");
          if (!inFile2->IsOpen()) {
            throw "ERROR: Second data file not found!";
          }
          // inHist2 = std::unique_ptr<TH3D>(dynamic_cast<TH3D*>(inFile2->Get(hist2MC.c_str())));
          // if (!inHist2) {
          //   throw "ERROR: Second data histogram not found!";
          // }
          // inHist2->SetDirectory(0);
          for (int i = 4-hist_number+1; i <= 4; i++) {
            histtemp = hist2MC;
            size_t pos = histtemp.find(sign); //works for ohp, mixing, we use ppmm (or sign)
            histtemp.insert(pos, to_string(i) + "_");
            std::unique_ptr<TH3D> inHisttemp(dynamic_cast<TH3D*>(inFile2->Get(histtemp.c_str())));
            inHisttemp->SetDirectory(0);
            inHistVec2.push_back(std::move(inHisttemp));
          }
          if (!inHistVec2[0]) {
            throw "ERROR: Second data histogram not found!";
          }
          inFile2->Close();
        }
        for (int i = 0; i < hist_number; i++) {
          cout<<"integral mc: "<<inHistVec[i]->Integral()<<" integral mc ref: "<<inHistVec2[i]->Integral()<<endl;
          integral = inHistVec[inHistVec.size()-1]->Integral(1, 5, 1, 3, 1, inHistVec[inHistVec.size()-1]->GetNbinsZ());
          integral2 = inHistVec2[inHistVec2.size()-1]->Integral(1, 5, 1, 3, 1, inHistVec2[inHistVec2.size()-1]->GetNbinsZ());
          inHistVec[i]->Divide(inHistVec[i].get(), inHistVec2[i].get(),
                              1./integral, 1./integral2);
                              // 1./inHistVec[inHistVec.size()-1]->Integral(), 1./inHistVec2[inHistVec2.size()-1]->Integral());
                              // 1./inHist->Integral(), 1./inHist2->Integral());
          cout<<"integral mc_ratio: "<<inHistVec[i]->Integral()<<endl;
          string histName = inHistVec[i]->GetName();
          histName += "_over_";
          histName += inHistVec2[i]->GetName();
          inHistVec[i]->SetName(histName.c_str());
        }
        inHistVec2.clear();
        // inHist->Divide(inHist.get(), inHist2.get(),
        //                1./inHist->Integral(), 1./inHist2->Integral());

        // string histName = inHist->GetName();
        // histName += "_over_";
        // histName += inHist2->GetName();
        // inHist->SetName(histName.c_str());
        // inHist2.reset();
      }

      // inHistMC = std::make_unique<TH3D>(*inHist);
      // inHistMC = std::move(inHist);
      // if (verbose) inHistMC->Print();
      if (verbose) inHistVec[0]->Print();
      for (auto& inHisttemp : inHistVec) {
        inHistMCVec.push_back(std::move(inHisttemp));
      }
      inHistVec.clear();
      inFile->Close();
    }

    // if (inHistData && inHistMC) {
    if (!inHistDataVec.empty() && !inHistMCVec.empty()) {
      // // string histName = inHistData->GetName();
      // string histName = inHistDataVec[0]->GetName();
      // histName += "_over_";
      // // histName += inHistMC->GetName();
      // histName += inHistMCVec[0]->GetName();
      // inHistRatio.push_back(std::make_unique<TH3D>(*inHistData));
      // inHistRatio.push_back(std::make_unique<TH3D>(*inHistDataVec[0]));
      // if (!hist2Data.empty() || !hist2MC.empty()) {
      //   inHistRatio[0]->Divide(inHistData.get(), inHistMC.get());
      // } else {
      //   inHistRatio[0]->Divide(inHistData.get(), inHistMC.get(),
      //                       1./inHistData->Integral(), 1./inHistMC->Integral());
      // }
      hist = new TH3D(*(inHistDataVec[1].get()));
      hist->Divide(inHistDataVec[1].get(), inHistMCVec[1].get());
      histLarge = new TH3D(*(inHistDataVec[2].get()));
      histLarge->Divide(inHistDataVec[2].get(), inHistMCVec[2].get());
                // 1. / inHistDataVec[0]->Integral(), 1. / inHistMCVec[0]->Integral());

      // inHistData.reset();
      // inHistMC.reset();
      // hist = inHistRatio[0].get();
      // hist = new TH3D(*(inHistRatio[0].get()));
      //clear InHistRatio
      // inHistRatio.clear();
      // HisttoData(inHistDataVec, inHistMCVec, inHistRatio, DataStruc);
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      std::cout << "DataStruc size: " << DataStruc.Size() << std::endl;
      
      r2 = true;
    // } else if (!inHistData && inHistMC) {
    } else if (inHistDataVec.empty() && !inHistMCVec.empty()) {
      // hist = new TH3D(*(inHistMC.get()));
      hist = new TH3D(*(inHistMCVec[1].get()));
      histLarge = new TH3D(*(inHistMCVec[2].get()));
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      r2 = false;
      c2_type = "mc";
    // } else if (inHistData && !inHistMC) {
    } else if (!inHistDataVec.empty() && inHistMCVec.empty()) {
      // hist = new TH3D(*(inHistData.get()));
      hist = new TH3D(*(inHistDataVec[1].get()));
      histLarge = new TH3D(*(inHistDataVec[2].get()));
      // HisttoData(inHistDataVec, inHistMCVec, DataStruc_0, qMin3D_0);
      HisttoData(inHistDataVec, inHistMCVec, DataStruc, qMin3D_02);
      r2 = false;
      c2_type = "data";
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
    cout<<"projMax "<<projMax<<endl;
    //get binwidth of the histogram hist
    double binWidth = hist->GetXaxis()->GetBinWidth(1);
    if (hist_number == 3) {
      projMax = int(projMax/2*0.04/binWidth);
    } else if (hist_number == 2) {
      projMax = int(projMax/6*0.12/binWidth);
    }
    cout<<"projMax after change"<<projMax<<endl;

    // if (zMin > hist->GetMinimum()) {
    //   zMin = hist->GetMinimum() + 0.02;
    //   // yMin = hist->GetMinimum() + 0.02;
    // }
    // if (zMax < hist->GetMaximum()) {
    //   zMax = hist->GetMaximum() + 0.02;
    //   yMax = hist->GetMaximum() + 0.02;
    // }
    // hist->SetMinimum(zMin);
    // hist->SetMaximum(zMax);
    double zMinFirst3Bins = zMin;
    double zMaxFirst3Bins = zMax;

    for (int ix = 1; ix <= 3; ++ix) { // Loop over the first 3 bins in X
        for (int iy = 1; iy <= 3; ++iy) { // Loop over the first 3 bins in Y
            for (int iz = 1; iz <= 3; ++iz) { // Loop over the first 3 bins in Z
                double binContent = hist->GetBinContent(ix, iy, iz);
                if (binContent < zMinFirst3Bins) {
                    zMinFirst3Bins = binContent;
                }
                if (binContent > zMaxFirst3Bins) {
                    zMaxFirst3Bins = binContent;
                }
                // cout<<"binContent "<<binContent<<" zMinFirst3Bins "<<zMinFirst3Bins<<" zMaxFirst3Bins "<<zMaxFirst3Bins<<endl;
            }
        }
    }

    // Update zMin and zMax based on the first 3 bins
    if (projMax < 5) {
      if (zMin > zMinFirst3Bins) {
          zMin = zMinFirst3Bins + 0.02;
          // yMin = zMinFirst3Bins + 0.02; // Uncomment if yMin needs to be updated similarly
      }
      if (zMax < zMaxFirst3Bins) {
          zMax = zMaxFirst3Bins + 0.02;
          yMax = zMaxFirst3Bins + 0.02; // Update yMax as well
      }
    }

    if (r2) {
      c2_func = "r2_";
    } else {
      c2_func = "c2_" + c2_type + "_";
    }

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

void Fit() {
  cout << "INFO: Fitting..." << endl;
  std::vector<int> qNBins3D = {14, 10, 7, 0};
  std::vector<double> qMin3D = qMin3D_02;
  // Calculate the number of bins to remove from the second-to-last segment
  double lastEdge = 4.22;
  double binWidthlast = 0.6;
  int qNBins3D_lessbins = static_cast<int>((lastEdge - qMin3D[qMin3D.size() - 1] + 0.001) / binWidthlast);

  // Subtract the calculated number of bins from the second-to-last segment
  qNBins3D[qNBins3D.size() - 2] -= qNBins3D_lessbins;

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
    case 7: c2 = new C2_7();
            break;
    case 8: c2 = new C2_8();
            break;
    case 9: c2 = new C2_9();
            break;
    case 10: c2 = new C2_10();
            break;
    default: throw "ERROR: C2 function index not found!";
             break;
  }

  c2->rejFrom = rejFrom;
  // c2->rejFromOsl = rejFromOsl;
  c2->rejTo = rejTo;
  // c2->rejToOsl = rejToOsl;
  c2->rej2From = rej2From;
  c2->rej2To = rej2To;
  c2->alpha = alpha;
  c2->fixC0 = fixC0;
  c2->useEps = useEps;
  c2->with_errors = with_errors;
  nParams = c2->nParams;
  cout<<"nParams "<<nParams<<endl;
  
  func = new TF3(("C2_" + std::to_string(c2Index)).c_str(), c2,
                 qMin, qMax,
                 qMin, qMax,
                 qMin, qMax,
                 c2->nParams);
  c2->setup(func);
  funcInit = new TF3(("C2_" + std::to_string(c2Index) + "_init").c_str(), c2,
                 qMin, qMax,
                 qMin, qMax,
                 qMin, qMax,
                 c2->nParams);
  c2->setup(funcInit);

  ROOT::Fit::DataOptions DataStruc_options;
  ROOT::Fit::DataRange DataStruc_range(3);

  DataStruc_options.fUseRange = true;
  double DataStruc_qmin = qMin > 0 ? qMin : 0;
  double DataStruc_qmax = qMax > 0 ? qMax : 4.2;
  // DataStruc_range.SetRange(DataStruc_qmin, DataStruc_qmax, DataStruc_qmin, DataStruc_qmax, DataStruc_qmin, DataStruc_qmax);
  DataStruc_range.SetRange(0, DataStruc_qmin, DataStruc_qmax); // X range
  DataStruc_range.SetRange(1, DataStruc_qmin, DataStruc_qmax); // Y range
  DataStruc_range.SetRange(2, DataStruc_qmin, DataStruc_qmax); // Z range
  cout << "DataStruc_qmin: " << DataStruc_qmin << " DataStruc_qmax: " << DataStruc_qmax << endl;


  ROOT::Fit::BinData DataStruc_copy = ROOT::Fit::BinData (DataStruc_options, DataStruc_range, 100000, 3);
  // auto DataStruc_copy = std::make_shared<ROOT::Fit::BinData>(DataStruc_options, DataStruc_range, 100000, 3);
  for (size_t i = 0; i < DataStruc.Size(); ++i) {
    double value, invError;
    const double* coords = DataStruc.GetPoint(i, value, invError);
    DataStruc_copy.Add(coords, value, 1.0/invError);
    // DataStruc_copy->Add(coords, value, 1.0/invError);
  }

  std::vector<std::tuple<std::vector<double>, double, double>> DataStruc_points;
    for (size_t i = 0; i < DataStruc.Size(); ++i) {
      double value, invError;
      const double* coords = DataStruc.GetPoint(i, value, invError);
      std::vector<double> coordinates(coords, coords + DataStruc.NDim());
      DataStruc_points.emplace_back(coordinates, value, 1.0/invError);
    }
  ROOT::Fit::Fitter fitterMigradInit, fitterMigrad;
  ROOT::Math::WrappedMultiTF1 wrappedFunc(*func, 3);
  std::cout << "wrappedFunc Dimensions: " << wrappedFunc.NDim() << std::endl;
  std::cout << "wrappedFunc Parameters: " << wrappedFunc.NPar() << std::endl;
  // make fit with Migrad initially
  fitterMigradInit.SetFunction(wrappedFunc);
  fitterMigradInit.Config().MinimizerOptions().SetMinimizerAlgorithm("Migrad");
  fitterMigradInit.Config().MinimizerOptions().SetMaxIterations(10000);
  fitterMigradInit.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
  // std::vector<double> defaultParams(nParams, 0.0);
  for (int i = 0; i < nParams; i++) {
    double min, max;
    func->GetParLimits(i, min, max);
    fitterMigradInit.Config().ParSettings(i).SetLimits(min, max);
    // defaultParams[i] = func->GetParameter(i);
  }
  std::vector<double> initialParams(nParams, 0.0);
  bool useInit = false;
  bool retMigradInit = false;
  if (useInit) {
    retMigradInit = fitterMigradInit.Fit(DataStruc_copy);
    const ROOT::Fit::FitResult & resMigradInit = fitterMigradInit.Result();
    resMigradInit.Print(std::cout);
    funcInit->SetFitResult(resMigradInit);
    for (int i = 0; i < nParams; ++i) {
      initialParams[i] = resMigradInit.Parameter(i);
    }
  } else {
    // Set default parameters for the initial fit
    for (int i = 0; i < nParams; ++i) {
      initialParams[i] = func->GetParameter(i);
    }
  }
  // Retrieve the fitted parameters from the initial fit


  // Print the fitted parameters from the initial fit
  std::cout << "Fitted parameters from initial fit: ";
  for (const auto& param : initialParams) {
      std::cout << param << " ";
  }
  std::cout << std::endl;
  std::vector<double> params = initialParams;
  // use default params for now
  // std::vector<double> params = defaultParams;
  // CustomChi2FCN chi2Fcn(DataStruc_copy, wrappedFunc);
  CustomChi2FCN chi2Fcn(DataStruc_copy, *func, qMin3D, qNBins3D);
  bool annealing = false;
  ROOT::Fit::Fitter fitterSimulatedAnnealing;
  if (!mean_fit) {
    ROOT::Math::WrappedMultiTF1 wf(*func, 3);
    fitterMigrad.SetFunction(wf);
    fitterMigrad.Config().MinimizerOptions().SetMinimizerAlgorithm("Migrad");
    // fitterMigrad.Config().SetMinimizer("Minuit","Migrad");

    if (annealing) {
      ROOT::Math::WrappedMultiTF1 wf2(*func, 3);
      fitterSimulatedAnnealing.SetFunction(wf2);
      fitterSimulatedAnnealing.Config().MinimizerOptions().SetMinimizerAlgorithm("SimulatedAnnealing");
      // fitterSimulatedAnnealing.Config().MinimizerOptions().SetMaxIterations(10000);
      // fitterSimulatedAnnealing.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
    }
    fitterMigrad.Config().MinimizerOptions().SetMaxIterations(10000);
    fitterMigrad.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
    // fitterSimulatedAnnealing.Config().MinimizerOptions().SetPrintLevel(2);
  } else {
        // Ensure DataStruc is correctly initialized
    std::cout << "DataStruc size: " << DataStruc.Size() << std::endl;
    if (DataStruc.Size() == 0) {
        throw std::runtime_error("DataStruc is not initialized or empty");
    }
    // for (int i = 0; i < nParams; ++i) {
    //     params[i] = func->GetParameter(i);
    // }

    // Print initialized parameters
    std::cout << "Initialized params: ";
    for (const auto& param : params) {
        std::cout << param << " ";
    }
    std::cout << std::endl;
    // Create the custom chi-squared function
    // auto chi_now = ComputeChi2(DataStruc_points, *func, qMin3D, qNBins3D);
    // std::cout << "Initial chi2: " << chi_now << std::endl;
    // CustomChi2FCN chi2Fcn(DataStruc_copy, *func);//, qMin3D, qNBins3D);
    std::cout << "CustomChi2FCN object created" << std::endl;
    // fitterMigrad.Config().SetMinimizer("Minuit2", "Migrad");
    fitterMigrad.Config().MinimizerOptions().SetMinimizerAlgorithm("Migrad");
    fitterMigrad.Config().MinimizerOptions().SetMaxIterations(10000);
    fitterMigrad.Config().MinimizerOptions().SetMaxFunctionCalls(10000);
    // fitterMigrad.Config().SetFunction(chi2Fcn, func->GetNpar());
    // fitterMigrad.SetFCN(nParams, chi2Fcn, nullptr, 0, true);
    // fitterMigrad.SetFCN(nParams, chi2Fcn, params.data(), 0, true);
    // fitterMigrad.SetFCN(chi2Fcn, params.data(), DataStruc_copy.Size(), 0);
    fitterMigrad.SetFCN(chi2Fcn, params.data(), DataStruc_copy.Size(), true);
    // fitterMigrad.SetFCN(chi2Fcn, params.data(), DataStruc_copy->Size(), 0);
    std::cout << "FCN set" << std::endl;

  }
  const int first_bins = 5;
  auto chi_first = ComputeChi2ForFirstBins(DataStruc_points, *func, qMin3D, qNBins3D, first_bins);
  chi_firstBins.first = chi_first;
  cout<< "Chi2 for "<< first_bins << " bins: "<< chi_firstBins.first << endl;

  for (int i = 0; i < nParams; i++) {
    double min, max;
    func->GetParLimits(i, min, max);
    fitterMigrad.Config().ParSettings(i).SetLimits(min, max);
    if (annealing) {
      fitterSimulatedAnnealing.Config().ParSettings(i).SetLimits(min, max);
    }
  }
  bool retSimulatedAnnealing = false;
  ROOT::Fit::FitResult resSimulatedAnnealing; 
  if (annealing) {
    retSimulatedAnnealing = fitterSimulatedAnnealing.Fit(DataStruc);
    resSimulatedAnnealing = fitterSimulatedAnnealing.Result();
  }
  // bool retMigrad = fitterMigrad.Fit(DataStruc);
  // const ROOT::Fit::FitResult & resMigrad = fitterMigrad.Result();

  // bool ret = false;
  // const ROOT::Fit::FitResult* res = nullptr;

  bool ret = false;
  bool retMigrad = false;
  const ROOT::Fit::FitResult* res = nullptr;

  std::cout << "Starting fit..." << std::endl;
  if (mean_fit) {
    for (int i = 0; i < nParams; i++) {
        auto settings = fitterMigrad.Config().ParSettings(i);
        std::cout << "Parameter " << i << ": "
                  << "Lower Limit = " << settings.LowerLimit() 
                  << ", Upper Limit = " << settings.UpperLimit() 
                  << std::endl;
    }
      retMigrad = fitterMigrad.FitFCN();
  } else {
      retMigrad = fitterMigrad.Fit(DataStruc);
  }
  std::cout << "Fit completed." << std::endl;

  double chi2_before = chi2Fcn(initialParams.data());
  double chi2_after = fitterMigrad.Result().MinFcnValue();
  std::cout << "Chi2 before: " << chi2_before <<" Chi2 after: " << chi2_after <<" Chi2 gain: " << chi2_before - chi2_after << std::endl;
  const ROOT::Fit::FitResult& resMigrad = fitterMigrad.Result();

  // Compare the results
  // if (retMigrad || retSimulatedAnnealing) {
  //   std::cout<<"Migrad chi2 "<<resMigrad.Chi2()<<" Simulated Annealing chi2 "<<resSimulatedAnnealing.Chi2()<<std::endl;
  //   if (resMigrad.Chi2() == resSimulatedAnnealing.Chi2()) {
  //     std::cout << "Both fits are equally good" << std::endl;
  //     ret = retMigrad;
  //     res = &resMigrad;
  //   } else {
  //     std::cout<<"THERE IS DIFFERENCE IN CHI2: "<<resMigrad.Chi2()-resSimulatedAnnealing.Chi2()<<std::endl;
  //     if (resMigrad.Chi2() < resSimulatedAnnealing.Chi2()) {
  //       std::cout << "Migrad fit is better" << std::endl;
  //       ret = retMigrad;
  //       res = &resMigrad;
  //     } else {
  //       std::cout << "Simulated Annealing fit is better" << std::endl;
  //       ret = retSimulatedAnnealing;
  //       res = &resSimulatedAnnealing;
  //     }

  //   }
  if (retMigrad) {
    std::cout << "Migrad fit is successfil" << std::endl;
    ret = retMigrad;
    res = &resMigrad;
  } else {
    std::cout << "Migrad fit failed" << std::endl;
  }
  if (annealing) {
    if (retSimulatedAnnealing) {
      std::cout << "Simulated Annealing fit is successful" << std::endl;
      ret = retSimulatedAnnealing;
      res = &resSimulatedAnnealing;
    } else {
      std::cout << "Simulated Annealing fit failed" << std::endl;
    }
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

    // gMinuit->mnmatu(1);
    // TFitResult t_res (res);



    //matrix print
    TMatrixDSym correlationMatrix(res->NPar());
    res->GetCorrelationMatrix(correlationMatrix);
    // //skip parameters with fixed values
    std::vector<int> fixedIndices;
    std::vector<std::string> parameterNames;
    for (int i = 0; i < res->NPar(); ++i) {
      if (fitterMigrad.Config().ParSettings(i).IsFixed()) {
        fixedIndices.push_back(i);
      } else {
        std::string parname = fitterMigrad.Config().ParSettings(i).Name();
        //check lenght of parname
        if (parname.length() < 8) {
          parname = parname + "    ";
        }
        parameterNames.push_back(parname);
      }
    }

    TMatrixDSym tempMatrix(res->NPar() - fixedIndices.size());
    int ii = 0;
    for (int row = 0; row < correlationMatrix.GetNrows(); ++row) {
      if (std::find(fixedIndices.begin(), fixedIndices.end(), row) != fixedIndices.end()) continue;
      int jj = 0;
      for (int col = 0; col < correlationMatrix.GetNcols(); ++col) {
        if (std::find(fixedIndices.begin(), fixedIndices.end(), col) != fixedIndices.end()) continue;
        tempMatrix(ii, jj) = correlationMatrix(row, col);
        ++jj;
      }
      ++ii;
    }
    correlationMatrix.ResizeTo(tempMatrix.GetNrows(), tempMatrix.GetNcols());
    correlationMatrix = tempMatrix;
    std::cout << "\t";
    for (int i = 0; i < correlationMatrix.GetNcols(); ++i) {
      std::cout << "\t" << parameterNames[i];
    }
    std::cout << std::endl;
    for (int i = 0; i < correlationMatrix.GetNrows(); ++i) {
      std::cout << parameterNames[i] << "\t";
      for (int j = 0; j < correlationMatrix.GetNcols(); ++j) {
          std::cout << std::fixed << std::setprecision(4) << correlationMatrix(i, j) << "\t\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    TMatrixDSymEigen eigensolver(correlationMatrix);
    TVectorD eigenValues = eigensolver.GetEigenValues();
    TMatrixD eigenVectors = eigensolver.GetEigenVectors();

    // Print the parameter names
    std::cout << "eigenval";
    for (int i = 0; i < correlationMatrix.GetNcols(); ++i) {
      std::cout << "\t" << parameterNames[i];
    }
    std::cout << std::endl;

    // Print out the eigenvalues and eigenvectors
    for (int i = 0; i < eigenValues.GetNrows(); ++i) {
      std::cout << i << ": " << eigenValues(i);
      for (int j = 0; j < correlationMatrix.GetNcols(); ++j) {
        std::cout << "\t" << eigenVectors(i, j) << "\t";
      }
      std::cout << std::endl;
    }
    // correlationMatrix.Print();
    fitData.parameterNames = parameterNames;
    fitData.correlationMatrix.ResizeTo(correlationMatrix.GetNrows(), correlationMatrix.GetNcols());
    fitData.eigenValues.ResizeTo(eigenValues.GetNrows());
    fitData.eigenVectors.ResizeTo(eigenVectors.GetNrows(), eigenVectors.GetNcols());

    fitData.correlationMatrix = correlationMatrix;
    fitData.eigenValues = eigenValues;
    fitData.eigenVectors = eigenVectors;
  }
  else
    Error("Fit", "3D fit failed");
}


void Output() {
  cout << "INFO: Saving..." << endl;

  outFilePathBase = "output/" + Today() + "/" + c2_func + std::to_string(c2Index);
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

  string titleR2 = "";
  if (r2) {
    titleR2 = "R_{2}";
  } else {
    titleR2 = "C_{2}";
  }

  // TFile* outFile = new TFile(outFilePathRoot.c_str(), "RECREATE");
  std::unique_ptr<TFile> outFile(TFile::Open(outFilePathRoot.c_str(), "RECREATE"));
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Error: Failed to open file " << outFilePathRoot << std::endl;
    return;
  }
  cout<<"outFilePathRoot "<<outFilePathRoot<<"plotName "<<plotName<<endl;
  hist->Write(plotName.c_str());
  // func->Write(("fitFunc_" + plotName).c_str());
  //save InHistRatio histograms
  gStyle->SetOptStat(1100); 
  for (size_t i = 0; i < InHistRatio.size(); ++i) {
    configure3DAxis(InHistRatio[i].get(), true);
    auto canvas = std::make_unique<TCanvas>(("canvas_" + std::to_string(i)).c_str(), ("canvas_" + std::to_string(i)).c_str(), 800, 600);
    // InHistRatio[i]->Write(("InHistRatio_" + std::to_string(i) + "_" + plotName).c_str());
    canvas->SetRightMargin(0.12); // Increase right margin to make space for the palette
    InHistRatio[i]->Draw("box2z");
    canvas->Write(("InHistRatio_" + std::to_string(i) + "_" + plotName).c_str());
  }

  std::map<std::tuple<double, double, double>, int> coordinates;

  for (unsigned int i = 0; i < DataStruc.Size(); ++i) {
      double value;
      const double* coords = DataStruc.GetPoint(i, value);

      coordinates[std::make_tuple(coords[0], coords[1], coords[2])] = i;
      // cout<<coords[0]<<"  "<<coords[1]<<"  "<<coords[2]<<"bin "<<i<<endl;
  }
  // copy hist to new hist
  std::vector<double> qMin3D = qMin3D_02;
  // std::vector<int> qNBins3D = {14, 10, 7, 0};
  std::vector<bool> reject_region = rej_OSL();
  std::vector<int> qNBins3D = {14, 15, 10, 7, 0};

  // Calculate the number of bins to remove from the second-to-last segment
  double lastEdge = 4.22;
  double binWidthlast = 0.6;
  int qNBins3D_lessbins = static_cast<int>((lastEdge - qMin3D[qMin3D.size() - 1] + 0.001) / binWidthlast);

  // Subtract the calculated number of bins from the second-to-last segment
  qNBins3D[qNBins3D.size() - 2] -= qNBins3D_lessbins;

  //reduce qNBins3D from the left to fit the size of qMin3D
  if (qNBins3D.size() > qMin3D.size()) {
    // Calculate the number of elements to remove from the beginning
    size_t elementsToRemove = qNBins3D.size() - qMin3D.size();
    // Remove elements from the beginning
    qNBins3D.erase(qNBins3D.begin(), qNBins3D.begin() + elementsToRemove);
  }
  //cout qnbins3d contents
  for (int i = 0; i < qNBins3D.size(); i++) {
    cout<<qNBins3D[i]<<" ";
  }
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

  int nBins = 0;
  vector<double> binEdges;
  // binEdges.push_back(0);
  for (int i = 0; i < qMin3D.size()-1; i++) {
    int bins = qNBins3D[i]*((qMin3D[i+1]-qMin3D[i])/(qMin3D[i+1]-qMin3D[0]));
    nBins += bins;
    double binWidth = (qMin3D[i+1]-qMin3D[0]) / qNBins3D[i];
    // std::cout<<"bins "<<bins<<" binWidth "<<binWidth<<std::endl;
    for (int j = 0; j < bins; j++) {
      binEdges.push_back(qMin3D[i] + j*binWidth);
      std::cout<<qMin3D[i] + j*binWidth<<" ";
    }
  }
  binEdges.push_back(qMin3D[qMin3D.size() - 1]);
  cout<<"number of bins: "<<binEdges.size()-1<<endl;
  //cout binEdges content
  for (int i = 0; i < binEdges.size()+1; i++) {
    cout<<binEdges[i]<<" ";
  }
  int Q_values = 10000;
  // std::unique_ptr<TH1D> hist_Q_1D(new TH1D("hist_Q_1D", "hist_Q_1D", nBins, newEdges.data()));
  // std::unique_ptr<TH1D> hist_Q_1D(new TH1D("hist_Q_1D", "hist_Q_1D", Q_values/50, 0, 4.2));
  // std::unique_ptr<TH1D> hist_Q_1D_p(new TH1D("hist_Q_1D_p", "hist_Q_1D_p", Q_values/50, 0, 4.2));
  // std::unique_ptr<TH1D> hist_Q_1D_rej(new TH1D("hist_Q_1D_rej", "hist_Q_1D_rej", Q_values/50, 0, 4.2));
  // std::unique_ptr<TH1D> hist_Q_1D_norej(new TH1D("hist_Q_1D_norej", "hist_Q_1D_norej", Q_values/50, 0, 4.2));
  auto hist_Q_1D = std::make_unique<TH1D>("hist_Q_1D", "hist_Q_1D", binEdges.size()-1, binEdges.data());
  auto hist_Q_1D_p = std::make_unique<TH1D>("hist_Q_1D_p", "hist_Q_1D_p", binEdges.size()-1, binEdges.data());
  auto hist_Q_1D_diff = std::make_unique<TH1D>("hist_Q_1D_diff", "hist_Q_1D_diff", binEdges.size()-1, binEdges.data());
  std::vector<std::unique_ptr<TH1D>> hists_p, hists_m;
  // auto hist_Q_1D_rej = std::make_unique<TH1D>("hist_Q_1D_rej", "hist_Q_1D_rej", binEdges.size()-1, binEdges.data());
  std::vector<std::unique_ptr<TH1D>> hist_Q_1D_rej;
  // auto hist_Q_1D_norej = std::make_unique<TH1D>("hist_Q_1D_norej", "hist_Q_1D_norej", binEdges.size()-1, binEdges.data());

  for (int i=0; i<reject_region.size(); i++) {
    hist_Q_1D_rej.push_back(std::make_unique<TH1D>(Form("hist_Q_1D_rej_%d", i), Form("hist_Q_1D_rej_%d", i), binEdges.size()-1, binEdges.data()));
  }
  for (int i=0; i<reject_region.size()+1; i++) {
    hists_p.push_back(std::make_unique<TH1D>(Form("hist_p_rej_%d", i), Form("hist_p_rej_%d", i), binEdges.size()-1, binEdges.data()));
    hists_m.push_back(std::make_unique<TH1D>(Form("hist_m_rej_%d", i), Form("hist_m_rej_%d", i), binEdges.size()-1, binEdges.data()));
    hists_p[i]->SetDirectory(0);
    hists_m[i]->SetDirectory(0);
  }

  hist_Q_1D->SetDirectory(0);
  hist_Q_1D_p->SetDirectory(0);
  hist_Q_1D_diff->SetDirectory(0);
  // hist_Q_1D_rej->SetDirectory(0);
  // hist_Q_1D_norej->SetDirectory(0);
  for (auto& hist : hist_Q_1D_rej) {
    hist->SetDirectory(0);
  }

  std::vector<double> chi_values(Q_values);
  std::vector<double> top_Q(Q_values);
  std::vector<double> top_Q_p(Q_values);
  // std::vector<double> top_Q_rej(Q_values);
  // std::vector<bool> top_Q_rej(Q_values);
  // std::vector<double> top_Q_norej(Q_values);
  // std::vector<bool> top_Q_norej(Q_values);
  std::vector<std::vector<double>> top_Q_rej(reject_region.size(), std::vector<double>(Q_values));
  for (int i = 0; i < Q_values; ++i) {
    chi_values[i] = 0;
    top_Q[i] = 0;
    top_Q_p[i] = 0;
    // top_Q_rej[i] = 0;
    // top_Q_norej[i] = 0;
    for (size_t j = 0; j < reject_region.size(); ++j) {
      top_Q_rej[j][i] = 0;
    }
  }
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
  //define rejection region as Qout, Qside, Qlong
  double rejFromOsl1d[3] = {0., 0.38, 0.};
  double rejToOsl1d[3] = {0.46, 0.62, 0.34};
  double norejFromOsl1d[3] = {0.18, 0.02, 0.02};
  double norejToOsl1d[3] = {0.54, 0.14, 0.06};
  double min_chi = 0;
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

      // min_chi = *std::min_element(chi_values.begin(), chi_values.end());
      // int min_index = std::distance(chi_values.begin(), std::find(chi_values.begin(), chi_values.end(), min_chi));
      
      chi = sqrt(pow(value - expected, 2));
      // if ((chi > 0.2)&&(sqrt(pow(coords[0], 2) + pow(coords[1], 2) + pow(coords[2], 2) < 0.6)))
      if ((chi > 0.5))
        std::cout<<"value "<<value<<" expected "<<expected<<" chi "<<chi<< " coords "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
      // if (chi > min_chi) {
      //     chi_values[min_index] = chi;
      //     top_Q[min_index] = 0;
      //     for (int j = 0; j < 3; ++j) {
      //         top_Q[min_index] += pow(coords[j], 2);
      //     }
      //     top_Q[min_index] = sqrt(top_Q[min_index]);
      //     top_Q_rej[min_index] = true;
      //     top_Q_norej[min_index] = true;
      //     //check if coords are in rejection region
      //     for (int j = 0; j < 3; ++j) {
      //         if (!(coords[j] >= rejFromOsl1d[j] && coords[j] <= rejToOsl1d[j])) {
      //             top_Q_rej[min_index] = false;
      //         }
      //         if (!(coords[j] >= norejFromOsl1d[j] && coords[j] <= norejToOsl1d[j])) {
      //             top_Q_norej[min_index] = false;
      //         }
      //         if (!top_Q_rej[min_index] && !top_Q_norej[min_index]) {
      //             break;
      //         }
      //     }
      // }
      // cout<<"chi "<<chi<<" value "<<value<<" expected "<<expected<<" coords "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<endl;
      if ((value - expected) > 0) {
        DataStruc_chi_p.Add(coords, chi, error);
        DataStruc_chi_m.Add(coords, 0);
      } else if ((value - expected) < 0) {
        DataStruc_chi_m.Add(coords, chi, error);
        DataStruc_chi_p.Add(coords, 0);
      }
      // Add the chi-squared value to the new BinData object
      avgChi += chi;
      if (chi > maxChi) {
        maxChi = chi;
      }
  }
  // fill hist_Q_1D with chi values
  // for (int i = 0; i < Q_values; ++i) {
  //     hist_Q_1D->Fill(top_Q[i]);
  // }
  // // Draw the histogram
  // // hist_Q_1D->Draw("BOX2Z");
  // hist_Q_1D->Write(("hist_Q_1D_" + plotName).c_str());
  // cout<<"peak_value_all "<<peak_value_all<<" peak_fit_all "<<peak_fit_all<<" bins_peak_all "<<bins_peak_all<<" ratio "<<(peak_value_all-peak_fit_all)/bins_peak_all<<endl;
  avgChi /= DataStruc.Size();
  cout<<"avgChi "<<avgChi<<" maxChi "<<maxChi<<endl;

  // c1->Write(("hist_Q_1D_" + plotName).c_str());
  // scale average chi
  // avgChi *= 3;
  avgChi *= 1;
  std::vector<std::unique_ptr<TH3D>> Hist_p;
  std::vector<std::unique_ptr<TH3D>> Hist_m;
  std::vector<std::unique_ptr<TH3D>> Hist_p_rej;
  std::vector<std::unique_ptr<TH3D>> Hist_m_rej;

  Hist_p.reserve(InHistRatio.size());
  Hist_m.reserve(InHistRatio.size());
  Hist_p_rej.reserve(InHistRatio.size());
  Hist_m_rej.reserve(InHistRatio.size());

  for (const auto& histtemp : InHistRatio) {
    if (!histtemp) {
        std::cerr << "Error: histtemp is a null pointer\n";
        continue;
    }

    auto newHist = std::make_unique<TH3D>(*histtemp);
    newHist->SetDirectory(0);
    newHist->Reset();
    Hist_p.push_back(std::move(newHist));

    auto newHist2 = std::make_unique<TH3D>(*histtemp);
    newHist2->SetDirectory(0);
    newHist2->Reset();
    Hist_m.push_back(std::move(newHist2));

    auto newHistrej = std::make_unique<TH3D>(*histtemp);
    newHistrej->SetDirectory(0);
    newHistrej->Reset();
    Hist_p_rej.push_back(std::move(newHistrej));

    auto newHist2rej = std::make_unique<TH3D>(*histtemp);
    newHist2rej->SetDirectory(0);
    newHist2rej->Reset();
    Hist_m_rej.push_back(std::move(newHist2rej));
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
  // std::vector<bool> reject_region = rej_OSL();
  // for(int i=4-hist_number;i<4;i++) {
  for(int i=0;i<hist_number;i++) {
    // int ixStart = InHistRatio[i]->GetXaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    // int iyStart = InHistRatio[i]->GetYaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    // int izStart = InHistRatio[i]->GetZaxis()->FindBin(qMin3D[i-(4-hist_number)]);
    int ixStart = InHistRatio[i]->GetXaxis()->FindBin(qMin3D[i]);
    int iyStart = InHistRatio[i]->GetYaxis()->FindBin(qMin3D[i]);
    int izStart = InHistRatio[i]->GetZaxis()->FindBin(qMin3D[i]);
    cout<<"ix starting from "<<ixStart<<" total bins "<<InHistRatio[i]->GetNbinsX()<<endl;
    // cout<<"ix starting from "<<InHistRatio[4-hist_number]->GetXaxis()->FindBin(qMin3D[0])<<" total bins "<<InHistRatio[i]->GetNbinsX()<<endl;
    // for (int ix=InHistRatio[4-hist_number]->GetXaxis()->FindBin(qMin3D[0]); ix<=InHistRatio[i]->GetNbinsX(); ++ix) {
    //   for (int iy=InHistRatio[4-hist_number]->GetYaxis()->FindBin(qMin3D[0]); iy<=InHistRatio[i]->GetNbinsY(); ++iy) {
    //     for (int iz=InHistRatio[4-hist_number]->GetZaxis()->FindBin(qMin3D[0]); iz<=InHistRatio[i]->GetNbinsZ(); ++iz) {
    for (int ix=InHistRatio[0]->GetXaxis()->FindBin(qMin3D[0]); ix<=InHistRatio[i]->GetNbinsX(); ++ix) {
      for (int iy=InHistRatio[0]->GetYaxis()->FindBin(qMin3D[0]); iy<=InHistRatio[i]->GetNbinsY(); ++iy) {
        for (int iz=InHistRatio[0]->GetZaxis()->FindBin(qMin3D[0]); iz<=InHistRatio[i]->GetNbinsZ(); ++iz) {
          double binCoord[3] ={InHistRatio[i]->GetXaxis()->GetBinCenter(ix), 
                               InHistRatio[i]->GetYaxis()->GetBinCenter(iy),
                               InHistRatio[i]->GetZaxis()->GetBinCenter(iz)};
          if (((ix<ixStart)&&(iy<iyStart)&&(iz<izStart))
            ||(binCoord[0] > qMax)||(binCoord[1] > qMax)||(binCoord[2] > qMax)
            ||(binCoord[0] < qMin)||(binCoord[1] < qMin)||(binCoord[2] < qMin)) {
            // InHistRatio[i]->SetBinContent(ix, iy, iz, 0);
            // InHistRatio[i]->SetBinError(ix, iy, iz, 0);
            Hist_p[i]->SetBinContent(ix, iy, iz, 0);
            Hist_p[i]->SetBinError(ix, iy, iz, 0);
            Hist_m[i]->SetBinContent(ix, iy, iz, 0);
            Hist_m[i]->SetBinError(ix, iy, iz, 0);
            continue;
          }
          min_chi = *std::min_element(chi_values.begin(), chi_values.end());
          int min_index = std::distance(chi_values.begin(), std::find(chi_values.begin(), chi_values.end(), min_chi));
          int bin = coordinates[std::make_tuple(binCoord[0], binCoord[1], binCoord[2])];
          double value_p, value_m, inv_error_p, inv_error_m;
          // get bin number for coordinates
          // int bin_test = coordinates[std::make_tuple(-2, -2, -2)];
          // const double* coords_test = DataStruc_chi.GetPoint(bin_test, value, inv_error);
          // cout<<"value "<<value<<"bin "<<bin_test<<"<- coords "<<coords_test[0]<<"  "<<coords_test[1]<<"  "<<coords_test[2]<<endl;
          const double* coords_p = DataStruc_chi_p.GetPoint(bin, value_p, inv_error_p);
          const double* coords_m = DataStruc_chi_m.GetPoint(bin, value_m, inv_error_m);
          // double binContent = inHistRatio->GetBinContent(ix, iy, iz);
          // double binError = inHistRatio->GetBinError(ix, iy, iz); 
          // if (value_p*pow(qBins3D_width[i]/qBins3D_width[4-hist_number],3)*inv_error_p > min_chi) {
          if (i<hist_number-1) {
            if ((value_p*inv_error_p > min_chi)&&(value_p*inv_error_p > 0.001)) {
              // chi_values[min_index] = value_p*pow(qBins3D_width[i]/qBins3D_width[4-hist_number],3)*inv_error_p;
              chi_values[min_index] = value_p*inv_error_p;
              // cout<<binCoord[0]<<"  "<<binCoord[1]<<"  "<<binCoord[2]<<"\t value "<<value_m<<"  error "<<1.0/inv_error_p<<"  chi "<<chi_values[min_index];
              top_Q_p[min_index] = 0;
              for (int j = 0; j < 3; ++j) {
                  top_Q_p[min_index] += pow(coords_p[j], 2);
              }
              top_Q_p[min_index] = sqrt(top_Q_p[min_index]);
              // cout<<"  top_Q_p "<<top_Q_p[min_index]<<endl;
            }
            // if (value_m*pow(qBins3D_width[i]/qBins3D_width[4-hist_number],3)*inv_error_m > min_chi) {
            if ((value_m*inv_error_m > min_chi)&&(value_m*inv_error_m > 0.001)) {
              // chi_values[min_index] = value_m*pow(qBins3D_width[i]/qBins3D_width[4-hist_number],3)*inv_error_m;
              chi_values[min_index] = value_m*inv_error_m;
              // cout<<binCoord[0]<<"  "<<binCoord[1]<<"  "<<binCoord[2]<<"\t value "<<value_m<<"  error "<<1.0/inv_error_m<<"  chi "<<-chi_values[min_index];
              top_Q[min_index] = 0;
              for (int j = 0; j < 3; ++j) {
                  top_Q[min_index] += pow(coords_m[j], 2);
              }
              top_Q[min_index] = sqrt(top_Q[min_index]);
              // cout<<"  top_Q "<<top_Q[min_index]<<endl;
              // top_Q_rej[min_index] = 1;
              // top_Q_norej[min_index] = 1;
              // //check if coords are not in rejection region
              // for (int j = 0; j < 3; ++j) {
              //   if (reject_region[0])
              //     if (!(coords_m[j] >= rejFromOsl[0][j] && coords_m[j] <= rejToOsl[0][j])) {
              //       top_Q_rej[min_index] = 0;
              //     }
              //   if (reject_region[1])
              //     if (!(coords_m[j] >= rejFromOsl[1][j] && coords_m[j] <= rejToOsl[1][j])) {
              //       top_Q_norej[min_index] = 0;
              //     }
              //   if (!top_Q_rej[min_index] && !top_Q_norej[min_index]) {
              //     break;
              //   }
              // }
              
              // for (int j = 0; j < 3; ++j) {
              //     if (!(coords_m[j] >= rejFromOsl1d[j] && coords_m[j] <= rejToOsl1d[j])) {
              //         top_Q_rej[min_index] = 0;
              //     }
              //     if (!(coords_m[j] >= norejFromOsl1d[j] && coords_m[j] <= norejToOsl1d[j])) {
              //         top_Q_norej[min_index] = 0;
              //     }
              //     if (!top_Q_rej[min_index] && !top_Q_norej[min_index]) {
              //         break;
              //     }
              // }
            }
            // Initialize top_Q_rej for all rejection regions
            if (((value_m*inv_error_m > min_chi)&&(value_m*inv_error_m > 0.001))||((value_p*inv_error_p > min_chi)&&(value_p*inv_error_p > 0.001))) {
              for (size_t k = 0; k < reject_region.size(); ++k) {
                top_Q_rej[k][min_index] = std::max(top_Q[min_index], top_Q_p[min_index]);
              }
              // cout<<"min_chi value: "<<min_chi<<" min index "<<min_index<<" top_Q_p "<<top_Q_p[min_index]<<" top_Q "<<top_Q[min_index]<<endl;
              // Check if coords are not in rejection regions
              for (int j = 0; j < 3; ++j) {
                for (size_t k = 0; k < reject_region.size(); ++k) {
                  if (reject_region[k]) {
                    if (!(coords_m[j] >= rejFromOsl[k][j] && coords_m[j] <= rejToOsl[k][j])) {
                      top_Q_rej[k][min_index] = 0;
                    }
                  }
                }
              }
            }
          }
          bool reject = false;
          for (int k = 0; k < reject_region.size(); k++) {
            if (reject_region[k]) {
              if ((binCoord[0] >= rejFromOsl[k][0]) && (binCoord[0] <= rejToOsl[k][0]) &&
                  (binCoord[1] >= rejFromOsl[k][1]) && (binCoord[1] <= rejToOsl[k][1]) &&
                  (binCoord[2] >= rejFromOsl[k][2]) && (binCoord[2] <= rejToOsl[k][2])) {
                reject = true;
                Hist_p[i]->SetBinContent(ix, iy, iz, 0);
                Hist_p[i]->SetBinError(ix, iy, iz, 0);
                Hist_m[i]->SetBinContent(ix, iy, iz, 0);
                Hist_m[i]->SetBinError(ix, iy, iz, 0);
                Hist_p_rej[i]->SetBinContent(ix, iy, iz, value_p);
                Hist_p_rej[i]->SetBinError(ix, iy, iz, 0);
                Hist_m_rej[i]->SetBinContent(ix, iy, iz, value_m);
                Hist_m_rej[i]->SetBinError(ix, iy, iz, 0);
                break;
              }
            }
          }
          if (reject) continue;
          // if ((binCoord[0] >= rejFromOsl[0]) && (binCoord[0] <= rejToOsl[0]) &&
          //       (binCoord[1] >= rejFromOsl[1]) && (binCoord[1] <= rejToOsl[1]) &&
          //       (binCoord[2] >= rejFromOsl[2]) && (binCoord[2] <= rejToOsl[2])) {
          //   Hist_p[i]->SetBinContent(ix, iy, iz, 0);
          //   Hist_p[i]->SetBinError(ix, iy, iz, 0);
          //   Hist_m[i]->SetBinContent(ix, iy, iz, 0);
          //   Hist_m[i]->SetBinError(ix, iy, iz, 0);
          //   continue;
          // }

          //error = 1./inv_error;
          //data.Add(binCoord, binContent, binError); //binCoordErr,
          // Hist[i]->SetBinContent(ix, iy, iz, value);
          // if only bins with chi 3x bigger than average chi
          float bin_scale=1;//pow(qWidth3D[i]/qWidth3D[0], 1.5);
          // InHistRatio[i]->SetBinContent(ix, iy, iz, value/bin_scale);



          // if(value_p > (avgChi))
          Hist_p[i]->SetBinContent(ix, iy, iz, value_p/bin_scale);
          // else
          //   Hist_p[i]->SetBinContent(ix, iy, iz, 0);
          // if(value_m > (avgChi))
          Hist_m[i]->SetBinContent(ix, iy, iz, value_m/bin_scale);
          // else
          //   Hist_m[i]->SetBinContent(ix, iy, iz, 0);


          
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
  draw_hist_chi(Hist_p_rej, "p_rej");
  draw_hist_chi(Hist_m_rej, "m_rej");


  //   // fill hist_Q_1D with chi values
  // for (int i = 0; i < Q_values; ++i) {
  //   // if (chi_values[i] > avgChi) {
  //     if (top_Q_p[i] > 0) {
  //       hist_Q_1D_p->Fill(top_Q_p[i], chi_values[i]);
  //       cout<<"top_Q_p "<<top_Q_p[i]<<" chi "<<chi_values[i]<<endl;
  //     }
  //     if (top_Q[i] > 0) {
  //       hist_Q_1D->Fill(top_Q[i], chi_values[i]);
  //       cout<<"top_Q "<<top_Q[i]<<" chi "<<chi_values[i]<<endl;
  //     }
  //     // hist_Q_1D->Fill(top_Q[i], chi_values[i]);
  //     // if (top_Q_rej[i]) {
  //     //   hist_Q_1D_rej->Fill(top_Q[i], chi_values[i]);
  //     // }
  //     // if (top_Q_norej[i]) {
  //     //   hist_Q_1D_norej->Fill(top_Q[i], chi_values[i]);
  //     // }
  //     for (size_t k = 0; k < reject_region.size(); ++k) {
  //       if (top_Q_rej[k][i]) {
  //         hist_Q_1D_rej[k]->Fill(top_Q[i], chi_values[i]);
  //       }
  //     }
  //   // }
  // }
  // std::vector<std::unique_ptr<TH1D>> hists;
  // // if (hist_Q_1D) {
  // //   std::unique_ptr<TH1D> histtemp_1D;
  // //   auto histtemp = std::make_unique<TH1D>(*hist_Q_1D);
  // //   histtemp_1D = std::make_unique<TH1D>(*hist_Q_1D);
  // //   histtemp_1D->SetDirectory(0);
  // //   hists.push_back(std::move(histtemp_1D));
  // //   histtemp_1D = std::make_unique<TH1D>(*hist_Q_1D_p);
  // //   histtemp_1D->SetDirectory(0);
  // //   hists.push_back(std::move(histtemp_1D));
  // //   histtemp_1D = std::make_unique<TH1D>(*hist_Q_1D_rej);
  // //   histtemp_1D->SetDirectory(0);
  // //   hists.push_back(std::move(histtemp_1D));
  // //   histtemp_1D = std::make_unique<TH1D>(*hist_Q_1D_norej);
  // //   histtemp_1D->SetDirectory(0);
  // //   hists.push_back(std::move(histtemp_1D));
  // //   // hists.push_back(std::make_unique<TH1D>());
  // //   // hists.push_back(std::make_unique<TH1D>());
  // //   // hists.push_back(std::make_unique<TH1D>());
  // //   // hists.push_back(std::make_unique<TH1D>());
  // // }
  // hists.push_back(std::move(hist_Q_1D_p));
  // hists.push_back(std::move(hist_Q_1D));
  // // hists.push_back(std::move(hist_Q_1D_rej));
  // // hists.push_back(std::move(hist_Q_1D_norej));
  // for (auto& hist_rej : hist_Q_1D_rej) {
  //   hists.push_back(std::move(hist_rej));
  // }
  // auto histtemp = std::make_unique<TH1D>(*hists[0]);
  // histtemp->SetDirectory(0);
  // hists.push_back(std::move(histtemp));
  // // hists.push_back(std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(hists[0]->Clone())));
  // // std::vector<TH1D*> hists = {hist_Q_1D.get(), hist_Q_1D_p.get(), hist_Q_1D_rej.get(), hist_Q_1D_norej.get()};
  // std::vector<std::vector<double>> top_Qs = {top_Q_p, top_Q};//, top_Q_rej, top_Q_norej};
  // for (const auto& top_Q_rej_vec : top_Q_rej) {
  //   top_Qs.push_back(top_Q_rej_vec);
  // }

  // // std::vector<std::vector<std::pair<double, int>>> binData(hists.size());
  // std::vector<std::vector<double>> binavg(hists.size(), std::vector<double>(binEdges.size()-1, 0.0));
  // std::vector<std::vector<int>> binnum(hists.size(), std::vector<int>(binEdges.size()-1, 0));
  // // cout<<"creating binData of size"<<hists.size()<<endl;

  // // for (int j = 0; j < hists.size(); ++j) {
  // //   // binData[j].resize(hists[j]->GetNbinsX(), {0, 0});
  // // }

  // for (int i = 0; i < Q_values; ++i) {
  //   // if (chi_values[i] > avgChi) {
  //     for (int j = 0; j < hists.size()-1; ++j) {
  //       if (top_Qs[j][i] > 0) {
  //         int bin = hists[j]->FindBin(top_Qs[1][i]);
  //         if (j == 0) 
  //           bin = hists[j]->FindBin(top_Qs[0][i]);
  //         // cout<<"bin "<<bin<<" chi "<<chi_values[i]<<"in hist "<<j<<" top_Q "<<top_Qs[0][i]<<endl;
  //         // binData[j][bin-1].first += chi_values[i];
  //         binavg[j][bin-1] += chi_values[i];
  //         // binData[j][bin-1].second++;
  //         binnum[j][bin-1]++;
  //         // if (j == 1)
  //         //   cout<<"bin "<<bin<<" chi "<<chi_values[i]<<"in hist "<<j<<" top_Q_p "<<top_Qs[1][i]<<endl;
  //         // else
  //         //   cout<<"bin "<<bin<<" chi "<<chi_values[i]<<"in hist "<<j<<" top_Q "<<top_Qs[0][i]<<endl;
  //         // cout<<"bindata size"<<binData[j].size()<<endl;
  //       }
  //       // cout<<"ALWAYS bin "<<" chi "<<chi_values[i]<<"in hist "<<j<<" top_Q "<<top_Qs[0][i]<<endl;
  //     }
  //   // }
  // }
  // // cout<<"filling hists"<<endl;

  // // Fill the histograms with the average value for each bin
  // for (int j = 0; j < hists.size()-1; ++j) {
  //   for (int bin = 0; bin < hists[j]->GetNbinsX(); bin++) {
  //     // if (binData[j][bin].second > 0) {  // Check if the count is greater than 0
  //     if (binnum[j][bin] > 0) {  // Check if the count is greater than 0
  //       // double average = binData[j][bin].first / binData[j][bin].second;
  //       double average = binavg[j][bin] / binnum[1][bin];
  //       if (j == 0) {
  //         average = binavg[j][bin] / binnum[0][bin];
  //       }
  //       hists[j]->SetBinContent(bin+1, average);
  //       cout<<"bin "<<bin<<" average "<<average<<"in hist "<<j<<endl;
  //     }
  //   }
  // }

  // Define and initialize histograms
  std::vector<std::unique_ptr<TH1D>> hists;
  hists.push_back(std::move(hist_Q_1D_p));
  hists.push_back(std::move(hist_Q_1D));

  for (auto& hist_rej : hist_Q_1D_rej) {
      hists.push_back(std::move(hist_rej));
  }
  // auto histtemp = std::make_unique<TH1D>(*hists[0]);
  // histtemp->SetDirectory(0);
  // hists.push_back(std::move(histtemp));
  hists.push_back(std::move(hist_Q_1D_diff));

  // Define top_Qs as a vector of different top_Q vectors
  std::vector<std::vector<double>> top_Qs = {top_Q_p, top_Q};
  for (const auto& top_Q_rej_vec : top_Q_rej) {
      top_Qs.push_back(top_Q_rej_vec);
  }

  // Initialize vectors to store bin averages and counts
  std::vector<std::vector<double>> binavg(hists.size(), std::vector<double>(binEdges.size(), 0.0));
  std::vector<std::vector<int>> binnum(hists.size(), std::vector<int>(binEdges.size(), 0));
  std::vector<std::vector<double>> binavg_rej_p(reject_region.size()+1, std::vector<double>(binEdges.size(), 0.0));
  std::vector<std::vector<double>> binavg_rej_m(reject_region.size()+1, std::vector<double>(binEdges.size(), 0.0));
  // Fill binavg and binnum arrays
  for (int i = 0; i < Q_values; ++i) {
      for (int j = 0; j < top_Qs.size(); ++j) {
          if (top_Qs[j][i] > 0) {
              int bin = hists[j]->FindBin(top_Qs[j][i]);
              if (bin > binEdges.size()) {
                  continue;
              }
              binavg[j][bin - 1] += chi_values[i];
              binnum[j][bin - 1]++;
              if (j < 2) {
                binavg[hists.size()-1][bin - 1] += chi_values[i]*(1-j*2);
                // cout<<"binavg "<<binavg[hists.size()-1][bin - 1]<<" chi "<<chi_values[i]<<" j "<<j<<endl;
                binnum[hists.size()-1][bin - 1]++;
              }
              if (j > 1) {
                if (top_Qs[0][i] > 0) {
                  binavg_rej_p[j-2][bin - 1] += chi_values[i];
                }
                if (top_Qs[1][i] > 0) {
                  binavg_rej_m[j-2][bin - 1] += chi_values[i];
                }
              }
          }
      }
  }

  // Fill the histograms with the average value for each bin
  for (int j = 0; j < hists.size()-1; ++j) {
      for (int bin = 0; bin < hists[j]->GetNbinsX(); ++bin) {
          if (binnum[j][bin] > 0) {
              double average = binavg[j][bin] / binnum[j][bin];
              hists[j]->SetBinContent(bin + 1, average);
              // cout << "bin " << bin << "binall "<<binavg[j][bin]<<" binnum "<<binnum[j][bin]<< " average " << average << " in hist " << j << endl;
          }
      }
  }
  // std::vector<std::unique_ptr<TH1D>> hists_p, hists_m;
  // for (int i=0; i<reject_region.size()+1; i++) {
  //   hists_p.push_back(std::make_unique<TH1D>(*hists[0]));
  //   hists_m.push_back(std::make_unique<TH1D>(*hists[1]));
  //   hists_p[i]->SetDirectory(0);
  //   hists_m[i]->SetDirectory(0);
  // }
  //copy contents of hists[0] and hists[1] to existing hists_p and hists_m
  hists_p[0]->Add(hists[0].get());
  hists_m[0]->Add(hists[1].get());
  for (int j = 0; j < reject_region.size(); ++j) {
    for (int bin = 0; bin < hists_p[j]->GetNbinsX(); ++bin) {
      if ((binnum[0][bin] > 0)&&(binnum[1][bin] > 0)) {
        double average_p = binavg_rej_p[j][bin] / binnum[0][bin];
        double average_m = binavg_rej_m[j][bin] / binnum[1][bin];
        hists_p[j+1]->SetBinContent(bin + 1, average_p);
        hists_m[j+1]->SetBinContent(bin + 1, average_m);
        // cout << "bin " << bin << "binall "<<binavg_rej_p[j][bin]<<" binnum "<<binnum[0][bin]<< " average " << average_p << " in hist " << j << endl;
      }
    }
  }
  //fill the last histogram EVERYWHERE HERE SHOULD BE >GetNbinsX() +1, 
  //however does not matter as its the end of the histogram
  for (int bin = 0; bin < hists[hists.size()-1]->GetNbinsX(); ++bin) {
    if (binnum[hists.size()-1][bin] > 0) {
      double average = binavg[hists.size()-1][bin] / binnum[hists.size()-1][bin];
      hists[hists.size()-1]->SetBinContent(bin + 1, average);
      // cout << "bin " << bin << "binall "<<binavg[hists.size()-1][bin]<<" binnum "<<binnum[hists.size()-1][bin]<< " average " << average << " in hist " << hists.size()-1 << endl;
    } 
    // else {
    //   hists[hists.size()-1]->SetBinContent(bin + 1, 0);
      // cout << "bin " << bin << "binall "<<binavg[hists.size()-1][bin]<<" binnum "<<binnum[hists.size()-1][bin]<< " average " << 0 << " in hist " << hists.size()-1 << endl;
    // }
  }
  // fill in seperate histograms for 
  // hists[4]->Add(hists[0].get(), hists[1].get(), 1, -1);
  // hists.back()->Add(hists[0].get(), hists[1].get(), 1, -1);

  // Draw the histogram
  // hist_Q_1D->Draw("BOX2Z");
  // Create a canvas
  // auto c1 = std::make_unique<TCanvas>("c1", "hist_Q_1D", 800, 600);
  // // auto hs = std::make_unique<THStack>("hs","hist_Q_1D");
  // // gStyle->SetPalette(kOcean);

  // // Draw the histograms on the canvas
  // hist_Q_1D->SetFillColor(kBlue);
  // hist_Q_1D->SetFillStyle(3005);
  // hist_Q_1D_p->SetFillColor(kCyan);
  // hist_Q_1D_p->SetFillStyle(3007);
  // hist_Q_1D_rej->SetFillColor(kRed);
  // hist_Q_1D_rej->SetFillStyle(3021);
  // hist_Q_1D_norej->SetFillColor(kGreen);
  // hist_Q_1D_norej->SetFillStyle(3009);
  // // hs->Add(hist_Q_1D.get());
  // // hs->Add(hist_Q_1D_rej.get());
  // // hs->Draw("pfc nostack");
  // // c1->Update();
  // hist_Q_1D->Draw("hist");
  // hist_Q_1D_p->Draw("hist same");
  // hist_Q_1D_rej->Draw("hist same");
  // hist_Q_1D_norej->Draw("hist same");

  int lastNonZeroBin = hists[0]->GetNbinsX();
  for (int bin = hists[0]->GetNbinsX(); bin >= 1; --bin) {
      if (hists[0]->GetBinContent(bin) > 0) {
          lastNonZeroBin = bin;
          break;
      }
  }
  double xMin = hists[0]->GetXaxis()->GetXmin();
  double xMax = hists[0]->GetXaxis()->GetBinLowEdge(lastNonZeroBin) + hists[0]->GetXaxis()->GetBinWidth(lastNonZeroBin);
  for (auto& hist : hists) {
      hist->GetXaxis()->SetRangeUser(xMin, xMax);
  }
  for (auto& hist : hists_p) {
      hist->GetXaxis()->SetRangeUser(xMin, xMax);
  }
  for (auto& hist : hists_m) {
      hist->GetXaxis()->SetRangeUser(xMin, xMax);
  }
  std::vector<int> colors = {600, 632, 432, 416, 880, 921, 616, 800};//{kBlue, kRed, kCyan, kGreen, kViolet, kPink, kMagenta, kYellow}; // in numbers: {600, 632, 432, 416, 880, 900, 616, 800}
  std::vector<int> styles = {3545, 3554, 3021, 3009, 3013, 3017, 3025, 3033};
  std::vector<int> markers = {20, 21, 22, 23, 24, 25, 26, 27};
  // Ensure the colors vector is large enough
  if (colors.size() < hists_p.size() + 1) {
    std::cerr << "Error: The colors vector is not large enough to accommodate all histograms." << hists_p.size()+1 << std::endl;
    return;
  }
  // // cut colors and styles to the number of histograms
  colors.resize(hists.size()-1);
  styles.resize(hists.size()-1);
  //add for last histogram color black and no infill
  colors.push_back(kBlack);
  styles.push_back(0);
  double hist_max = 0;
  // Apply styles and colors
  for (int j = 0; j < hists.size(); ++j) {
    hists[j]->SetStats(0);
    //set marker color
    if (j > 1) {
      hists[j]->SetMarkerColor(colors[j]);
      hists[j]->SetMarkerStyle(markers[j-2]);
    } else {
      hists[j]->SetFillColor(colors[j]);
      hists[j]->SetFillStyle(styles[j]);
      hists[j]->SetLineColor(colors[j]);
    }
    //get the maximal value from all histograms
    if (hists[j]->GetMaximum() > hist_max) {
      hist_max = hists[j]->GetMaximum();
    }
    std::cout<<"hist "<<j<<" color "<<colors[j]<<" style "<<styles[j]<<std::endl;
  }
  for (int j = 0; j < hists.size()-1; ++j) {
    hists[j]->SetMaximum(hist_max*1.1);
  }
  vector<string> ref_region;
  for (int i = 0; i < reject_region.size(); ++i) {
    if (i < 3) {
      ref_region.push_back("ULS");
    } else {
      ref_region.push_back("OHP");
    }
  }
  // Create a canvas
  auto c1 = std::make_unique<TCanvas>("c1", "hist_Q_1D_chi", 800, 600);
  // TCanvas* c1 = new TCanvas("c1", "hist_Q_1D_chi", 800, 600);
  c1->cd();
  auto pad1 = std::make_unique<TPad>("pad1", "pad1", 0, 0.25, 1, 1);
  auto pad2 = std::make_unique<TPad>("pad2", "pad2", 0, 0.01, 1, 0.25);
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  hists[0]->SetTitle("");
  hists[0]->GetYaxis()->SetTitle(("|"+ titleR2 + " - fit| regions").c_str());

  // Draw the first two histograms as histograms
  for (int j = 0; j < 2; ++j) {
      std::string option = (j == 0) ? "HIST" : "HIST SAME";
      hists[j]->Draw(option.c_str());
  }
  cout<<"drawn rej p,m"<<endl;
  for (int j = 2; j < hists.size()-1; ++j) {
    // // std::string option = (j == 0) ? "hist" : "hist same";
    // std::string option = (j == 0) ? "P" : "P SAME";
    // // std::string option = "hist same";
    // hists[j]->Draw(option.c_str());
    hists[j]->Draw("P SAME");
  }
  cout<<"drawn rej regions"<<endl;
  // Ensure that the histogram and canvas are correctly allocated
  if (!hists.back()) {
    std::cerr << "Error: Last histogram is a null pointer" << std::endl;
    return;
  }

  if (!c1) {
    std::cerr << "Error: Canvas c1 is a null pointer" << std::endl;
    return;
  }
  //draw legend
  auto leg = std::make_unique<TLegend>(0.5, 0.6, 0.9, 0.9);
  leg->AddEntry(hists[0].get(), "Difference above fit", "f");
  leg->AddEntry(hists[1].get(), "Difference below fit", "f");
  for (int k = 0; k < reject_region.size(); ++k) {
    string region_range = "Rej region " + std::to_string(k+1) + " " + ref_region[k] + ": ";
    for (int l = 0; l < 3; ++l) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(2) 
             << rejFromOsl[k][l] << "-" << rejToOsl[k][l];
      region_range += stream.str() + "; ";
    }
    leg->AddEntry(hists[k+2].get(), region_range.c_str(), "p");
  }
  leg->Draw();
  pad2->cd();
  //dont draw stats for last histogram and title
  hists.back()->SetTitle("");
  hists.back()->SetStats(0);
  hists.back()->GetXaxis()->SetLabelSize(0.15);
  hists.back()->GetXaxis()->SetTitle("Q[GeV]");
  hists.back()->GetXaxis()->SetTitleSize(0.15);
  hists.back()->GetYaxis()->SetLabelSize(0.15);
  hists.back()->GetYaxis()->SetTitle((titleR2 + " - fit").c_str());
  hists.back()->GetYaxis()->SetTitleSize(0.15);
  hists.back()->GetYaxis()->SetTitleOffset(0.2);
  cout<<"drawning last histogram"<<endl;
  hists.back()->Draw("hist");
  //get contents of last histogram
  // for (int i = 1; i < hists.back()->GetNbinsX()+1; ++i) {
  //   cout<<"bin center "<<hists.back()->GetBinCenter(i)<<" content "<<hists.back()->GetBinContent(i)<<endl;
  // }
  std::unique_ptr<TLine> line = std::make_unique<TLine>(xMin, 0, xMax, 0);
  line->SetLineStyle(2);
  line->Draw();
  cout<<"drawn line"<<endl;
  c1->Write(("hist_Q_1D_" + plotName).c_str());
  // hists.clear();
  // Check if the histogram can be drawn
  std::cout << "Histograms and canvas drawn successfully." << std::endl;

  //draw histograms for rejection regions p and m
  double hist_max_p = 0;
  double hist_max_m = 0;
  for (int j = 0; j < hists_p.size(); ++j) {
    hists_p[j]->SetStats(0);
    hists_m[j]->SetStats(0);
    //set marker color
    if (j < 1) {
      hists_p[0]->SetFillColor(colors[0]);
      hists_p[0]->SetFillStyle(styles[0]);
      hists_p[0]->SetLineColor(colors[0]);
      hists_m[0]->SetFillColor(colors[1]);
      hists_m[0]->SetFillStyle(styles[1]);
      hists_m[0]->SetLineColor(colors[1]);
      cout<<" hist color "<<hists_p[j]->GetFillColor()<<endl;
    } else {
      // string title = string(hists_p[j]->GetTitle()) + " rej " + std::to_string(j);
      // hists_p[j]->SetTitle(title.c_str());
      // hists_m[j]->SetTitle(title.c_str());
      // hists_p[j]->SetName(title.c_str());
      // hists_m[j]->SetName(title.c_str());
      hists_p[j]->SetFillColor(colors[j+1]);
      hists_p[j]->SetFillStyle(styles[j+1]);
      hists_p[j]->SetLineColor(colors[j+1]);
      hists_m[j]->SetFillColor(colors[j+1]);
      hists_m[j]->SetFillStyle(styles[j+1]);
      hists_m[j]->SetLineColor(colors[j+1]);
      // hists_p[j]->SetMarkerColor(colors[j+1]);
      // hists_p[j]->SetMarkerStyle(markers[j-1]);
      // hists_m[j]->SetMarkerColor(colors[j+1]);
      // hists_m[j]->SetMarkerStyle(markers[j-1]);
      cout<<" hist "<<j<<" color "<<hists_p[j]->GetFillColor()<<endl;
    }
    //get the maximal value from all histograms
    if (hists_p[j]->GetMaximum() > hist_max_p) {
      hist_max_p = hists_p[j]->GetMaximum();
    }
    if (hists_m[j]->GetMaximum() > hist_max_m) {
      hist_max_m = hists_m[j]->GetMaximum();
    }
    // std::cout<<"hist "<<j<<" color "<<colors[j+1]<<" style "<<styles[j]<<std::endl;
    //get color from the first histogram

  }
  // check colors of histograms
  for (int j = 0; j < hists_p.size(); ++j) {
    cout<<" hist "<<j<<" color "<<hists_p[j]->GetFillColor()<<endl;
  }
  // hists_p[1]->SetFillColor(432);
  // hists_p[2]->SetFillColor(416);
  // hists_p[3]->SetFillColor(880);
  // hists_p[4]->SetFillColor(921);
  // hists_m[1]->SetFillColor(432);
  // hists_m[2]->SetFillColor(416);
  // hists_m[3]->SetFillColor(880);
  // hists_m[4]->SetFillColor(921);
  for (int j = 0; j < hists_p.size(); ++j) {
    hists_p[j]->SetMaximum(hist_max_p*1.1);
    hists_m[j]->SetMaximum(hist_max_m*1.1);
  }
  // Create 2 canvas, we can in loop, we dont need two pads now
  auto c2_p = std::make_unique<TCanvas>("c2_p", "hist_Q_1D_chi_p_rej", 800, 600);
  c2_p->cd();
  hists_p[0]->SetTitle("");
  hists_p[0]->GetYaxis()->SetTitle(("|" + titleR2 + " - fit| regions").c_str());
  hists_p[0]->GetXaxis()->SetTitle("Q[GeV]");
  hists_p[0]->Draw("hist");
  auto hs = std::make_unique<THStack>("hs", "stacked_rej");
  if (reject_region.size() > 0) {
    for (int j = 1; j < hists_p.size(); ++j) {
      hs->Add(hists_p[j].get());
    }
    hs->Draw("hist same");
  }
  //draw legend with individual histograms and ranges of rejection regions
  auto legend = std::make_unique<TLegend>(0.5, 0.6, 0.9, 0.9);
  legend->AddEntry(hists_p[0].get(), "Difference above fit", "f");
      
  // add to the title the rejection region Qosl numbers
  for (int k = 0; k < reject_region.size(); ++k) {
    string region_range = "Rej region " + std::to_string(k+1) + " " + ref_region[k] + ": ";
    for (int l = 0; l < 3; ++l) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(2) 
             << rejFromOsl[k][l] << "-" << rejToOsl[k][l];
      region_range += stream.str() + "; ";
    }
    legend->AddEntry(hists_p[k+1].get(), region_range.c_str(), "f");
  }
  legend->Draw();
  // c2_p->Update();
  c2_p->Write(("hist_Q_1D_rej_p_" + plotName).c_str());

  auto c2_m = std::make_unique<TCanvas>("c2_m", "hist_Q_1D_chi_m_rej", 800, 600);
  c2_m->cd();
  hists_m[0]->SetTitle("");
  hists_m[0]->GetYaxis()->SetTitle(("|" + titleR2 + " - fit| regions").c_str());
  hists_m[0]->GetXaxis()->SetTitle("Q[GeV]");
  hists_m[0]->Draw("hist");
  auto hs2 = std::make_unique<THStack>("hs2", "stacked_rej2");
  if (reject_region.size() > 0) {
    for (int j = 1; j < hists_m.size(); ++j) {
      hs2->Add(hists_m[j].get());
    }
    hs2->Draw("hist same");
  }
  //draw legend with individual histograms and ranges of rejection regions
  auto legend2 = std::make_unique<TLegend>(0.5, 0.6, 0.9, 0.9);
  legend2->AddEntry(hists_m[0].get(), "Difference below fit", "f");
  // add to the title the rejection region Qosl numbers
  for (int k = 0; k < reject_region.size(); ++k) {
    string region_range = "Rej region " + std::to_string(k+1) + " " + ref_region[k] + ": ";
    for (int l = 0; l < 3; ++l) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(2) 
             << rejFromOsl[k][l] << "-" << rejToOsl[k][l];
      region_range += stream.str() + "; ";
    }
    legend2->AddEntry(hists_m[k+1].get(), region_range.c_str(), "f");
  }
  legend2->Draw();
  // c2_m->Update();
  c2_m->Write(("hist_Q_1D_rej_m_" + plotName).c_str());
  cout<<"saved 1d rejection histograms"<<endl;
  
  // delete c1;
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
  size_t pos = plotName.find("c2");
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
  cout<<"writing to txt file"<<endl;
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
  
  cout << "INFO: Saving output files..." << endl;
  outFile->Write();
  outFile->Close();
  // delete outFile;

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
  tableOutput << " & ";
  tableOutput << chi_firstBins.first;
  tableOutput << " \\\\";
  tableOutput << endl;
  tableOutput.close();


  /**
   * CSV output
   */
  string outFilePathCsv = outFilePathBase + "/csv";
  system(("mkdir -p " + outFilePathCsv).c_str());
  string outFilePathCsvCorr = outFilePathCsv + "/correlation";
  system(("mkdir -p " + outFilePathCsvCorr).c_str());
  outFilePathCsv += "/";
  outFilePathCsv += plotName;
  outFilePathCsv += ".csv";
  outFilePathCsvCorr += "/";
  outFilePathCsvCorr += plotName;
  outFilePathCsvCorr += ".csv";
  size_t posi = outFilePathCsvCorr.find_last_of(".");
  outFilePathCsvCorr.insert(posi, "_corr");

  if (verbose) {
    cout << "INFO: Saving CSV file here:" << endl;
    cout << "      " << outFilePathCsv << endl;
  }
//outputCSV(hist, outFilePathCsv);
  outputCSV(func, outFilePathCsv);
  outputCSV(fitData, outFilePathCsvCorr);
}


void Plot() {
  cout << "INFO: Plotting..." << endl;

  string outFilePathBase = "output";
  string today = Today();
  c2_func += std::to_string(c2Index);
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
    system(("mkdir -p " + outFilePath + "/correlation").Data());
    if (verbose) {
      cout << "INFO: Saving " << formatVec.at(i) << " files here:" << endl;
      cout << "      " << outFilePath << endl;
    }
  }

  // TCanvas* canvas = new TCanvas("canvas", "canvas", 50, 50, 600, 600);
  auto canvas = std::make_unique<TCanvas>("canvas_1D", "canvas_1D", 50, 50, 600, 600);

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
  bool useInit = false;

  /**
   * Q_out
   */
  // TH1D* histX;
  // if (histData.find("_cube_") != std::string::npos) {
  //   TH1D* histXpart1 = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                      (histName + "_out_part1").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histXpart2 = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                      (histName + "_out_part2").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histXpart3 = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                      (histName + "_out_part3").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histXpart4 = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                      (histName + "_out_part4").c_str(),
  //                      1, 27, 1, 27));

  //   histX = new TH1D((histName + "_out").c_str(),
  //                    hist->GetTitle(),
  //                    varQbinsSize, &varQbins[0]);

  //   for (int i = 1; i <= 15; ++i) {
  //     histX->SetBinContent(i, histXpart1->GetBinContent(i) / (9 * 9));
  //     histX->SetBinError(i, histXpart1->GetBinError(i) / (9 * 9));
  //   }

  //   histX->SetBinContent(16, histXpart2->GetBinContent(17) / (3 * 3));
  //   histX->SetBinError(16, histXpart2->GetBinError(17) / (3 * 3));
  //   histX->SetBinContent(17, histXpart2->GetBinContent(20) / (3 * 3));
  //   histX->SetBinError(17, histXpart2->GetBinError(20) / (3 * 3));
  //   histX->SetBinContent(18, histXpart2->GetBinContent(23) / (3 * 3));
  //   histX->SetBinError(18, histXpart2->GetBinError(23) / (3 * 3));
  //   histX->SetBinContent(19, histXpart2->GetBinContent(26) / (3 * 3));
  //   histX->SetBinError(19, histXpart2->GetBinError(26) / (3 * 3));

  //   histX->SetBinContent(20, histXpart3->GetBinContent(32));
  //   histX->SetBinError(20, histXpart3->GetBinError(32));
  //   histX->SetBinContent(21, histXpart3->GetBinContent(41));
  //   histX->SetBinError(21, histXpart3->GetBinError(41));
  //   histX->SetBinContent(22, histXpart3->GetBinContent(50));
  //   histX->SetBinError(22, histXpart3->GetBinError(50));

  //   histX->SetBinContent(23, histXpart4->GetBinContent(68));
  //   histX->SetBinError(23, histXpart4->GetBinError(68));
  //   histX->SetBinContent(24, histXpart4->GetBinContent(95));
  //   histX->SetBinError(24, histXpart4->GetBinError(95));
  //   histX->SetBinContent(25, histXpart4->GetBinContent(122));
  //   histX->SetBinError(25, histXpart4->GetBinError(122));
  //   histX->SetBinContent(26, histXpart4->GetBinContent(149));
  //   histX->SetBinError(26, histXpart4->GetBinError(149));
  //   histX->SetBinContent(27, histXpart4->GetBinContent(176));
  //   histX->SetBinError(27, histXpart4->GetBinError(176));

  //   delete histXpart1;
  //   delete histXpart2;
  //   delete histXpart3;
  //   delete histXpart4;
  // } else {
  //   histX = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                 (histName + "_out").c_str(),
  //                 projMin, projMax,
  //                 projMin, projMax));
  //   histX->Scale(1. / (projRange * projRange));
  //   /*histDataStruc = dynamic_cast<TH1D*>(hist->ProjectionX(
  //                 (histName + "_out").c_str(),
  //                 projMin, projMax,
  //                 projMin, projMax));
  //   histDataStruc->Reset("ICESM");*/
  // }
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

  // histX->Print();
  //copy hist into histDataStruc
  // TH1D* histDataStruc = dynamic_cast<TH1D*>(histX->Clone());
  // // TH1D* histDataStruc = makeProjectionHist(DataStruc_0, "out", projQmax, projRange, qMin3D_0.size()-1);
  vector<std::unique_ptr<TLatex>> lexVec;
  std::vector<bool> reject_region = rej_OSL();


  std::vector<std::string> projections = {"out", "side", "long"};
  std::vector<std::string> projections_text = {"out", "side", "long"};
  if (!reject_final) {//add to all projection names "_{NOREJ}" 
    for (auto& proj : projections_text) {
      proj += "_{NOREJ}";
    }
  }
  std::vector<std::string> qTexts = {qOutText.Data(), qSideText.Data(), qLongText.Data()};

  for (size_t i = 0; i < projections.size(); ++i) {
    auto histDataStruc = makeProjectionHist(DataStruc, projections[i], projQmax, qMax, qMin3D_02);
    // histDataStruc->GetXaxis()->SetRangeUser(qMin, qMax);
    histDataStruc->SetLineColor(4);
    histDataStruc->SetMarkerColor(4);
    histDataStruc->SetMinimum(yMin);
    histDataStruc->SetMaximum(yMax);
    histDataStruc->SetTitle("");
    histDataStruc->GetXaxis()->SetTitle(("Q_{" + projections_text[i] + "} [GeV]").c_str());
    if (r2) {
      histDataStruc->GetYaxis()->SetTitle(("R_{2} (Q_{" + projections_text[i] + "})").c_str());
    } else {
      histDataStruc->GetYaxis()->SetTitle(("C_{2} (Q_{" + projections_text[i] + "})").c_str());
    }
    histDataStruc->GetXaxis()->SetTitleSize(0.05);
    histDataStruc->GetYaxis()->SetTitleSize(0.05);

    std::unique_ptr<TGraph> graphInit;
    std::unique_ptr<TGraph> graph;
    if (histData.find("_cube_") != std::string::npos) {
      if (useInit)
      //   graphInit = std::unique_ptr<TGraph>(makeProjection(hist, funcInit, projections[i], qMin, qMax, 1, 9));
      // graph = std::unique_ptr<TGraph>(makeProjection(hist, func, projections[i], qMin, qMax, 1, 9));
        graphInit = std::unique_ptr<TGraph>(makeProjection(hist, funcInit, projections[i], qMin, qMax, 1, 9));
      graph = std::unique_ptr<TGraph>(makeProjection(hist, func, projections[i], qMin, qMax, 1, 9));
    } else {
      if (useInit)
      //   graphInit = std::unique_ptr<TGraph>(makeProjection(hist, funcInit, projections[i], qMin, qMax, projMin, projMax));
      // graph = std::unique_ptr<TGraph>(makeProjection(hist, func, projections[i], qMin, qMax, projMin, projMax));
        graphInit = std::unique_ptr<TGraph>(makeProjection(hist, funcInit, projections[i], qMin, qMax, projMin, projMax));
      graph = std::unique_ptr<TGraph>(makeProjection(hist, func, projections[i], qMin, qMax, projMin, projMax));
    }
    if (useInit) {
      graphInit->SetLineColor(3);
      graphInit->SetLineWidth(2);
    }
    graph->SetLineColor(2);
    graph->SetLineWidth(2);
    auto textInit = std::make_unique<TPaveText>(.33, .63, .63, .93, "brNDC");
    std::unique_ptr<TPaveText> text;
    //add text "Final fit" above fit result of func as a separate text box
    auto textFinal = std::make_unique<TPaveText>(.63, .93, .93, .98, "brNDC");
    if (useInit) {
      textFinal->SetFillStyle(0);
      textFinal->SetFillColor(0);
      textFinal->SetBorderSize(0);
      textFinal->SetTextColor(1);
      textFinal->SetTextFont(42);
      textFinal->SetTextAlign(11);
      textFinal->SetTextSize(0.035);
      textFinal->AddText("Final fit");
      textInit->SetFillStyle(0);
      textInit->SetFillColor(0);
      textInit->SetBorderSize(0);
      textInit->SetTextColor(1);
      textInit->SetTextFont(42);
      textInit->SetTextAlign(11);
      textInit->SetTextSize(0.035);
      addFitResults(textInit.get(), funcInit);

      text = std::make_unique<TPaveText>(.63, .63, .93, .93, "brNDC");
      text->SetFillStyle(0);
      text->SetFillColor(0);
      text->SetBorderSize(0);
      text->SetTextColor(1);
      text->SetTextFont(42);
      text->SetTextAlign(11);
      text->SetTextSize(0.035);
      addFitResults(text.get(), func);
    } else {
      text = std::make_unique<TPaveText>(.53, .53, .93, .93, "brNDC");
      text->SetFillStyle(0);
      text->SetFillColor(0);
      text->SetBorderSize(0);
      text->SetTextColor(1);
      text->SetTextFont(42);
      text->SetTextAlign(11);
      text->SetTextSize(0.045);
      addFitResults(text.get(), func);
    }

    auto text2 = std::make_unique<TPaveText>(.53, .14, .93, .38, "brNDC");
    text2->SetFillStyle(0);
    text2->SetFillColor(0);
    text2->SetBorderSize(0);
    text2->SetTextColor(1);
    text2->SetTextFont(42);
    text2->SetTextAlign(11);
    text2->SetTextSize(0.045);
    for (const auto& comment : commentVec) {
      text2->AddText(comment.c_str());
    }
    text2->AddText(c2IndexText);
    if (histData.find("_cube_") != std::string::npos) {
      text2->AddText("0 #leq Q_{side} < 180 MeV");
      text2->AddText("0 #leq Q_{long} < 180 MeV");
    } else {
      text2->AddText(qTexts[(i + 1) % 3].c_str());
      text2->AddText(qTexts[(i + 2) % 3].c_str());
    }

    auto line = std::make_unique<TLine>(histDataStruc->GetXaxis()->GetXmin(), 1., histDataStruc->GetXaxis()->GetXmax(), 1.);
    line->SetLineStyle(7);

    histDataStruc->Draw();
    if (useInit)
      graphInit->Draw("Csame");
    graph->Draw("Csame");
    line->Draw();
    if (useInit) {
      textFinal->Draw();
      textInit->Draw();
    }
    text->Draw();
    text2->Draw();

    std::vector<std::unique_ptr<TArrow>> arrows; // Store arrows to manage their lifetime
    if (rejTo > rejFrom) {
      double ar_y = yMin + 0.14 * (yMax - yMin);
      arrows.push_back(std::make_unique<TArrow>(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>"));
      arrows.back()->Draw();
      lexc.SetTextSize(0.041);
      lexc.SetTextAlign(23);
      ar_y = yMin + 0.12 * (yMax - yMin);
      lexc.DrawLatex(0.5 * (rejFrom + rejTo), ar_y, "excluded");
    }

    for (size_t j = 0; j < reject_region.size(); ++j) {
      if (reject_region.at(j)) {
        double ar_y = yMin + 0.14 * (yMax - yMin) + 0.01*(reject_region.size()-1-j);
        if (rejToOsl[j][i] > 0.85)
          ar_y = yMin + 0.01;
        arrows.push_back(std::make_unique<TArrow>(rejFromOsl[j][i], ar_y, rejToOsl[j][i], ar_y, 0.015, "<>"));
        arrows.back()->Draw();
        lexVec.push_back(std::make_unique<TLatex>(0.5 * (rejFromOsl[j][i] + rejToOsl[j][i]), ar_y+0.01, ("rej " + std::to_string(j+1)).c_str()));
        lexVec.back()->SetTextSize(0.02);
        lexVec.back()->SetTextAlign(23);
        lexVec.back()->Draw();
      }
    }
    lexc.SetTextSize(0.041);
    lexc.SetTextAlign(23);
    ar_y = yMin + 0.12 * (yMax - yMin);
    // lexc.DrawLatex(0.5 * (rejFromOsl[j][i] + rejToOsl[j][i]), ar_y, "excluded");
    lexc.DrawLatex(qMin+0.3, ar_y, "excluded");

    for (const auto& format : formatVec) {
      std::string outFilePath = Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(), 
                                                            today.c_str(), 
                                                            c2_func.c_str(), 
                                                            format.c_str(), 
                                                            plotName.c_str(), 
                                                            projections_text[i].c_str(), 
                                                            format.c_str());
      canvas->Print(outFilePath.c_str());
    }
    canvas->Clear();
  }

  // std::unique_ptr<TH1D> histDataStruc = makeProjectionHist(DataStruc, "out", projQmax, qMax, qMin3D_02);
  // histDataStruc->GetXaxis()->SetRangeUser(qMin, qMax-qMin);
  // // histDataStruc->GetXaxis()->SetLimits(qMin, qMax);
  // // TH1D* projectionHist = makeProjectionHist(DataStruc, "out", projQmax, qMin, qMin3D_02.size()-1);
  // // for (int i = 1; i <= histDataStruc->GetNbinsX(); ++i) {
  // //   histDataStruc->SetBinContent(i, projectionHist->GetBinContent(i));
  // //   histDataStruc->SetBinError(i, projectionHist->GetBinError(i));
  // // }
  // histDataStruc->SetLineColor(4);
  // histDataStruc->SetMarkerColor(4);
  // // histX->SetLineColor(4);
  // // histX->SetMarkerColor(4);

  // // histX->SetMinimum(yMin);
  // // histX->SetMaximum(yMax);

  // // histX->SetTitle("");
  // // histX->GetXaxis()->SetTitle("Q_{out} [GeV]");
  // // if (r2) {
  // //   histX->GetYaxis()->SetTitle("R_{2} (Q_{out})");
  // // } else {
  // //   histX->GetYaxis()->SetTitle("C_{2} (Q_{out})");
  // // }

  // // histX->GetXaxis()->SetTitleSize(0.05);
  // // histX->GetYaxis()->SetTitleSize(0.05);

  // histDataStruc->SetMinimum(yMin);
  // histDataStruc->SetMaximum(yMax);

  // histDataStruc->SetTitle("");
  // histDataStruc->GetXaxis()->SetTitle("Q_{out} [GeV]");
  // if (r2) {
  //   histDataStruc->GetYaxis()->SetTitle("R_{2} (Q_{out})");
  // } else {
  //   histDataStruc->GetYaxis()->SetTitle("C_{2} (Q_{out})");
  // }

  // histDataStruc->GetXaxis()->SetTitleSize(0.05);
  // histDataStruc->GetYaxis()->SetTitleSize(0.05);

  // TGraph* graphX;
  // if (histData.find("_cube_") != std::string::npos) {
  //   graphX = makeProjection(hist, func, "out", qMin, qMax, 1, 9);
  // } else {
  //   graphX = makeProjection(hist, func, "out", qMin, qMax, projMin, projMax);
  // }
  // graphX->SetLineColor(2);
  // graphX->SetLineWidth(2);

  // TPaveText *textX = new TPaveText(.53, .53, .93, .93, "brNDC");
  // textX->SetFillStyle(0);
  // textX->SetFillColor(0);
  // textX->SetBorderSize(0);
  // textX->SetTextColor(1);
  // textX->SetTextFont(42);
  // textX->SetTextAlign(11);
  // textX->SetTextSize(0.045);
  // addFitResults(textX, func);

  // TPaveText *text2X = new TPaveText(.53, .14, .93, .38, "brNDC");
  // text2X->SetFillStyle(0);
  // text2X->SetFillColor(0);
  // text2X->SetBorderSize(0);
  // text2X->SetTextColor(1);
  // text2X->SetTextFont(42);
  // text2X->SetTextAlign(11);
  // text2X->SetTextSize(0.045);
  // for (size_t i = 0; i < commentVec.size(); ++i) {
  //   text2X->AddText(commentVec.at(i).c_str());
  // }
  // text2X->AddText(c2IndexText);
  // if (histData.find("_cube_") != std::string::npos) {
  //   text2X->AddText("0 #leq Q_{side} < 180 MeV");
  //   text2X->AddText("0 #leq Q_{long} < 180 MeV");
  // } else {
  //   text2X->AddText(qSideText.Data());
  //   text2X->AddText(qLongText.Data());
  // }

  // // TLine* lineX = new TLine(histX->GetXaxis()->GetXmin(), 1.,
  // //                          histX->GetXaxis()->GetXmax(), 1.);
  // TLine* lineX = new TLine(histDataStruc->GetXaxis()->GetXmin(), 1.,
  //                          histDataStruc->GetXaxis()->GetXmax(), 1.);
  // lineX->SetLineStyle(7);

  // // histX->Draw();
  // histDataStruc->Draw();
  // graphX->Draw("Csame");
  // lineX->Draw();
  // textX->Draw();
  // text2X->Draw();

  // if (rejTo > rejFrom) {
  //   ar_y = yMin + 0.14 * (yMax - yMin);
  //   TArrow *ar1 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
  //   ar1 -> Draw();

  //   lexc.SetTextSize(0.041);
  //   lexc.SetTextAlign(23);
  //   ar_y = yMin + 0.12 * (yMax - yMin);
  //   lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  // }
  // for (size_t i = 0; i < reject_region.size(); ++i) {
  //   if (reject_region.at(i)) {
  //     std::cout << "rejected region: " << rejFromOsl[i][0] << " - " << rejToOsl[i][0] << std::endl;
  //     ar_y = yMin + 0.14 * (yMax - yMin);
  //     TArrow *ar2 = new TArrow(rejFromOsl[i][0], ar_y, rejToOsl[i][0], ar_y, 0.015, "<>");
  //     ar2->Draw();

  //     lexc.SetTextSize(0.041);
  //     lexc.SetTextAlign(23);
  //     ar_y = yMin + 0.12 * (yMax - yMin);
  //     lexc.DrawLatex(0.5 * (rejFromOsl[i][0] + rejToOsl[i][0]) + 0.1, ar_y, "excluded");
  //   }
  // }

  // for (size_t i = 0; i < formatVec.size(); ++i) {
  //   outFilePath.Form("%s/%s/%s/%s/%s_out.%s", outFilePathBase.c_str(),
  //                                          today.c_str(),
  //                                          c2_func.c_str(),
  //                                          formatVec.at(i).c_str(),
  //                                          plotName.c_str(),
  //                                          formatVec.at(i).c_str());
  //   canvas->Print(outFilePath);
  // }

  // // delete histX;
  // delete graphX;
  // delete textX;
  // delete text2X;
  // delete lineX;

  // /**
  //  * Q_side
  //  */
  // TH1D* histY;
  // if (histData.find("_cube_") != std::string::npos) {
  //   TH1D* histYpart1 = dynamic_cast<TH1D*>(hist->ProjectionY(
  //                      (histName + "_side_part1").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histYpart2 = dynamic_cast<TH1D*>(hist->ProjectionY(
  //                      (histName + "_side_part2").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histYpart3 = dynamic_cast<TH1D*>(hist->ProjectionY(
  //                      (histName + "_side_part3").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histYpart4 = dynamic_cast<TH1D*>(hist->ProjectionY(
  //                      (histName + "_side_part4").c_str(),
  //                      1, 27, 1, 27));

  //   histY = new TH1D((histName + "_side").c_str(),
  //                    hist->GetTitle(),
  //                    varQbinsSize, &varQbins[0]);

  //   for (int i = 1; i <= 15; ++i) {
  //     histY->SetBinContent(i, histYpart1->GetBinContent(i) / (9 * 9));
  //     histY->SetBinError(i, histYpart1->GetBinError(i) / (9 * 9));
  //   }

  //   histY->SetBinContent(16, histYpart2->GetBinContent(17) / (3 * 3));
  //   histY->SetBinError(16, histYpart2->GetBinError(17) / (3 * 3));
  //   histY->SetBinContent(17, histYpart2->GetBinContent(20) / (3 * 3));
  //   histY->SetBinError(17, histYpart2->GetBinError(20) / (3 * 3));
  //   histY->SetBinContent(18, histYpart2->GetBinContent(23) / (3 * 3));
  //   histY->SetBinError(18, histYpart2->GetBinError(23) / (3 * 3));
  //   histY->SetBinContent(19, histYpart2->GetBinContent(26) / (3 * 3));
  //   histY->SetBinError(19, histYpart2->GetBinError(26) / (3 * 3));

  //   histY->SetBinContent(20, histYpart3->GetBinContent(32));
  //   histY->SetBinError(20, histYpart3->GetBinError(32));
  //   histY->SetBinContent(21, histYpart3->GetBinContent(41));
  //   histY->SetBinError(21, histYpart3->GetBinError(41));
  //   histY->SetBinContent(22, histYpart3->GetBinContent(50));
  //   histY->SetBinError(22, histYpart3->GetBinError(50));

  //   histY->SetBinContent(23, histYpart4->GetBinContent(68));
  //   histY->SetBinError(23, histYpart4->GetBinError(68));
  //   histY->SetBinContent(24, histYpart4->GetBinContent(95));
  //   histY->SetBinError(24, histYpart4->GetBinError(95));
  //   histY->SetBinContent(25, histYpart4->GetBinContent(122));
  //   histY->SetBinError(25, histYpart4->GetBinError(122));
  //   histY->SetBinContent(26, histYpart4->GetBinContent(149));
  //   histY->SetBinError(26, histYpart4->GetBinError(149));
  //   histY->SetBinContent(27, histYpart4->GetBinContent(176));
  //   histY->SetBinError(27, histYpart4->GetBinError(176));

  //   delete histYpart1;
  //   delete histYpart2;
  //   delete histYpart3;
  //   delete histYpart4;
  // } else {
  //   histY = dynamic_cast<TH1D*>(hist->ProjectionY(
  //                 (histName + "_side").c_str(),
  //                 projMin, projMax,
  //                 projMin, projMax));
  //   histY->Scale(1. / (projRange * projRange));
  // }

  // // histDataStruc = makeProjectionHist(DataStruc_0, "side", projQmax, projRange, qMin3D_0.size()-1);
  // histDataStruc = makeProjectionHist(DataStruc, "side", projQmax, qMax, qMin3D_02);
  // histDataStruc->SetLineColor(4);
  // histDataStruc->SetMarkerColor(4);
  // // histY->SetLineColor(4);
  // // histY->SetMarkerColor(4);

  // // histY->SetMinimum(yMin);
  // // histY->SetMaximum(yMax);

  // // histY->SetTitle("");
  // // histY->GetXaxis()->SetTitle("Q_{side} [GeV]");
  // // if (r2) {
  // //   histY->GetYaxis()->SetTitle("R_{2} (Q_{side})");
  // // } else {
  // //   histY->GetYaxis()->SetTitle("C_{2} (Q_{side})");
  // // }
  // // histY->GetXaxis()->SetTitleOffset(1.0);

  // // histY->GetXaxis()->SetTitleSize(0.05);
  // // histY->GetYaxis()->SetTitleSize(0.05);
  // histDataStruc->SetMinimum(yMin);
  // histDataStruc->SetMaximum(yMax);
  // histDataStruc->GetXaxis()->SetRangeUser(qMin, qMax-qMin);

  // histDataStruc->SetTitle("");
  // histDataStruc->GetXaxis()->SetTitle("Q_{side} [GeV]");
  // if (r2) {
  //   histDataStruc->GetYaxis()->SetTitle("R_{2} (Q_{side})");
  // } else {
  //   histDataStruc->GetYaxis()->SetTitle("C_{2} (Q_{side})");
  // }

  // histDataStruc->GetXaxis()->SetTitleSize(0.05);
  // histDataStruc->GetYaxis()->SetTitleSize(0.05);

  // TGraph* graphY;
  // if (histData.find("_cube_") != std::string::npos) {
  //   graphY = makeProjection(hist, func, "side", qMin, qMax, 1, 9);
  // } else {
  //   graphY = makeProjection(hist, func, "side", qMin, qMax, projMin, projMax);
  // }
  // graphY->SetLineColor(2);
  // graphY->SetLineWidth(2);

  // TPaveText *textY = new TPaveText(.53, .53, .93, .93, "brNDC");
  // textY->SetFillStyle(0);
  // textY->SetFillColor(0);
  // textY->SetBorderSize(0);
  // textY->SetTextColor(1);
  // textY->SetTextFont(42);
  // textY->SetTextAlign(11);
  // textY->SetTextSize(0.045);
  // addFitResults(textY, func);

  // TPaveText *text2Y = new TPaveText(.53, .14, .93, .38, "brNDC");
  // text2Y->SetFillStyle(0);
  // text2Y->SetFillColor(0);
  // text2Y->SetBorderSize(0);
  // text2Y->SetTextColor(1);
  // text2Y->SetTextFont(42);
  // text2Y->SetTextAlign(11);
  // text2Y->SetTextSize(0.045);
  // for (size_t i = 0; i < commentVec.size(); ++i) {
  //   text2Y->AddText(commentVec.at(i).c_str());
  // }
  // text2Y->AddText(c2IndexText);
  // if (histData.find("_cube_") != std::string::npos) {
  //   text2Y->AddText("0 #leq Q_{out} < 180 MeV");
  //   text2Y->AddText("0 #leq Q_{long} < 180 MeV");
  // } else {
  //   text2Y->AddText(qOutText.Data());
  //   text2Y->AddText(qLongText.Data());
  // }

  // TLine* lineY = new TLine(histY->GetXaxis()->GetXmin(), 1.,
  //                          histY->GetXaxis()->GetXmax(), 1.);
  // lineY->SetLineStyle(7);

  // // histY->Draw();
  // histDataStruc->Draw();
  // graphY->Draw("Csame");
  // lineY->Draw();
  // textY->Draw();
  // text2Y->Draw();

  // if (rejTo > rejFrom) {
  //   ar_y = yMin + 0.14 * (yMax - yMin);
  //   TArrow *ar2 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
  //   ar2 -> Draw();

  //   lexc.SetTextSize(0.041);
  //   lexc.SetTextAlign(23);
  //   ar_y = yMin + 0.12 * (yMax - yMin);
  //   lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  // }
  // for (size_t i = 0; i < reject_region.size(); ++i) {
  //   if (reject_region.at(i)) {
  //     ar_y = yMin + 0.14 * (yMax - yMin);
  //     TArrow *ar2 = new TArrow(rejFromOsl[i][1], ar_y, rejToOsl[i][1], ar_y, 0.015, "<>");
  //     ar2->Draw();

  //     lexc.SetTextSize(0.041);
  //     lexc.SetTextAlign(23);
  //     ar_y = yMin + 0.12 * (yMax - yMin);
  //     lexc.DrawLatex(0.5 * (rejFromOsl[i][1] + rejToOsl[i][1]) + 0.1, ar_y, "excluded");
  //   }
  // }

  // for (size_t i = 0; i < formatVec.size(); ++i) {
  //   outFilePath.Form("%s/%s/%s/%s/%s_side.%s", outFilePathBase.c_str(),
  //                                           today.c_str(),
  //                                           c2_func.c_str(),
  //                                           formatVec.at(i).c_str(),
  //                                           plotName.c_str(),
  //                                           formatVec.at(i).c_str());
  //   canvas->Print(outFilePath);
  // }

  // delete histY;
  // delete graphY;
  // delete textY;
  // delete text2Y;
  // delete lineY;

  // /**
  //  * Q_long
  //  */
  // TH1D* histZ;
  // if (histData.find("_cube_") != std::string::npos) {
  //   TH1D* histZpart1 = dynamic_cast<TH1D*>(hist->ProjectionZ(
  //                      (histName + "_long_part1").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histZpart2 = dynamic_cast<TH1D*>(hist->ProjectionZ(
  //                      (histName + "_long_part2").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histZpart3 = dynamic_cast<TH1D*>(hist->ProjectionZ(
  //                      (histName + "_long_part3").c_str(),
  //                      1, 9, 1, 9));
  //   TH1D* histZpart4 = dynamic_cast<TH1D*>(hist->ProjectionZ(
  //                      (histName + "_long_part4").c_str(),
  //                      1, 27, 1, 27));

  //   histZ = new TH1D((histName + "_long").c_str(),
  //                    hist->GetTitle(),
  //                    varQbinsSize, &varQbins[0]);

  //   for (int i = 1; i <= 15; ++i) {
  //     histZ->SetBinContent(i, histZpart1->GetBinContent(i) / (9 * 9));
  //     histZ->SetBinError(i, histZpart1->GetBinError(i) / (9 * 9));
  //   }

  //   histZ->SetBinContent(16, histZpart2->GetBinContent(17) / (3 * 3));
  //   histZ->SetBinError(16, histZpart2->GetBinError(17) / (3 * 3));
  //   histZ->SetBinContent(17, histZpart2->GetBinContent(20) / (3 * 3));
  //   histZ->SetBinError(17, histZpart2->GetBinError(20) / (3 * 3));
  //   histZ->SetBinContent(18, histZpart2->GetBinContent(23) / (3 * 3));
  //   histZ->SetBinError(18, histZpart2->GetBinError(23) / (3 * 3));
  //   histZ->SetBinContent(19, histZpart2->GetBinContent(26) / (3 * 3));
  //   histZ->SetBinError(19, histZpart2->GetBinError(26) / (3 * 3));

  //   histZ->SetBinContent(20, histZpart3->GetBinContent(32));
  //   histZ->SetBinError(20, histZpart3->GetBinError(32));
  //   histZ->SetBinContent(21, histZpart3->GetBinContent(41));
  //   histZ->SetBinError(21, histZpart3->GetBinError(41));
  //   histZ->SetBinContent(22, histZpart3->GetBinContent(50));
  //   histZ->SetBinError(22, histZpart3->GetBinError(50));

  //   histZ->SetBinContent(23, histZpart4->GetBinContent(68));
  //   histZ->SetBinError(23, histZpart4->GetBinError(68));
  //   histZ->SetBinContent(24, histZpart4->GetBinContent(95));
  //   histZ->SetBinError(24, histZpart4->GetBinError(95));
  //   histZ->SetBinContent(25, histZpart4->GetBinContent(122));
  //   histZ->SetBinError(25, histZpart4->GetBinError(122));
  //   histZ->SetBinContent(26, histZpart4->GetBinContent(149));
  //   histZ->SetBinError(26, histZpart4->GetBinError(149));
  //   histZ->SetBinContent(27, histZpart4->GetBinContent(176));
  //   histZ->SetBinError(27, histZpart4->GetBinError(176));

  //   delete histZpart1;
  //   delete histZpart2;
  //   delete histZpart3;
  //   delete histZpart4;
  // } else {
  //   histZ = dynamic_cast<TH1D*>(hist->ProjectionZ(
  //                 (histName + "_long").c_str(),
  //                 projMin, projMax,
  //                 projMin, projMax));
  //   histZ->Scale(1. / (projRange * projRange));
  // }

  // // histDataStruc = makeProjectionHist(DataStruc_0, "long", projQmax, projRange, qMin3D_0.size()-1);
  // histDataStruc = makeProjectionHist(DataStruc, "long", projQmax, qMax, qMin3D_02);
  // histDataStruc->SetLineColor(4);
  // histDataStruc->SetMarkerColor(4);
  // // histZ->SetLineColor(4);
  // // histZ->SetMarkerColor(4);

  // // histZ->SetMinimum(yMin);
  // // histZ->SetMaximum(yMax);

  // // histZ->SetTitle("");
  // // histZ->GetXaxis()->SetTitle("Q_{long} [GeV]");
  // // if (r2) {
  // //   histZ->GetYaxis()->SetTitle("R_{2} (Q_{long})");
  // // } else {
  // //   histZ->GetYaxis()->SetTitle("C_{2} (Q_{long})");
  // // }

  // // histZ->GetXaxis()->SetTitleSize(0.05);
  // // histZ->GetYaxis()->SetTitleSize(0.05);
  // histDataStruc->SetMinimum(yMin);
  // histDataStruc->SetMaximum(yMax);
  // histDataStruc->GetXaxis()->SetRangeUser(qMin, qMax-qMin);
  
  // histDataStruc->SetTitle("");
  // histDataStruc->GetXaxis()->SetTitle("Q_{long} [GeV]");
  // if (r2) {
  //   histDataStruc->GetYaxis()->SetTitle("R_{2} (Q_{long})");
  // } else {
  //   histDataStruc->GetYaxis()->SetTitle("C_{2} (Q_{long})");
  // }

  // histDataStruc->GetXaxis()->SetTitleSize(0.05);
  // histDataStruc->GetYaxis()->SetTitleSize(0.05);
  // TGraph* graphZ;
  // if (histData.find("_cube_") != std::string::npos) {
  //   graphZ = makeProjection(hist, func, "long", qMin, qMax, projMin, projMax);
  // } else {
  //   graphZ = makeProjection(hist, func, "long", qMin, qMax, projMin, projMax);
  // }
  // graphZ->SetLineColor(2);
  // graphZ->SetLineWidth(2);

  // TPaveText *textZ = new TPaveText(.53, .53, .93, .93, "brNDC");
  // textZ->SetFillStyle(0);
  // textZ->SetFillColor(0);
  // textZ->SetBorderSize(0);
  // textZ->SetTextColor(1);
  // textZ->SetTextFont(42);
  // textZ->SetTextAlign(11);
  // textZ->SetTextSize(0.045);
  // addFitResults(textZ, func);

  // TPaveText *text2Z = new TPaveText(.53, .14, .93, .38, "brNDC");
  // text2Z->SetFillStyle(0);
  // text2Z->SetFillColor(0);
  // text2Z->SetBorderSize(0);
  // text2Z->SetTextColor(1);
  // text2Z->SetTextFont(42);
  // text2Z->SetTextAlign(11);
  // text2Z->SetTextSize(0.045);
  // for (size_t i = 0; i < commentVec.size(); ++i) {
  //   text2Z->AddText(commentVec.at(i).c_str());
  // }
  // text2Z->AddText(c2IndexText);
  // if (histData.find("_cube_") != std::string::npos) {
  //   text2Z->AddText("0 #leq Q_{out} < 180 MeV");
  //   text2Z->AddText("0 #leq Q_{side} < 180 MeV");
  // } else {
  //   text2Z->AddText(qOutText.Data());
  //   text2Z->AddText(qSideText.Data());
  // }

  // TLine* lineZ = new TLine(histZ->GetXaxis()->GetXmin(), 1.,
  //                          histZ->GetXaxis()->GetXmax(), 1.);
  // lineZ->SetLineStyle(7);

  // // histZ->Draw();
  // histDataStruc->Draw();
  // graphZ->Draw("Csame");
  // lineZ->Draw();
  // textZ->Draw();
  // text2Z->Draw();

  // if (rejTo > rejFrom) {
  //   ar_y = yMin + 0.14 * (yMax - yMin);
  //   TArrow *ar3 = new TArrow(rejFrom, ar_y, rejTo, ar_y, 0.015, "<>");
  //   ar3 -> Draw();

  //   lexc.SetTextSize(0.041);
  //   lexc.SetTextAlign(23);
  //   ar_y = yMin + 0.12 * (yMax - yMin);
  //   lexc.DrawLatex(0.5 * (rejFrom + rejTo) - 0.1, ar_y, "excluded");
  // }
  // for (size_t i = 0; i < reject_region.size(); ++i) {
  //   if (reject_region.at(i)) {
  //     ar_y = yMin + 0.14 * (yMax - yMin);
  //     TArrow *ar3 = new TArrow(rejFromOsl[i][2], ar_y, rejToOsl[i][2], ar_y, 0.015, "<>");
  //     ar3->Draw();

  //     lexc.SetTextSize(0.041);
  //     lexc.SetTextAlign(23);
  //     ar_y = yMin + 0.12 * (yMax - yMin);
  //     lexc.DrawLatex(0.5 * (rejFromOsl[i][2] + rejToOsl[i][2]) + 0.1, ar_y, "excluded");
  //   }
  // }

  // for (size_t i = 0; i < formatVec.size(); ++i) {
  //   outFilePath.Form("%s/%s/%s/%s/%s_long.%s", outFilePathBase.c_str(),
  //                                           today.c_str(),
  //                                           c2_func.c_str(),
  //                                           formatVec.at(i).c_str(),
  //                                           plotName.c_str(),
  //                                           formatVec.at(i).c_str());
  //   canvas->Print(outFilePath);
  // }

  // delete histZ;
  // delete graphZ;
  // delete textZ;
  // delete text2Z;
  // delete lineZ;
  // delete canvas;

  /**
   * 2D Plots
   */
  std::string ifReject = "";
  if (!reject_final) {
    ifReject = "_{NOREJ}";
  }
  // Q_out_side
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          ("out_side" + ifReject).c_str(),
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "out_side", commentVec, outFilePath,
                     projMin, projMax);
    size_t pos = outFilePath.Last('.');
    if (pos != kNPOS) { // kNPOS is used to check if the character was not found
        outFilePath.Insert(pos, "_L");
    }
    plot2Dprojection(histLarge, func, "out_side", commentVec, outFilePath,
                     1, 1);
  }

  // Q_side_long
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          ("side_long" + ifReject).c_str(),
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "side_long", commentVec, outFilePath,
                     projMin, projMax);
    size_t pos = outFilePath.Last('.');
    if (pos != kNPOS) { // kNPOS is used to check if the character was not found
        outFilePath.Insert(pos, "_L");
    }
    plot2Dprojection(histLarge, func, "side_long", commentVec, outFilePath,
                     1, 1);
  }

  // Q_out_long
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          plotName.c_str(),
                                          ("out_long" + ifReject).c_str(),
                                          formatVec.at(i).c_str());
    plot2Dprojection(hist, func, "out_long", commentVec, outFilePath,
                     projMin, projMax);
    size_t pos = outFilePath.Last('.');
    if (pos != kNPOS) { // kNPOS is used to check if the character was not found
        outFilePath.Insert(pos, "_L");
    }
    plot2Dprojection(histLarge, func, "out_long", commentVec, outFilePath,
                     1, 1);
  }

  // Correlation matrix
  for (size_t i = 0; i < formatVec.size(); ++i) {
    outFilePath.Form("%s/%s/%s/%s/%s/%s_%s.%s", outFilePathBase.c_str(),
                                          today.c_str(),
                                          c2_func.c_str(),
                                          formatVec.at(i).c_str(),
                                          "correlation",
                                          plotName.c_str(),
                                          ("correlation" + ifReject).c_str(),
                                          formatVec.at(i).c_str());
    plotCorrelationMatrix(fitData, outFilePath);
  }
}