// std
#include <iostream>
#include <fstream>
#include <iomanip>
// ROOT
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH3.h>
#include <TH2.h>
#include <TF3.h>
#include <TList.h>
#include <TArrayD.h>
#include <TCanvas.h>
#include <TPaveText.h>
// Fit
#include "utils.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"


/**
 * Custom defined chi2 function for fitting 3D histograms.
 */

// CustomChi2FCN::CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D)
//     : fData(data), fFunc(func), fQMin3D(qMin3D), fQNBins3D(qNBins3D),
//       fWrappedFunc(std::make_unique<ROOT::Math::WrappedMultiTF1>(fFunc, 3)),
//       fIntegrator(ROOT::Math::IntegrationMultiDim::kADAPTIVE) {}

// double CustomChi2FCN::operator()(const double* params) {
//     fFunc.SetParameters(params);
//     double chi2 = 0.0;

//     for (unsigned int i = 0; i < fData.Size(); ++i) {
//         double value, invError;
//         const double* coords = fData.GetPoint(i, value, invError);
//         double binWidth = computeBinWidth(coords);

//         double binMin[3], binMax[3];
//         setBinRange(coords, binWidth, binMin, binMax);

//         double integral = fIntegrator.Integral(*fWrappedFunc, binMin, binMax);
//         double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
//         double meanValue = integral / binVolume;

//         double diff = value - meanValue;
//         chi2 += (diff * diff) * (invError * invError);
//     }

//     return chi2;
// }

// double CustomChi2FCN::computeBinWidth(const double* coords) {
//     double maxCoord = *std::max_element(coords, coords + 3);
//     double binWidth = 0.0;
//     for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
//         if (fQMin3D[k] <= maxCoord && maxCoord < fQMin3D[k + 1]) {
//             binWidth = (fQMin3D[k + 1] - fQMin3D[k]) / fQNBins3D[k];
//             break;
//         }
//     }
//     return binWidth;
// }

// void CustomChi2FCN::setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) {
//     for (int j = 0; j < 3; ++j) {
//         binMin[j] = coords[j] - binWidth / 2.0;
//         binMax[j] = coords[j] + binWidth / 2.0;
//     }
// }

// Constructor
CustomChi2FCN::CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D)
    : fData(data), fFunc(func), fQMin3D(qMin3D), fQNBins3D(qNBins3D),
      fWrappedFunc(std::make_unique<ROOT::Math::WrappedMultiTF1>(fFunc, 3)),
      fIntegrator(ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
        std::cout << "CustomChi2FCN constructor" << std::endl;
      }

// Function evaluation
double CustomChi2FCN::operator()(const double* params) {
    fFunc.SetParameters(params);
    fWrappedFunc->SetParameters(params);
    //cout the parameters
    // for (int i = 0; i < fFunc.GetNpar(); i++) {
    //     std::cout << "Par " << i << " = " << params[i]<<";";
    // }
    double chi2 = 0.0;

    for (unsigned int i = 0; i < fData.Size(); ++i) {
        double value, invError;
        const double* coords = fData.GetPoint(i, value, invError);
        double binWidth = computeBinWidth(coords);

        double binMin[3], binMax[3];
        setBinRange(coords, binWidth, binMin, binMax);

        double integral = fIntegrator.Integral(*fWrappedFunc, binMin, binMax);
        double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
        double meanValue = integral / binVolume;

        double diff = value - meanValue;
        // std::cout << "binvalue: "<< value << " meanvalue: "<< meanValue << " diff: "<< diff << " invError: "<< invError << std::endl;
        chi2 += (diff * diff) * (invError * invError);
    }
    std::cout << std::fixed << std::setprecision(6);
    // std::cout <<"chi2: "<< chi2<< std::endl;
    return chi2;
}

// Helper methods
double CustomChi2FCN::computeBinWidth(const double* coords) const {
    double maxCoord = *std::max_element(coords, coords + 3);
    double binWidth = 0.0;
    for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
        if (fQMin3D[k] <= maxCoord && maxCoord < fQMin3D[k + 1]) {
            binWidth = (fQMin3D[k + 1] - fQMin3D[0]) / fQNBins3D[k];
            break;
        }
    }
    return binWidth;
}

void CustomChi2FCN::setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) const {
    for (int j = 0; j < 3; ++j) {
        binMin[j] = coords[j] - binWidth / 2.0;
        binMax[j] = coords[j] + binWidth / 2.0;
    }
}

// CustomChi2FCN::CustomChi2FCN(const ROOT::Fit::BinData& data, const IModelFunction& func)
//     : BaseObjFunction(func.NPar(), data.Size()), fData(data), fFunc(func), fNEffPoints(0) {
//       std::cout << "CustomChi2FCN constructor" << std::endl;
// }

// ROOT::Math::BasicFitMethodFunction<ROOT::Math::IMultiGenFunction>::BaseFunction* CustomChi2FCN::Clone() const {
//     return new CustomChi2FCN(fData, fFunc);
// }

// double CustomChi2FCN::DoEval(const double* x) const {
//     this->UpdateNCalls();
//     // if (!fData.HaveCoordErrors()) {
//     // if (fData.HaveCoordErrors()) {
//         // return ROOT::Fit::FitUtil::EvaluateChi2(fFunc, fData, x, fNEffPoints, ROOT::Fit::kDefault);
//     // } else {
//         return ROOT::Fit::FitUtil::EvaluateChi2Effective(fFunc, fData, x, fNEffPoints);
//     // }
// }

// double CustomChi2FCN::DataElement(const double* x, unsigned int i, double* g, double* h, bool fullHessian) const {
//     return ROOT::Fit::FitUtil::EvaluateChi2Residual(fFunc, fData, x, i, g);
// }

// CustomChi2FCN::CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& binEdges, const std::vector<int>& binCounts)
//     : fData(data), fFunc(func), fBinEdges(binEdges), fBinCounts(binCounts) {
//     fWrappedFunc = std::make_unique<ROOT::Math::WrappedMultiTF1>(fFunc, 3);
//     fIntegrator.SetFunction(*fWrappedFunc);
// }

// unsigned int CustomChi2FCN::NDim() const {
//     return fFunc.GetNpar();
// }

// double CustomChi2FCN::DoEval(const double* params) const {
//     double chi2 = 0.0;

//     std::cout << "DoEval, before loop" << std::endl;
//     for (unsigned int i = 0; i < fData.Size(); ++i) {
//         double value, invError;
//         const double* coords = fData.GetPoint(i, value, invError);

//         double binMin[3], binMax[3];
//         double binWidth = determineBinWidth(coords);
//         for (int j = 0; j < 3; ++j) {
//             binMin[j] = coords[j] - binWidth / 2.0;
//             binMax[j] = coords[j] + binWidth / 2.0;
//         }

//         std::cout << "Bin " << i << ": coords = (" << coords[0] << ", " << coords[1] << ", " << coords[2] << "), "
//                   << "value = " << value << ", invError = " << invError << std::endl;

//         fFunc.SetParameters(params);
//         double integral = fIntegrator.Integral(binMin, binMax);
//         double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
//         double meanValue = integral / binVolume;

//         std::cout << "Integral = " << integral << ", binVolume = " << binVolume << ", meanValue = " << meanValue << std::endl;

//         double diff = value - meanValue;
//         chi2 += (diff * diff) * (invError * invError);
//     }

//     std::cout << "DoEval, after loop, chi2 = " << chi2 << std::endl;
//     return chi2;
// }

// CustomChi2FCN* CustomChi2FCN::Clone() const {
//     return new CustomChi2FCN(fData, fFunc, fBinEdges, fBinCounts);
// }

// double CustomChi2FCN::determineBinWidth(const double* coords) const {
//     double maxCoord = *std::max_element(coords, coords + 3);
//     for (size_t i = 0; i < fBinEdges.size() - 1; ++i) {
//         if (fBinEdges[i] <= maxCoord && fBinEdges[i + 1] > maxCoord) {
//             return (fBinEdges[i + 1] - fBinEdges[i]) / fBinCounts[i];
//         }
//     }
//     return 0.0;
// }
// CustomChi2FCN::CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D)
//     : fData(data), fFunc(func), fQMin3D(qMin3D), fQNBins3D(qNBins3D),
//       fWrappedFunc(std::make_unique<ROOT::Math::WrappedMultiTF1>(fFunc, 3)),
//       fIntegrator(ROOT::Math::IntegrationMultiDim::kADAPTIVE) {
//         std::cout << "CustomChi2FCN constructor" << std::endl;
//       }

// double CustomChi2FCN::DoEval(const double* params) const {
//     // fFunc.SetParameters(params);
//     double chi2 = 0.0;
//     std::cout << "DoEval, before loop" << std::endl;
//     for (unsigned int i = 0; i < fData.Size(); ++i) {
//         double value, invError;
//         const double* coords = fData.GetPoint(i, value, invError);
//         double binWidth = computeBinWidth(coords);

//         double binMin[3], binMax[3];
//         setBinRange(coords, binWidth, binMin, binMax);

//         std::cout << "Bin " << i << ": coords = (" << coords[0] << ", " << coords[1] << ", " << coords[2] << "), "
//                   << "value = " << value << ", invError = " << invError << std::endl;

//         // double integral = fIntegrator.Integral(*fWrappedFunc, 3, binMin, binMax);
//         double integral = const_cast<ROOT::Math::IntegratorMultiDim&>(fIntegrator).Integral(*fWrappedFunc, 3, binMin, binMax);
//         double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
//         double meanValue = integral / binVolume;

//         std::cout << "Integral = " << integral << ", binVolume = " << binVolume << ", meanValue = " << meanValue << std::endl;

//         double diff = value - meanValue;
//         chi2 += (diff * diff) * (invError * invError);
//     }

//     return chi2;
// }

// unsigned int CustomChi2FCN::NDim() const {
//     return fFunc.GetNpar();
// }

// ROOT::Math::IBaseFunctionMultiDim* CustomChi2FCN::Clone() const {
//     return new CustomChi2FCN(fData, fFunc, fQMin3D, fQNBins3D);
// }

// double CustomChi2FCN::computeBinWidth(const double* coords) const {
//     double maxCoord = *std::max_element(coords, coords + 3);
//     double binWidth = 0.0;
//     for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
//         if (fQMin3D[k] <= maxCoord && maxCoord < fQMin3D[k + 1]) {
//             binWidth = (fQMin3D[k + 1] - fQMin3D[k]) / fQNBins3D[k];
//             break;
//         }
//     }
//     return binWidth;
// }

// void CustomChi2FCN::setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) const {
//     for (int j = 0; j < 3; ++j) {
//         binMin[j] = coords[j] - binWidth / 2.0;
//         binMax[j] = coords[j] + binWidth / 2.0;
//     }
// }


/**
 * \brief Make projection of 3D histogram fitting function to one dimension.
 */
TGraph* makeProjection(TH3D* hist, TF3* fit, const string& component,
                       double qMin, double qMax,
                       int projMin, int projMax) {
  bool projOut = false;
  bool projSide = false;
  bool projLong = false;

  vector<double> qBoud;
  vector<double> x;
  vector<double> y;
  size_t nSamples = 100;
  double qStep;
  double projQmin;
  double projQmax;
  double projQrange;
  //TF3* fit;
  TGraph* graph;

  if (component.compare("out") == 0) {
    projOut = true;
  } else if (component.compare("side") == 0) {
    projSide = true;
  } else if (component.compare("long") == 0) {
    projLong = true;
  } else {
    throw "ERROR: Projection not found!";
  }
  // std::cout<<"QMIN QMAX "<<qMin<<"\t"<<qMax<<std::endl;
  if(qMax < qMin) {
    double qmax = qMax;
    qMax = qMin;
    qMin = qmax;
  }
  qStep = (qMax - qMin) / nSamples;
  for (size_t i = 0; i < nSamples; ++i) {
    qBoud.emplace_back(qMin + i * qStep);
  }
  qBoud.emplace_back(qMax);
  // for (size_t i = 0; i < nSamples; ++i) {
  //   std::cout<<i<<"\t"<<qBoud.at(i)<<"\t";
  // }

  projQmin = hist->GetXaxis()->GetBinLowEdge(projMin);
  projQmax = hist->GetXaxis()->GetBinLowEdge(projMax) +
             hist->GetXaxis()->GetBinWidth(projMax);
  projQrange = projQmax - projQmin;

/*
  ROOT::Math::Integrator fit;
  fit.SetFunction(*dynamic_cast<TF3*>(hist->GetListOfFunctions()->First()));
  fit.SetRelTolerance(0.01);
  double low_lim,up_lim;
  int low, up, increment, i = 1;
  if (qBoud.at(0)<qBoud.at(qBoud.size() - 1)) {
    low_lim = qBoud.at(0);
    up_lim = qBoud.at(qBoud.size() - 1);
    low = 0;
    up = qBoud.size() - 1;
    increment = 1;
    i = low;
  } else {
    low_lim = qBoud.at(qBoud.size() - 1);
    up_lim = qBoud.at(0);
    low = qBoud.size() - 1;
    up = 0;
    increment = -1;
    i = low;
  }
TH1F *h1 = new TH1F("h1", "", 100, projQmin, projQmax);
  if (projOut) {
      std::cout<<"low "<<low<<" up "<<up<<" increment "<<increment<<" inc*i "<<increment*i<<std::endl;
      for (i; increment*i < up; i += increment) {
        std::cout<<"low "<<low<<" up "<<up<<" increment "<<increment<<" qBoud.at(i) "<<qBoud.at(i)<<" qBoud.at(i + increment) "<<qBoud.at(i + increment)<<" projmin "<<projQmin<<" projmax "<< projQmax<<std::endl;
        x.emplace_back((qBoud.at(i + increment) + qBoud.at(i)) / 2);
        y.emplace_back(fit.Integral({qBoud.at(i), qBoud.at(i + increment),
                                     projQmin, projQmax,
        projQmin, projQmax}) /
                                         (qStep * projQrange * projQrange));
        std::cout<<"x "<<x.at(0)<<" y "<<y.at(0)<<std::endl;
        
        for (int j = 1; j <= 100; j++) {
    double x = qBoud.at(i) + j * (qBoud.at(i + increment) - qBoud.at(i)) / 100.0;
    double y = projQmin + j * (projQmax - projQmin) / 100.0;
    double z = projQmin + j * (projQmax - projQmin) / 100.0;
    double val = fit.Eval(x, y, z);
    h1->Fill(val);
    std::cout << "x = " << x << ", y = " << y << ", z = " << z << ", f(x,y,z) = " << val << std::endl;
}

// draw histogram of evaluated function values
TCanvas *c15 = new TCanvas("c15", "", 800, 600);
h1->Draw();
        // Evaluate fit function at 100 points in integration limits
const int num_points = 100; // Number of points to generate
double x, y, z;

for (int j = 0; j < num_points; j++) {
    x = qBoud.at(i) + j * (qBoud.at(i + increment) - qBoud.at(i)) / (num_points - 1);
    for (int k = 0; k < num_points; k++) {
        y = projQmin + k * (projQmax - projQmin) / (num_points - 1);
        z = projQmin + k * (projQmax - projQmin) / (num_points - 1);
        double result = fit.Eval(x, y, z);
        std::cout << "x = " << x << ", y = " << y << ", z = " << z << ", f(x,y,z) = " << result << std::endl;
    }
}


      }
    } 
    
    */
/*
else if (projSide) {
      for (i; increment*i < up; i += increment) {
        std::cout<<"HEREEEEE"<<qBoud.at(i)<<"\t"<<qBoud.at(i+increment)<<"\t"<<projQmin<<"\t"<<projQmax<<std::endl;
        x.emplace_back((qBoud.at(i + increment) + qBoud.at(i)) / 2);
        y.emplace_back(fit.Integral({projQmin, projQmax,
                                      qBoud.at(i), qBoud.at(i + increment),
        projQmin, projQmax}) /
                                          (qStep * projQrange * projQrange));
      }
    } else if (projLong) {
      for (i; increment*i < up; i += increment) {
      //      qBoud.at(i) -= 1e-9;
      //  qBoud.at(i + 1) += 1e-9;
        std::cout<<"HEREEEEE"<<qBoud.at(i)<<"\t"<<qBoud.at(i+increment)<<"\t"<<projQmin<<"\t"<<projQmax<<std::endl;
        x.emplace_back((qBoud.at(i + increment) + qBoud.at(i)) / 2);
        y.emplace_back(fit.Integral({projQmin, projQmax,
                                     projQmin, projQmax,
        qBoud.at(i), qBoud.at(i + increment)}) /
                                         (qStep * projQrange * projQrange));
      }
    }*/


  //fit = dynamic_cast<TF3*>(hist->GetListOfFunctions()->First());
  
  if (projOut) {
    for (size_t i = 0; i < (qBoud.size() - 1); ++i) {
      x.emplace_back((qBoud.at(i + 1) + qBoud.at(i)) / 2.0);
      y.emplace_back(fit->Integral(qBoud.at(i), qBoud.at(i + 1),
                                   projQmin, projQmax,
                                   projQmin, projQmax) /
                                       (1.0* qStep * projQrange * projQrange));
    }
  } else if (projSide) {
    for (size_t i = 0; i < (qBoud.size() - 1); ++i) {
      x.emplace_back((qBoud.at(i + 1) + qBoud.at(i)) / 2.0);
      y.emplace_back(fit->Integral(projQmin, projQmax,
                                    qBoud.at(i), qBoud.at(i + 1),
                                    projQmin, projQmax) /
                                        (1.0* qStep * projQrange * projQrange));
    }
  } else if (projLong) {
    for (size_t i = 0; i < (qBoud.size() - 1); ++i) {
      x.emplace_back((qBoud.at(i + 1) + qBoud.at(i)) / 2.0);
      y.emplace_back(fit->Integral(projQmin, projQmax,
                                   projQmin, projQmax,
                                   qBoud.at(i), qBoud.at(i + 1)) /
                                       (1.0* qStep * projQrange * projQrange));
    }
  }
    
  // std::cout<<"projQmin "<<projQmin<<" projQmax "<<projQmax<<std::endl;
  // std::cout<<"x"<<std::endl;
  // std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, " / "));
  // std::cout<<std::endl<<"y"<<std::endl;
  // std::copy(y.begin(), y.end(), std::ostream_iterator<double>(std::cout, " / "));
  // std::cout<<std::endl;
  graph = new TGraph(nSamples, &x[0], &y[0]);

  return graph;
}

TGraph2D* make2DProjection(TF3* fit, const std::string& axes, double qMin, double qMax, double projQmin, double projQmax) {
  int nSamples = 15;
  double projQrange = projQmax - projQmin;
  // Define qStep and qBoud
  double qStep = (qMax - qMin) / nSamples;
  std::vector<double> qBoud;
  for (int i = 0; i <= nSamples; ++i) {
      qBoud.emplace_back(qMin + i * qStep);
  }

  // Create a TGraph2D to store the function values
  TGraph2D* graph2D = new TGraph2D();
  int pointIndex = 0;
  // Determine the projection axes
  std::string xAxisTitle, yAxisTitle;
  if (axes == "out_side") {
      xAxisTitle = "Q_{out} [GeV]";
      yAxisTitle = "Q_{side} [GeV]";
      for (int i = 0; i < qBoud.size() - 1; ++i) {
          for (int j = 0; j < qBoud.size() - 1; ++j) {
              double integral = fit->Integral(qBoud.at(i), qBoud.at(i + 1),
                                      qBoud.at(j), qBoud.at(j + 1),
                                      projQmin, projQmax) / 
                                      (1.0 * qStep * qStep * projQrange);
              double xCenter = (qBoud.at(i) + qBoud.at(i + 1)) / 2.0;
              double yCenter = (qBoud.at(j) + qBoud.at(j + 1)) / 2.0;
              // std::cout<< "Qout "<<qBoud.at(i)<<" Qside "<<qBoud.at(j)<<" integral "<<integral<<std::endl;
              graph2D->SetPoint(pointIndex++, xCenter, yCenter, integral);
          }
      }
  } else if (axes == "side_long") {
      xAxisTitle = "Q_{side} [GeV]";
      yAxisTitle = "Q_{long} [GeV]";
      for (int i = 0; i < qBoud.size() - 1; ++i) {
          for (int j = 0; j < qBoud.size() - 1; ++j) {
              double integral = fit->Integral(projQmin, projQmax,
                                      qBoud.at(i), qBoud.at(i + 1),
                                      qBoud.at(j), qBoud.at(j + 1)) / 
                                      (1.0 * projQrange * qStep * qStep);
              double xCenter = (qBoud.at(i) + qBoud.at(i + 1)) / 2.0;
              double yCenter = (qBoud.at(j) + qBoud.at(j + 1)) / 2.0;
              graph2D->SetPoint(pointIndex++, xCenter, yCenter, integral);
          }
      }
  } else if (axes == "out_long") {
      xAxisTitle = "Q_{out} [GeV]";
      yAxisTitle = "Q_{long} [GeV]";
      for (int i = 0; i < qBoud.size() - 1; ++i) {
          for (int j = 0; j < qBoud.size() - 1; ++j) {
              double integral = fit->Integral(qBoud.at(i), qBoud.at(i + 1),
                                      projQmin, projQmax,
                                      qBoud.at(j), qBoud.at(j + 1)) / 
                                      (1.0 * qStep * projQrange * qStep);
              double xCenter = (qBoud.at(i) + qBoud.at(i + 1)) / 2.0;
              double yCenter = (qBoud.at(j) + qBoud.at(j + 1)) / 2.0;
              graph2D->SetPoint(pointIndex++, xCenter, yCenter, integral);
          }
      }
  } else {
      throw std::invalid_argument("Invalid axes for 2D projection");
  }

  graph2D->SetTitle("");
  graph2D->GetXaxis()->SetTitle(xAxisTitle.c_str());
  graph2D->GetYaxis()->SetTitle(yAxisTitle.c_str());
  graph2D->GetZaxis()->SetTitle("Function Value");

  return graph2D;
}

/**
 * \brief Make projection of 3D histogram to one dimension with correct binning.
*/
std::unique_ptr<TH1D> makeProjectionHist(ROOT::Fit::BinData& datastruc, const string& component,
                       double projQmax, double qMax, std::vector<double> qMin3D) {
  bool projOut = false;
  bool projSide = false;
  bool projLong = false;
  unsigned int n = 0;
  int maxn = 100000;
  double projQmaxGeV = projQmax / 1000;
  double val = 1;
  double invErr = 1;
  double value = 1;
  double qMin = qMin3D[0];
  // double qMin = 0.0;
  if (component.compare("out") == 0) {
    projOut = true;
  } else if (component.compare("side") == 0) {
    projSide = true;
  } else if (component.compare("long") == 0) {
    projLong = true;
  } else {
    throw "ERROR: Projection not found!";
  }
  int cycle;
  if (projOut)
    cycle = 0;
  if (projSide)
    cycle = 1;
  if (projLong)
    cycle = 2;
  // define edges and number of bins for qMin
  // std::vector<double> qMin3D = {0.0, 0.28, 0.6, 1.2, 4.2};
  // if (hist_number == 3) {
  //   qMin3D = {0.0, 0.6, 1.2, 4.2};
  // } else if (hist_number == 2) {
  //   qMin3D = {0.0, 1.2, 4.2};
  // }
  // for (auto& value : qMin3D) {
  //   value += qMin;
  // }
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
  // if (hist_number == 3) {
  //   qNBins3D = {15, 10, 7, 0};
  // } else if (hist_number == 2) {
  //   qNBins3D = {10, 7, 0};
  // }

  int nBins = 0;
  std::vector<int> binsUpToQmax;
  vector<double> binVal;
  vector<double> binErr;
  vector<double> newEdges;
  binVal.push_back(0);
  binErr.push_back(0);
  // newEdges.push_back(0);
  for (int i = 0; i < qMin3D.size() - 1; i++) {
    int bins = qNBins3D[i] * ((qMin3D[i + 1] - qMin3D[i]) / (qMin3D[i + 1] - qMin));
    nBins += bins;
    double binWidth = (qMin3D[i + 1] - qMin) / qNBins3D[i];
    for (int j = 0; j < bins; j++) {
      newEdges.push_back(qMin3D[i] + j * binWidth);
      std::cout<<"nbins"<<nBins<<"bins "<<bins<<" newEdges "<<newEdges.back()<<std::endl;
      binVal.push_back(0);
      binErr.push_back(0);
    }
    if (projQmaxGeV > qMin3D[0] && projQmaxGeV <= qMin3D[1]) {
      int binsQmax = static_cast<int>((projQmaxGeV - qMin) / binWidth + 0.001);
      if (binsQmax == 0) {
        binsUpToQmax.push_back(1);
      } else {
        binsUpToQmax.push_back(binsQmax);
      }
    } else {
      int binsQmax = static_cast<int>(0.6 / binWidth);
      binsUpToQmax.push_back(binsQmax);
    }
  }
  //add the last edge
  newEdges.push_back(qMin3D[qMin3D.size() - 1]);
  std::cout<<" newEdges "<<newEdges.back()<<std::endl;
  binVal.push_back(0);
  binErr.push_back(0);
  // TH1D* histDataStruc = new TH1D("histDataStruc", "histDataStruc", nBins, newEdges.data());
  std::unique_ptr<TH1D> histDataStruc = std::make_unique<TH1D>("histDataStruc", "histDataStruc", nBins, newEdges.data());
  histDataStruc->Sumw2();

  // make vector of a size of the number of bins
  vector<double> binVec(nBins+1, 0);
  while (n <= datastruc.Size()) {
    bool skip = true;
    const double *point = datastruc.GetPoint(n, val, invErr);
    double x[3];
    for (int i = 0; i < 3; i++)
      x[i] = point[(i + cycle) % 3];
    n++;
    for (int i = 0; i < 3; i++) {
      if (x[i] < 0) {
        skip = false;
        break;
      }
    }
    if (!skip) continue;
    double binwidth = 0.0;
    for (size_t j = 0; j < qMin3D.size() - 1; j++) {
      if (qMin3D[j] <= x[0] && qMin3D[j + 1] > x[0]) {
        binwidth = (qMin3D[j + 1] - qMin) / qNBins3D[j];
        break;
      }
    }
    if (val == 0) {
      continue;
    }
    if ((x[1] + binwidth / 2 > projQmaxGeV +0.001 || x[2] + binwidth / 2 > projQmaxGeV +0.001) &&
      !(round((x[1] - (binwidth / 2 + qMin))*100) == 0 && round((x[2] - (binwidth / 2 + qMin))*100) == 0)) { 
        if (x[1] < projQmaxGeV && x[2] < projQmaxGeV)
          std::cout<<"Skipping point with x[0] "<<x[0]<<" x[1] "<<x[1]<<" x[2] "<<x[2]<<" val "<<val<<" binwidth / 2 "<<binwidth / 2<<std::endl;
      continue;
    }
    // int norm_val = 1;
    // for (size_t i = 0; i < qMin3D.size() - 1; i++) {
    //   if ((x[0] > qMin3D[i]) && (x[0] <= qMin3D[i + 1])) {
    //     norm_val = binsUpToQmax[i];
    //     break;
    //   }
    // }
    // std::cout<<"x0 "<<x[0]<<" val "<<val<<" norm_val "<<norm_val<<"normalized val "<<val/(norm_val*norm_val)<<" invErr "<<invErr<<" newErr "<<1.0/invErr/invErr/(norm_val*norm_val)<<std::endl;
    double newval = val; // / (norm_val * norm_val);
    double newerr = 1.0 / (invErr * invErr); // / (norm_val * norm_val*norm_val * norm_val);
    binVal.at(histDataStruc->FindBin(x[0]) - 1) += newval;
    binErr.at(histDataStruc->FindBin(x[0]) - 1) += newerr;
    // std::cout<<"binVal "<<binVal.at(histDataStruc->FindBin(x[0]) - 1)<<" binErr "<<binErr.at(histDataStruc->FindBin(x[0]) - 1)<<" bin "<<histDataStruc->FindBin(x[0]) - 1<<" x0 "<<x[0]<<std::endl;
    binVec.at(histDataStruc->FindBin(x[0]) - 1) ++;
  }

  for (int i = 0; i < histDataStruc->GetNbinsX(); i++) {
    if (binVec.at(i) > 0) {
      binVal.at(i) /= binVec.at(i);
      binErr.at(i) = sqrt(binErr.at(i)) / binVec.at(i);
    }
    std::cout<<"binVal "<<binVal.at(i)<<" binErr "<<binErr.at(i)<<" binnorm "<<binVec.at(i)<<std::endl;
    histDataStruc->SetBinContent(i+1, binVal.at(i));
    histDataStruc->SetBinError(i+1, binErr.at(i));
  }

  return histDataStruc;
}

/**
 *  \brief Increases uncertainties for Q values outside of the peak.
 */
void enlargeUncertainties(TH3D* hist) {
  double q_osl;
  double enlargement;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    for (int j = 1; j <= hist->GetNbinsY(); ++j) {
      for (int k = 1; k <= hist->GetNbinsZ(); ++k) {
        q_osl = sqrt(pow(hist->GetXaxis()->GetBinCenter(i), 2) +
        pow(hist->GetYaxis()->GetBinCenter(j), 2) +
        pow(hist->GetZaxis()->GetBinCenter(k), 2))/sqrt(3.);
        if (q_osl > 300.) {
          //enlargement = q_osl / 900. + 5 / 9.; //1 <error< 5x
          enlargement = 10;
          hist->SetBinError(hist->GetBin(i, j, k),
                            enlargement*hist->GetBinErrorLow(i, j, k));
        }
      }
    }
  }
}

/**
 * \brief Scale all histogram axes from MeVs to more readable GeVs.
 */
void scaleAxesMeVtoGeV(TH2D* hist) {
  if (!hist) {
    throw "ERROR: Histogram does no exist!";
  }

  scaleAxisMeVtoGeV(hist, "X");
  scaleAxisMeVtoGeV(hist, "Y");
}

/**
 * \brief Scale axis from MeVs to more readable GeVs.
 */
void scaleAxisMeVtoGeV(TH2D* hist, std::string axisName) {
  TAxis* axis;
  if (axisName.compare("X") == 0) {
    axis = hist->GetXaxis();
  } else if (axisName.compare("Y") == 0) {
    axis = hist->GetYaxis();
  } else if (axisName.compare("Z") == 0) {
    axis = hist->GetZaxis();
  } else {
    throw "ERROR: Axis name not recognized!";
  }

/*
  TArrayD axisArray(*(axis->GetXbins()));

  if (axisArray.GetSize() != 0) {
    for (int i = 0; i < axisArray.GetSize(); ++i) {
      axisArray.SetAt(axisArray.GetAt(i) / 1000., i);
    }
    axis->Set((axisArray.GetSize() - 1), axisArray.GetArray());
  } else {
    axis->SetLimits(axis->GetXmin() / 1000., axis->GetXmax() / 1000.);
  }
  */
}


/**
 * \brief Add fit results to TPaveText.
 */
void addFitResults(TPaveText* text, TF3* func) {
  TString paramText;
  for (int i = 0; i < func->GetNpar(); ++i) {
    if (i == 0 && !showC0) {
      continue;
    }
    if (!useEps && func->GetNpar() == i + 1) {
      continue;
    }

    paramText.Form("%s = %.3f #pm %.3f",
                   func->GetParName(i),
                   func->GetParameter(i),
                   func->GetParError(i));
    text->AddText(paramText.Data());
    paramText.Clear();
  }

  paramText.Form("#chi^{2} / ndf = %.3f ",
                 func->GetChisquare() / func->GetNDF());
  text->AddText(paramText.Data());
  paramText.Clear();
  if (chi_firstBins.first > 0) {
    paramText.Form("#chi^{2}_{%s} / ndf = %.3f",
                   chi_firstBins.second.c_str(),
                   chi_firstBins.first);
    text->AddText(paramText.Data());
    paramText.Clear();
  }
}

/**
 * \brief Output fitted results in csv format
 */
//void outputCSV(TH3D* hist, const std::string& filename) {
void outputCSV(TF3* func, const std::string& filename) {
  //TF3* func = dynamic_cast<TF3*>(hist->GetListOfFunctions()->First());
  if (!func) {
    throw "ERROR: Can't find fitting function!";
  }

  std::ofstream csvOutput;
  csvOutput.open(filename);
  for (int i = 0; i < func->GetNpar(); ++i) {
    if (i == 0 && !showC0) {
      continue;
    }
    if (!useEps && func->GetNpar() == i + 1) {
      continue;
    }
    csvOutput << func->GetParName(i) << ",";
    csvOutput << "#sigma_{" << func->GetParName(i) << "},";
  }
  csvOutput << "#chi^{2}";
  if (chi_firstBins.first > 0) {
    csvOutput << ",#chi^{2}_{" << chi_firstBins.second << "}";
  }
  csvOutput << std::endl;
  for (int i = 0; i < func->GetNpar(); ++i) {
    if (i == 0 && !showC0) {
      continue;
    }
    if (!useEps && func->GetNpar() == i + 1) {
      continue;
    }
    csvOutput << func->GetParameter(i) << ",";
    csvOutput << func->GetParError(i) << ",";
  }
  csvOutput << func->GetChisquare() / func->GetNDF();
  if (chi_firstBins.first > 0) {
    csvOutput << "," << chi_firstBins.first;
  }
  csvOutput << std::endl;
  csvOutput.close();
}

void outputCSV(const FitResultData& fitData, const std::string& filename) {
  std::ofstream csvOutput;
  csvOutput.open(filename);
  for (int i = 0; i < fitData.correlationMatrix.GetNcols(); ++i) {
    csvOutput << "," << fitData.parameterNames[i];
  }
  csvOutput << std::endl;
  for (int i = 0; i < fitData.correlationMatrix.GetNrows(); ++i) {
    csvOutput << fitData.parameterNames[i];
    for (int j = 0; j < fitData.correlationMatrix.GetNcols(); ++j) {
      csvOutput << "," << std::fixed << std::setprecision(4) << fitData.correlationMatrix(i, j);
    }
    csvOutput << std::endl;
    csvOutput << "," << fitData.parameterNames[i];
  }
  csvOutput << std::endl;

  // Print out the eigenvalues and eigenvectors
  for (int i = 0; i < fitData.eigenValues.GetNrows(); ++i) {
    csvOutput << i << ": " << fitData.eigenValues(i);
    for (int j = 0; j < fitData.correlationMatrix.GetNcols(); ++j) {
      csvOutput << "," << fitData.eigenVectors(i, j);
    }
    csvOutput << std::endl;
  }
  csvOutput.close();
}
  

/**
 * \brief Project 3D C2 function onto 2D plane.
 */
void plot2Dprojection(TH3D* inHist,
                      TF3* func,
                      const std::string& axes,
                      std::vector<std::string> commentVec,
                      const TString& outFilePath,
                      int projMin, int projMax) {
  TCanvas* canvas = new TCanvas("canvas", "canvas", 50, 50, 900, 600);
  canvas->SetRightMargin(.35);
  canvas->SetLeftMargin(.11);
  gPad->SetPhi(300);
  gPad->Update();

  if (projMin < 1) {
    projMin = 1;
    std::cout << "WARNING: Projection minimum out of range. "
              << "Using 1 instead!" << std::endl;
  }
  if (projMax > inHist->GetXaxis()->GetNbins() || projMax < 1) {
    projMax = inHist->GetXaxis()->GetNbins();
    std::cout << "WARNING: Projection maximum out of range. "
              << "Using maximal value!" << std::endl;
  }
  string option;
  if (axes.compare("out_side") == 0) {
    option = "yxe";
    inHist->GetZaxis()->SetRange(projMin, projMax);
  }
  if (axes.compare("side_long") == 0) {
    option = "zye";
    inHist->GetXaxis()->SetRange(projMin, projMax);
  }
  if (axes.compare("out_long") == 0) {
    option = "zxe";
    inHist->GetYaxis()->SetRange(projMin, projMax);
  }
  TH2D* hist = dynamic_cast<TH2D*>(inHist->Project3D(option.c_str()));

  inHist->GetXaxis()->UnZoom();
  inHist->GetYaxis()->UnZoom();
  inHist->GetZaxis()->UnZoom();

  std::string histname = inHist->GetName();
  histname += "_";
  histname += axes;
  hist->SetName(histname.c_str());

  hist->Scale(1. / (projMax - projMin + 1));

  // if (zMin > hist->GetMinimum()) {
  //   zMin = hist->GetMinimum() + 0.02;
  // }
  // if (zMax < hist->GetMaximum()) {
  //   zMax = hist->GetMaximum() + 0.02;
  // }
  hist->SetMinimum(zMin);
  hist->SetMaximum(zMax);
  string titleR2 = "";
  if (outFilePath.Contains("c2")) {
    titleR2 = "C_{2}";
  } else {
    titleR2 = "R_{2}";
  }
  hist->SetTitle("");
  if (axes.compare("out_side") == 0) {
    hist->GetXaxis()->SetTitle("Q_{out} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{side} [GeV]");
    hist->GetZaxis()->SetTitle((titleR2 + "(Q_{out}, Q_{side})").c_str());
  }
  if (axes.compare("side_long") == 0) {
    hist->GetXaxis()->SetTitle("Q_{side} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{long} [GeV]");
    hist->GetZaxis()->SetTitle((titleR2 + "(Q_{side}, Q_{long})").c_str());
  }
  if (axes.compare("out_long") == 0) {
    hist->GetXaxis()->SetTitle("Q_{out} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{long} [GeV]");
    hist->GetZaxis()->SetTitle((titleR2 + "(Q_{out}, Q_{long})").c_str());
  }
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);

  scaleAxesMeVtoGeV(hist);

  TPaveText *text = new TPaveText(.68, .5, 1.07, .9, "brNDC");
  text->SetFillStyle(0);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextColor(1);
  text->SetTextFont(42);
  text->SetTextAlign(11);
  addFitResults(text, func);

  TPaveText *text2 = new TPaveText(.68, .15, 1.07, .38, "brNDC");
  text2->SetFillStyle(0);
  text2->SetFillColor(0);
  text2->SetBorderSize(0);
  text2->SetTextColor(1);
  text2->SetTextFont(42);
  text2->SetTextAlign(11);
  for (size_t i = 0; i < commentVec.size(); ++i) {
    text2->AddText(commentVec.at(i).c_str());
  }
  text2->AddText(c2IndexText);

  double projValMin = hist->GetXaxis()->GetBinLowEdge(projMin);
  double projValMax = hist->GetXaxis()->GetBinLowEdge(projMax) +
                      hist->GetXaxis()->GetBinWidth(projMax);

  TString valRangeText;
  if (axes.compare("out_side") == 0) {
    if (projValMax < 1.) {
      valRangeText.Form("%.0f #leq Q_{long} < %.0f MeV",
                        projValMin*1000, projValMax*1000);
    } else {
      valRangeText.Form("%.0f #leq Q_{long} < %.2f GeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }
  if (axes.compare("side_long") == 0) {
    if (projValMax < 1.) {
      valRangeText.Form("%.0f #leq Q_{out} < %.0f MeV",
                        projValMin*1000, projValMax*1000);
    } else {
      valRangeText.Form("%.0f #leq Q_{out} < %.2f GeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }
  if (axes.compare("out_long") == 0) {
    if (projValMax < 1.) {
      valRangeText.Form("%.0f #leq Q_{side} < %.0f MeV",
                        projValMin*1000, projValMax*1000);
    } else {
      valRangeText.Form("%.0f #leq Q_{side} < %.2f GeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }

  // TH2D* histBelow = (TH2D*)hist->Clone("histBelow");
  // TH2D* histAbove = (TH2D*)hist->Clone("histAbove");

  // // Clear the contents of the cloned histograms
  // histBelow->Reset();
  // histAbove->Reset();

  // Iterate over the bins and split the contents
  // for (int i = 1; i <= hist->GetNbinsX(); ++i) {
  //     for (int j = 1; j <= hist->GetNbinsY(); ++j) {
  //         double binContent = hist->GetBinContent(i, j);
  //         double x = hist->GetXaxis()->GetBinCenter(i);
  //         double y = hist->GetYaxis()->GetBinCenter(j);
  //         double fitValue = func->Eval(x, y, 0); // Assuming z=0 for 2D projection

  //         if (binContent < fitValue) {
  //             histBelow->SetBinContent(i, j, binContent);
  //         } else {
  //             histAbove->SetBinContent(i, j, binContent);
  //         }
  //     }
  // }
  // Create the 2D projection of the function
  TGraph2D* funcProjection = make2DProjection(func, axes, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), projValMin, projValMax);
  funcProjection->SetMarkerStyle(20);
  funcProjection->SetMarkerColor(kRed);
  funcProjection->SetFillColor(kRed);
  funcProjection->SetFillColorAlpha(kRed, 0.1); // Set the fill color 
  // Set the axis ranges to match the histogram
  funcProjection->GetXaxis()->SetLimits(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  funcProjection->GetYaxis()->SetLimits(hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
  funcProjection->GetHistogram()->SetMinimum(hist->GetMinimum());
  funcProjection->GetHistogram()->SetMaximum(hist->GetMaximum());
  // std::cout << "Histogram X-axis: " << hist->GetXaxis()->GetXmin() << " to " << hist->GetXaxis()->GetXmax() << std::endl;
  // std::cout << "Histogram Y-axis: " << hist->GetYaxis()->GetXmin() << " to " << hist->GetYaxis()->GetXmax() << std::endl;
  // std::cout << "Histogram Z-axis: " << hist->GetMinimum() << " to " << hist->GetMaximum() << std::endl;
  
  // std::cout << "Graph X-axis: " << funcProjection->GetXaxis()->GetXmin() << " to " << funcProjection->GetXaxis()->GetXmax() << std::endl;
  // std::cout << "Graph Y-axis: " << funcProjection->GetYaxis()->GetXmin() << " to " << funcProjection->GetYaxis()->GetXmax() << std::endl;
  // std::cout << "Graph Z-axis: " << funcProjection->GetZaxis()->GetXmin() << " to " << funcProjection->GetZaxis()->GetXmax() << std::endl;

  // //cout bincontents and coords
  // for (int i = 0; i < funcProjection->GetN(); ++i) {
  //   double x, y, z;
  //   funcProjection->GetPoint(i, x, y, z);
  //   std::cout << "Point " << i << ": (" << x << ", " << y << ", " << z << ")" << std::endl;
  // }
  // //now contents of hist
  // for (int i = 0; i < hist->GetNbinsX(); ++i) {
  //   for (int j = 0; j < hist->GetNbinsY(); ++j) {
  //     double x = hist->GetXaxis()->GetBinCenter(i);
  //     double y = hist->GetYaxis()->GetBinCenter(j);
  //     double z = hist->GetBinContent(i, j);
  //     std::cout << "Histogram Bin (" << i << ", " << j << "): (" << x << ", " << y << ", " << z << ")" << std::endl;
  //   }
  // }

  canvas->SetSupportGL(true);
  hist->Draw("lego2");
  // histBelow->Draw("lego2");
  // histAbove->Draw("lego2 same");
  text->Draw();
  text2->Draw();
  canvas->Print(outFilePath);

  //print funcProjection->Draw("surf1 same"); on second canvas
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 50, 50, 900, 600);
  canvas2->SetRightMargin(.35);
  canvas2->SetLeftMargin(.11);
  gPad->SetPhi(300);
  gPad->Update();
  funcProjection->Draw("surf1");
  text->Draw();
  text2->Draw();
  //add to outFilePath "_func" before .
  std::string outFilePathFunc = outFilePath.Data();
  size_t dotPos = outFilePathFunc.find(".");
  if (dotPos != std::string::npos) {
    outFilePathFunc.insert(dotPos, "_func");
  } else {
    outFilePathFunc += "_func";
  }
  canvas2->Print(outFilePathFunc.c_str());

  delete hist;
  delete funcProjection;
  delete text;
  delete text2;
  delete canvas;
  delete canvas2;
}

void plotCorrelationMatrix(const FitResultData& fitData,
                           const TString& outFilePath) {
  TCanvas* canvas = new TCanvas("canvas", "canvas", 50, 50, 900, 600);
  canvas->SetRightMargin(.35);
  canvas->SetLeftMargin(.11);
  gPad->SetPhi(300);
  gPad->Update();

  TH2D* hist = new TH2D("hist", "hist", fitData.correlationMatrix.GetNcols(),
                        0, fitData.correlationMatrix.GetNcols(),
                        fitData.correlationMatrix.GetNrows(),
                        0, fitData.correlationMatrix.GetNrows());
  hist->SetMinimum(-1);
  hist->SetMaximum(1);
  hist->SetStats(0);
  hist->SetTitle("");
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();

  for (int i = 0; i < fitData.correlationMatrix.GetNrows(); ++i) {
    for (int j = 0; j < fitData.correlationMatrix.GetNcols(); ++j) {
      hist->SetBinContent(i + 1, j + 1, fitData.correlationMatrix(i, j));
    }
    // Remove spaces from the parameter names
    std::string nameWithoutSpaces = fitData.parameterNames[i];
    nameWithoutSpaces.erase(std::remove(nameWithoutSpaces.begin(), nameWithoutSpaces.end(), ' '), nameWithoutSpaces.end());
    // Set the bin labels to the parameter names without spaces
    hist->GetXaxis()->SetBinLabel(i + 1, nameWithoutSpaces.c_str());
    hist->GetYaxis()->SetBinLabel(i + 1, nameWithoutSpaces.c_str());
  }

  hist->Draw("COLZ");  // Draw the histogram

  for (int i = 0; i < fitData.correlationMatrix.GetNrows(); ++i) {
    for (int j = 0; j < fitData.correlationMatrix.GetNcols(); ++j) {
      // Print the bin content in each cell
      TText* text = new TText(j + 0.5, i + 0.5, Form("%.2f", fitData.correlationMatrix(i, j)));
      text->SetTextAlign(22);  // Center alignment
      text->SetTextSize(0.02);  // Set text size
      text->SetNDC(false);  // Use user coordinates instead of normalized coordinates
      text->Draw("SAME");  // Use the "SAME" option to draw on top of the histogram
    }
  }

  canvas->Print(outFilePath);

  delete hist;
  delete canvas;
}
/**
 * \brief Input parser.
 */
InputParser::InputParser (int& argc, char** argv) {
  for (int i = 1; i < argc; ++i)
    this->tokens.push_back(string(argv[i]));
}

const string& InputParser::getCmdOption(const string& option) const {
  vector<string>::const_iterator itr;
  itr = find(this->tokens.begin(), this->tokens.end(), option);
  if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
    return *itr;
  }
  static const string empty_string("");

  return empty_string;
}

bool InputParser::cmdOptionExists(const string& option) const {
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

string MonthNumToName(int month) {
  vector<string> months;
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

string GetCurrentYear() {
  tm now = GetCurrentTime();
  int year = now.tm_year + 1900;

  return std::to_string(year);
}

string GetCurrentMonth() {
  tm now = GetCurrentTime();
  int month = now.tm_mon;

  return MonthNumToName(month);
}

string GetCurrentDay() {
  tm now = GetCurrentTime();
  int day = now.tm_mday;

  string sDay;
  if (day < 10) {
    sDay += "0";
  }
  sDay += std::to_string(day);

  return sDay;
}

string GetCurrentHour() {
  tm now = GetCurrentTime();
  int hour = now.tm_hour;

  string sHour;
  if (hour < 10) {
    sHour += "0";
  }
  sHour += std::to_string(hour);

  return sHour;
}

string GetCurrentMinute() {
  tm now = GetCurrentTime();
  int minute = now.tm_min;

  string sMinute;
  if (minute < 10) {
    sMinute += "0";
  }
  sMinute += std::to_string(minute);

  return sMinute;
}

string GetCurrentSecond() {
  tm now = GetCurrentTime();
  int second = now.tm_sec;

  string sSecond;
  if (second < 10) {
    sSecond += "0";
  }
  sSecond += std::to_string(second);

  return sSecond;
}

string Now() {
  string now = GetCurrentHour();
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

string Today() {
  string today = GetCurrentYear();
  today += "_";
  today += GetCurrentMonth();
  today += "_";
  today += GetCurrentDay();

  return today;
}

