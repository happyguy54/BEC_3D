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
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include <Math/WrappedMultiTF1.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/IFunction.h>
#include <Math/IParamFunction.h>
#include <Math/FitMethodFunction.h>
#include <Fit/FitUtil.h>
#include <Math/FitMethodFunction.h>
#include <Math/Minimizer.h>
#include <Minuit2/FCNBase.h>


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
extern std::pair<double, std::string> chi_firstBins;

using std::vector;
using std::string;
using std::find;
using std::to_string;


// Projections
TGraph* makeProjection(TH3D*,TF3*, const string&, double, double, int, int);
std::unique_ptr<TH1D> makeProjectionHist(ROOT::Fit::BinData&, const string&, double, double, std::vector<double>);
void plot2Dprojection(TH3D*, TF3*, const std::string&,
                      std::vector<std::string>,
                      const TString&,
                      int, int);

// Histogram manipulation
void enlargeUncertainties(TH3D*);
void scaleAxisMeVtoGeV(TH2D*, std::string);
void scaleAxesMeVtoGeV(TH2D*);


class CustomChi2FCN : public ROOT::Math::IMultiGenFunction {
public:
    CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D);
    ~CustomChi2FCN() override = default;

    double operator()(const double* params); // No `override` here
    unsigned int NDim() const override { return fFunc.GetNpar(); }
    double DoEval(const double* params) const override { return const_cast<CustomChi2FCN*>(this)->operator()(params); }
    ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
        return new CustomChi2FCN(fData, fFunc, fQMin3D, fQNBins3D);
    }

private:
    ROOT::Fit::BinData& fData;
    TF3& fFunc;
    const std::vector<double>& fQMin3D;
    const std::vector<int>& fQNBins3D;
    std::unique_ptr<ROOT::Math::WrappedMultiTF1> fWrappedFunc;
    ROOT::Math::IntegratorMultiDim fIntegrator;

    double computeBinWidth(const double* coords) const;
    void setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) const;
};

// class CustomChi2FCN : public ROOT::Math::IMultiGenFunction {
// public:
//     CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D);
//     double operator()(const double* params);

// private:
//     ROOT::Fit::BinData& fData;
//     TF3& fFunc;
//     const std::vector<double>& fQMin3D;
//     const std::vector<int>& fQNBins3D;
//     std::unique_ptr<ROOT::Math::WrappedMultiTF1> fWrappedFunc;
//     ROOT::Math::IntegratorMultiDim fIntegrator;

//     double computeBinWidth(const double* coords);
//     void setBinRange(const double* coords, double binWidth, double* binMin, double* binMax);
// };
// class CustomChi2FCN : public ROOT::Math::IMultiGenFunction {
// public:
//     CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& qMin3D, const std::vector<int>& qNBins3D);
//     double DoEval(const double* x) const override;
//     unsigned int NDim() const override;
//     ROOT::Math::IBaseFunctionMultiDim* Clone() const override;

// private:
//     ROOT::Fit::BinData& fData;
//     TF3& fFunc;
//     const std::vector<double>& fQMin3D;
//     const std::vector<int>& fQNBins3D;
//     std::unique_ptr<ROOT::Math::WrappedMultiTF1> fWrappedFunc;
//     mutable ROOT::Math::IntegratorMultiDim fIntegrator;

//     double computeBinWidth(const double* coords) const;
//     void setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) const;
// };

// class CustomChi2FCN : public ROOT::Math::IBaseFunctionMultiDim {
// public:
//     CustomChi2FCN(ROOT::Fit::BinData& data, TF3& func, const std::vector<double>& binEdges, const std::vector<int>& binCounts);
//     unsigned int NDim() const override;
//     double DoEval(const double* x) const override;
//     CustomChi2FCN* Clone() const override;

// private:
//     ROOT::Fit::BinData& fData;
//     TF3& fFunc;
//     std::vector<double> fBinEdges;
//     std::vector<int> fBinCounts;
//     std::unique_ptr<ROOT::Math::WrappedMultiTF1> fWrappedFunc;
//     mutable ROOT::Math::IntegratorMultiDim fIntegrator;

//     double determineBinWidth(const double* coords) const;
// };
// class CustomChi2FCN : public ROOT::Math::BasicFitMethodFunction<ROOT::Math::IMultiGenFunction> {
// public:
//     typedef ROOT::Math::BasicFitMethodFunction<ROOT::Math::IMultiGenFunction> BaseObjFunction;
//     typedef BaseObjFunction::BaseFunction BaseFunction;
//     typedef ROOT::Math::IParamMultiFunction IModelFunction;
//     typedef BaseObjFunction::Type_t Type_t;

//     CustomChi2FCN(const ROOT::Fit::BinData& data, const IModelFunction& func);
//     virtual ~CustomChi2FCN() {}

//     virtual BaseFunction* Clone() const override;
//     virtual double DoEval(const double* x) const override;
//     virtual double DataElement(const double* x, unsigned int i, double* g = nullptr, double* h = nullptr, bool fullHessian = false) const override;

// private:
//     const ROOT::Fit::BinData& fData;
//     const IModelFunction& fFunc;
//     mutable unsigned int fNEffPoints;
// };


// class CustomChi2FCN : public ROOT::Math::IMultiGenFunction {
// public:
//     CustomChi2FCN(ROOT::Fit::BinData& data, const ROOT::Math::WrappedMultiTF1& wrappedFunc, ROOT::Math::Integrator& integrator)
//         : data_(data), wrappedFunc_(wrappedFunc), fIntegrator_(integrator) {}

//     double operator()(const double* params) const override {
//         double chi2 = 0;
//         for (unsigned int i = 0; i < data_.Size(); ++i) {
//             double value, invError;
//             const double* coords = data_.GetPoint(i, value, invError);

//             double binWidth = computeBinWidth(coords);
//             double binMin[3], binMax[3];
//             setBinRange(coords, binWidth, binMin, binMax);

//             double integral = fIntegrator_.Integral(*wrappedFunc_, binMin, binMax);
//             double binVolume = (binMax[0] - binMin[0]) * (binMax[1] - binMin[1]) * (binMax[2] - binMin[2]);
//             double meanValue = integral / binVolume;

//             double diff = value - meanValue;
//             chi2 += (diff * diff) * (invError * invError);
//         }
//         return chi2;
//     }

//     unsigned int NDim() const override {
//         return data_.Size();
//     }

//     double computeBinWidth(const double* coords) const {
//         double maxCoord = *std::max_element(coords, coords + 3);
//         double binWidth = 0.0;
//         for (size_t k = 0; k < fQMin3D.size() - 1; ++k) {
//             if (fQMin3D[k] <= maxCoord && maxCoord < fQMin3D[k + 1]) {
//                 binWidth = (fQMin3D[k + 1] - fQMin3D[k]) / fQNBins3D[k];
//                 break;
//             }
//         }
//         return binWidth;
//     }

//     void setBinRange(const double* coords, double binWidth, double* binMin, double* binMax) const {
//         for (int j = 0; j < 3; ++j) {
//             binMin[j] = coords[j] - binWidth / 2.0;
//             binMax[j] = coords[j] + binWidth / 2.0;
//         }
//     }

// private:
//     ROOT::Fit::BinData& data_;
//     const ROOT::Math::WrappedMultiTF1& wrappedFunc_;
//     ROOT::Math::Integrator& fIntegrator_;
// };


// Fit results
void addFitResults(TPaveText*, TF3*);
struct FitResultData {
  std::vector<std::string> parameterNames;
  TMatrixDSym correlationMatrix;
  TVectorD eigenValues;
  TMatrixD eigenVectors;
};
//void outputCSV(TH3D*, const std::string&);
void outputCSV(TF3*, const std::string&);
void outputCSV(const FitResultData&, const std::string&);

void plotCorrelationMatrix(const FitResultData&, const TString&); 

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
