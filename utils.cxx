// std
#include <iostream>
#include <fstream>
// ROOT
#include <TGraph.h>
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

/**
 * \brief Make projection of 3D histogram to one dimension with correct binning.
*/
TH1D* makeProjectionHist(ROOT::Fit::BinData& datastruc, const string& component,
                       double projQmax, int projRange, int hist_number) {
  bool projOut = false;
  bool projSide = false;
  bool projLong = false;
  unsigned int n=0;
  int maxn = 100000;
  double val = 1;
  double invErr = 1;
  double value = 1;
  double qMin = 0.02;
  if (component.compare("out") == 0) {
    projOut = true;
  } else if (component.compare("side") == 0) {
    projSide = true;
  } else if (component.compare("long") == 0) {
    projLong = true;
  } else {
    throw "ERROR: Projection not found!";
  }
  int cycle; //if "out" then cycle = 0
  if (projOut)
    cycle = 0;
  if (projSide)
    cycle = 1;
  if (projLong)
    cycle = 2;
  // define edges and number of bins for qMin
  std::vector<double> qMin3D = {0.0, 0.28, 0.6, 1.2, 4.2};
  // all components of qMin3D +=0.02
  for(auto& value : qMin3D) {
    value += qMin;
  }
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
      //std::cout<<qMin3D[i] + j*binWidth;
      binVal.push_back(0);
      binErr.push_back(0);
    }
  }
  TH1D* histDataStruc = new TH1D("histDataStruc", "histDataStruc", nBins, newEdges.data());
  histDataStruc->Sumw2();
    while (n<=datastruc.Size()) {
      bool skip = true;
      const double *point = datastruc.GetPoint(n, val, invErr);
      double x[3];
      for (int i = 0; i < 3; i++)
      x[i] = point[(i+cycle)%3];
      //x[2] = temp;
      n++;
      for (int i = 0; i < 3; i++) {
        if (x[i] < 0) {
            skip = false;
            break;
        }
      } 
      if (!skip) continue;
      if ((((x[1]>projQmax/1000)||(x[2]>projQmax/1000))&&(x[0]<0.6+qMin))||(val==0)) {
        /*std::cout<<"skipping point N= "<<n<<"value= "<<val;
        for (int i = 0; i < 3; i++) {
          std::cout << ", x[" << i << "] = " << x[i];
        }
        std::cout<<std::endl;*/
        continue;
      }
      if (((x[1]>0.6+qMin)||(x[2]>0.6+qMin))&&(x[0]>=0.6+qMin))
        continue;
      //if ((x[1]>projQmax)&&(x[2]>projQmax)) break;
      /*std::cout<<" N= "<<n<<"value= "<<val;
      for (int i = 0; i < 3; i++) {
        std::cout << ", x[" << i << "] = " << x[i];
      }
      std::cout<<std::endl;*/
      binVal.at(histDataStruc->FindBin(x[0])-1) += val;
      binErr.at(histDataStruc->FindBin(x[0])-1) += 1.0/(invErr*invErr);
      //std::cout<<"Binerror"<<binErr.at(histDataStruc->FindBin(x[0])-1)<<std::endl;
    //histDataStruc->Fill(x[0],val/(projRan * projRan));
    }
    
    //histDataStruc->Scale(1. / (projRange * projRange));
    //std::cout<<"N = "<<n<<std::endl;
    for (int i = 1; i<histDataStruc->GetNbinsX();i++) {
      for (int j = 0; j<qMin3D.size(); j++) {
        if ((qMin3D[j]<=histDataStruc->GetBinCenter(i))&&(qMin3D[j+1]>histDataStruc->GetBinCenter(i)))
          projRan = j;
      }
      if (hist_number == 3 && projRan < 1) {
        projRan = 1;
      } else if (hist_number == 2 && projRan < 2) {
        projRan = 2;
      }
      switch(projRan){
        case 0: projRan = projRange/(projRan+1); //projRange/(100*qMin3D[i+1]/qNBins3D[i])*2;
          break;
        case 1: projRan = floor(projRange/(projRan+1));
          break;
        case 2: projRan = 5;
          break;
        case 3: projRan = 1;
          break;
      }
      // std::cout<<"projRange "<<projRange<<" projRan "<<projRan<<std::endl;
      histDataStruc->SetBinContent(i, binVal.at(i-1)/(projRan * projRan));
      histDataStruc->SetBinError(i, sqrt(binErr.at(i-1))/(projRan * projRan));
      // std::cout<<"Bin "<<i<<", xaxis "<<histDataStruc->GetXaxis()->GetBinCenter(i)<<", value: "<<histDataStruc->GetBinContent(i)<<", error: "<<histDataStruc->GetBinError(i)<<std::endl;
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
  csvOutput << std::endl;
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

  hist->SetMinimum(zMin);
  hist->SetMaximum(zMax);

  hist->SetTitle("");
  if (axes.compare("out_side") == 0) {
    hist->GetXaxis()->SetTitle("Q_{out} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{side} [GeV]");
    hist->GetZaxis()->SetTitle("R_{2}(Q_{out}, Q_{side})");
  }
  if (axes.compare("side_long") == 0) {
    hist->GetXaxis()->SetTitle("Q_{side} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{long} [GeV]");
    hist->GetZaxis()->SetTitle("R_{2}(Q_{side}, Q_{long})");
  }
  if (axes.compare("out_long") == 0) {
    hist->GetXaxis()->SetTitle("Q_{out} [GeV]");
    hist->GetYaxis()->SetTitle("Q_{long} [GeV]");
    hist->GetZaxis()->SetTitle("R_{2}(Q_{out}, Q_{long})");
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
      valRangeText.Form("%.0f #leq Q_{long} < %.0f GeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }
  if (axes.compare("side_long") == 0) {
    if (projValMax < 1.) {
      valRangeText.Form("%.0f #leq Q_{out} < %.0f MeV",
                        projValMin*1000, projValMax*1000);
    } else {
      valRangeText.Form("%.0f #leq Q_{out} < %.0f GeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }
  if (axes.compare("out_long") == 0) {
    if (projValMax < 1.) {
      valRangeText.Form("%.0f #leq Q_{side} < %.0f GeV",
                        projValMin*1000, projValMax*1000);
    } else {
      valRangeText.Form("%.0f #leq Q_{side} < %.0f MeV",
                        projValMin, projValMax);
    }
    text2->AddText(valRangeText.Data());
  }

  hist->Draw("SURF1");
  text->Draw();
  text2->Draw();
  canvas->Print(outFilePath);

  delete hist;
  delete text;
  delete text2;
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

