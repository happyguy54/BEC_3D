// std
#include<iostream>
#include<cmath>
#include<algorithm>
// ROOT
#include<TFile.h>
#include<TF3.h>
// Fit
#include "c2.h"


using std::cout;
using std::endl;
using std::sqrt;
using std::exp;


C2::C2() {
  index = 0;
  nParams = 0;
  rejFrom = 1.;
  rejTo = -1.;
  rejFromOsl = {1., 1., 1.};
  rejToOsl = {-1., -1., -1.};
  rej2From = 1;
  rej2To = -1;
  useEps = false;
  fixC0 = false;
}

void C2::setup(TF3* func) {
  func->SetParName(0, "C_{0}");

  func->SetParName(2, "R_{out}");
  func->SetParName(3, "R_{side}");
  func->SetParName(4, "R_{long}");

  func->SetParameter(0, 1.);
  func->SetParameter(1, .5);
  func->SetParameter(2, 1.2);
  func->SetParameter(3, 1.2);
  func->SetParameter(4, 1.2);

  func->SetParLimits(0, 0., 100.);
  func->SetParLimits(1, 0., 1.);
  func->SetParLimits(2, 0., 100.);
  func->SetParLimits(3, 0., 100.);
  func->SetParLimits(4, 0., 100.);
}


C2_1::C2_1() : C2::C2() {
  index = 1;
  nParams = 6;
}

void C2_1::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#epsilon");
  func->SetParameter(5, 1e-3);
  func->SetParLimits(5, -100., 100.);
  if (!useEps) {
    func->FixParameter(5, 1e-3);
  }
  if (fixC0) {
    func->FixParameter(0, 1.);
  }
}

double C2_1::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double eps = 1.;
  if (useEps) {
      eps = 1. + (q_osl * p[5]);
  }

  double fitval = c0 * (1 + lambda * exp(- pow(r_out * q_out, 2)
                                         - pow(r_side * q_side, 2)
                                         - pow(r_long * q_long, 2))) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }

  return fitval;
}


C2_2::C2_2() : C2::C2() {
  index = 2;
  nParams = 6;
}

void C2_2::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#epsilon");
  func->SetParameter(5, 1e-3);
  func->SetParLimits(5, -100., 100.);
  if (!useEps) {
    func->FixParameter(5, 1e-3);
  }
  if (fixC0) {
    func->FixParameter(0, 1.);
  }
}

double C2_2::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double eps = 1.;
  if (useEps) {
      eps = 1. + (q_osl * p[5]);
  }

  double fitval = c0 * (1 + lambda * exp(- r_out * q_out
                                         - r_side * q_side
                                         - r_long * q_long)) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }

  return fitval;
}


C2_3::C2_3() : C2::C2() {
  index = 3;
  nParams = 7;
}

void C2_3::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#alpha");
  func->SetParameter(5, .5);
  func->SetParLimits(5, 0., 5.);

  func->SetParName(6, "#epsilon");
  func->SetParameter(6, 1e-3);
  func->SetParLimits(6, -100., 100.);
  if (!useEps) {
    func->FixParameter(6, 1e-3);
  }
}

double C2_3::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double alpha = p[5];
  double eps = 1.;
  if (useEps) {
    eps = 1. + (q_osl * p[6]);
  }

  double fitval = c0 * (1 + lambda * exp(- pow(r_out * q_out, alpha)
                                         - pow(r_side * q_side, alpha)
                                         - pow(r_long * q_long, alpha))) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }

  return fitval;
}

C2_4::C2_4() : C2::C2() {
  index = 4;
  nParams = 9;
}

void C2_4::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#alpha_{out}");
  func->SetParName(6, "#alpha_{side}");
  func->SetParName(7, "#alpha_{long}");
  func->SetParameter(5, 1.5);
  func->SetParameter(6, 1.5);
  func->SetParameter(7, 1.5);
  func->SetParLimits(5, 0., 5.);
  func->SetParLimits(6, 0., 5.);
  func->SetParLimits(7, 0., 5.);
  
  func->SetParName(8, "#epsilon");
  func->SetParameter(8, 1e-3);
  func->SetParLimits(8, -100., 100.);
  if (!useEps) {
    func->FixParameter(8, 1e-3);
  }
}

double C2_4::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double alpha_out = p[5];
  double alpha_side = p[6];
  double alpha_long = p[7];
  double eps = 1.;
  if (useEps) {
    eps = 1. + (q_osl * p[8]);
  }

  double fitval = c0 * (1 + lambda * exp(- pow(r_out * q_out, alpha_out)
                                         - pow(r_side * q_side, alpha_side)
                                         - pow(r_long * q_long, alpha_long))) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }

  return fitval;
}

C2_5::C2_5() : C2::C2() {
  index = 5;
  nParams = 7;
}

void C2_5::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#alpha_{out}");
  func->SetParameter(5, 1.5);
  func->SetParLimits(5, 0., 5.);
  
  func->SetParName(6, "#epsilon");
  func->SetParameter(6, 1e-3);
  func->SetParLimits(6, -100., 100.);
  if (!useEps) {
    func->FixParameter(6, 1e-3);
  }
}

double C2_5::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double alpha_out = p[5];
  double alpha_side = 2;
  double alpha_long = 1;
  double eps = 1.;
  if (useEps) {
    eps = 1. + (q_osl * p[6]);
  }

  double fitval = c0 * (1 + lambda * exp(- pow(r_out * q_out, alpha_out)
                                         - pow(r_side * q_side, alpha_side)
                                         - pow(r_long * q_long, alpha_long))) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }
  std::vector<int> diff(rejToOsl.size());
  std::transform(rejToOsl.begin(), rejToOsl.end(), rejFromOsl.begin(), diff.begin(), std::minus<int>());
  bool biggerOsl = std::all_of(diff.begin(), diff.end(), [](int i){ return i > 0; });

  if (biggerOsl) {
    bool rejOsl = false;
    for (size_t i = 0; i < rejFromOsl.size(); ++i) {
      if ((x[i] > rejFromOsl[i] && x[i] < rejToOsl[i])) {
        rejOsl = true;
        break;
      }
    }
    if (rejOsl) {
        TF3::RejectPoint();
    }
  }

  return fitval;
}

C2_6::C2_6() : C2::C2() {
  index = 6;
  nParams = 8;
}

void C2_6::setup(TF3* func) {
  C2::setup(func);

  func->SetParName(1, "#lambda");

  func->SetParName(5, "#alpha_{out}");
  func->SetParName(6, "#alpha_{side}");
  func->SetParName(7, "#alpha_{long}");
  func->SetParameter(5, 1.5);
  func->SetParameter(6, 1.5);
  func->SetParameter(7, 1.5);
  func->SetParLimits(5, 0., 5.);
  func->SetParLimits(6, 0., 5.);
  func->SetParLimits(7, 0., 5.);
  
  // func->SetParName(8, "#epsilon");
  // func->SetParameter(8, 1e-3);
  // func->SetParLimits(8, -10., 10.);
  // if (!useEps) {
  //   func->FixParameter(8, 1e-3);
  // }
  // func->FixParameter(0, 1);
  // func->FixParameter(7, 1);
}

double C2_6::operator() (double* x, double* p) {
  double q_out = x[0];
  double q_side = x[1];
  double q_long = x[2];
  double q_osl = sqrt(q_out * q_out +
                      q_side * q_side +
                      q_long * q_long) / sqrt(3.);
  double c0 = p[0];
  double lambda = p[1];
  double r_out = p[2] / 0.197;
  double r_side = p[3] / 0.197;
  double r_long = p[4] / 0.197;
  double alpha_out = p[5];
  double alpha_side = p[6];
  double alpha_long = p[7];
  double eps = 1.;
  // if (useEps) {
  //   eps = 1. + (q_osl * p[8]);
  // }

  double fitval = c0 * (1 + lambda * exp(- pow(r_out * q_out, alpha_out)
                                         - pow(r_side * q_side, alpha_side)
                                         - pow(r_long * q_long, alpha_long))) * eps;

  q_osl = q_osl;
  if (rejFrom < rejTo && q_osl > rejFrom && q_osl < rejTo) {
      TF3::RejectPoint();
  }
  if (rej2From < rej2To && q_osl > rej2From && q_osl < rej2To) {
      TF3::RejectPoint();
  }
  std::vector<int> diff(rejToOsl.size());
  std::transform(rejToOsl.begin(), rejToOsl.end(), rejFromOsl.begin(), diff.begin(), std::minus<int>());
  bool biggerOsl = std::all_of(diff.begin(), diff.end(), [](int i){ return i > 0; });
  if (biggerOsl) {
    bool rejOsl = false;
    // std::cout << "rejFromOsl.size(): " << rejFromOsl.size() << endl;
    for (size_t i = 0; i < rejFromOsl.size(); ++i) {
      if ((x[i] > rejFromOsl[i] && x[i] < rejToOsl[i])) {
        rejOsl = true;
        break;
      }
    }
    std::cout << "before TF3::RejectPoint(); rejOsl: " << rejOsl << endl;
    if (rejOsl) {
        TF3::RejectPoint();
    }
  }

  return fitval;
}

//par[0] * (1 + 2 * par[1] * (1-par[1]) * exp(-pow(rR * X, 2)) + pow(par[1], 2) * exp(-2 * pow(rR * X, 2)))