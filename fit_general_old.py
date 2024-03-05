#!/usr/bin/env python
# vim: sw=2 ts=8 si nocindent hlsearch

from sys import argv, exit
from os import execlp
from C2 import C2

def use():
  print 'Use:', argv[0], "[-y_min min] [-y_max max] [-h1 hist1] [-h2 hist2] [-mc_h1 mc_hist1] [-mc_h2 mc_hist2] [-title h_title] [-pic name] [-data file] [-data2 file] [-mc file] [-mc2 file] [-rej min] [+rej max] [-rej2 min] [+rej2 max] [-Qmin q] [-Qmax q] [-c2 n] [+-eps] [-mcev mcev] [-fitint]"
  exit(1)

def evalarg(i):
  if len(argv) > i:
    try:
      return eval(argv[i])
    except:
      pass
  use()

def strarg(i):
  if len(argv) > i:
    return argv[i]
  use()

data_hist1 = None
data_hist2 = None
mc_hist1 = None
mc_hist2 = None
title = "++--  MC / MC_opo; Q [MeV]; ++--  MC / MC_opo;"
picname = "ppmm_q_mc_mc_opo"
y_min = 0.80
y_max = 1.6

data_file = "../out/bec_trPt_100_0_Qcut_0.root"
mc_file   = "../out/bec_mc_trPt_100_0_Qcut_0.root"
data_file2 = None
mc_file2   = None
c2_n = 1
use_eps = 1
rej_from = 400.
rej_to = 800.
rej2_from = -1.
rej2_to = -1.
Qmin = 20.
Qmax = -1.
mcev = None
fitpar = "RE"

i = 1
while i < len(argv):
  if argv[i] == '-y_min':
    i += 1
    y_min = evalarg(i)
  elif argv[i] == '-y_max':
    i += 1
    y_max = evalarg(i)
  elif argv[i] == '-mcev':
    i += 1
    mcev = evalarg(i)
  elif argv[i] == '-Qmin':
    i += 1
    Qmin = evalarg(i)
  elif argv[i] == '-Qmax':
    i += 1
    Qmax = evalarg(i)
  elif argv[i] == '-rej':
    i += 1
    rej_from = evalarg(i)
  elif argv[i] == '+rej':
    i += 1
    rej_to = evalarg(i)
  elif argv[i] == '-rej2':
    i += 1
    rej2_from = evalarg(i)
  elif argv[i] == '+rej2':
    i += 1
    rej2_to = evalarg(i)
  elif argv[i] == '-c2':
    i += 1
    c2_n = evalarg(i)
  elif argv[i] == '-eps':
    use_eps = True
  elif argv[i] == '+eps':
    use_eps = False
  elif argv[i] == '-fitint':
    fitpar = "IRE"
  elif argv[i] == '-h1':
    i += 1
    data_hist1 = strarg(i)
  elif argv[i] == '-h2':
    i += 1
    data_hist2 = strarg(i)
  elif argv[i] == '-mc_h1':
    i += 1
    mc_hist1 = strarg(i)
  elif argv[i] == '-mc_h2':
    i += 1
    mc_hist2 = strarg(i)
  elif argv[i] == '-title':
    i += 1
    title = strarg(i)
  elif argv[i] == '-pic':
    i += 1
    picname = strarg(i)
  elif argv[i] == '-data':
    i += 1
    data_file = strarg(i)
  elif argv[i] == '-mc':
    i += 1
    mc_file = strarg(i)
  elif argv[i] == '-data2':
    i += 1
    data_file2 = strarg(i)
  elif argv[i] == '-mc2':
    i += 1
    mc_file2 = strarg(i)
  else:
    use()
  i += 1


c2fun = C2()

rejQ = "rej_from = %s, rej_to = %s" % (rej_from * 1.e-3, rej_to * 1.e-3)
rej2Q = "rej2_from = %s, rej2_to = %s" % (rej2_from * 1.e-3, rej2_to * 1.e-3)


ofile = file("fit.C", "w")
ofile.write("""\
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <iostream>
#include <stdio.h>
#include "AtlasStyle_red.C"
#include "AtlasLabels.C"



Double_t %s;
Double_t %s;
Double_t Qmin = %f;
Double_t Qmax = %f;

""" % (rejQ, rej2Q, Qmin, Qmax))

ofile.write(c2fun.get_fun(c2_n, use_eps))
ofile.write(c2fun.get_fun_red(c2_n, use_eps))

ofile.write("\nTH1D* mk_hist()\n{\n")
if data_hist1:
  hist = "h_data"
  ofile.write("  TFile* DataF = new TFile(\"%s\", \"READ\");\n" % data_file)
  if data_file2:
    ofile.write("  TFile* DataF2 = new TFile(\"%s\", \"READ\");\n" % data_file2)
  ofile.write("  TH1D *h_data;\n")
  if data_hist2:
    ofile.write("  TH1D *h_data1;\n")
    ofile.write("  TH1D *h_data2;\n")
  ofile.write("  DataF -> cd(\"/\");\n")
  if data_hist2:
    ofile.write("  TString dh1, dh2;\n")
    ofile.write("  dh1 = \"%s\";\n" % data_hist1)
    ofile.write("  dh2 = \"%s\";\n" % data_hist2)
    ofile.write("  h_data1 = (TH1D*) DataF -> Get(dh1.Data());\n")
    if data_file2:
      ofile.write("  h_data2 = (TH1D*) DataF2 -> Get(dh2.Data());\n")
    else:
      ofile.write("  h_data2 = (TH1D*) DataF -> Get(dh2.Data());\n")
    ofile.write("  h_data = new TH1D(*h_data1);\n")
    ofile.write("  h_data -> Divide(h_data1, h_data2, 1. / h_data1 -> Integral(), 1. / h_data2 -> Integral());\n")
  else:
    ofile.write("  TString dh;\n")
    ofile.write("  dh = \"%s\";\n" % data_hist1)
    ofile.write("  h_data = (TH1D*) DataF -> Get(dh.Data());\n")

if mc_hist1:
  hist = "h_mc"
  ofile.write("  TFile* McF = new TFile(\"%s\", \"READ\");\n" % mc_file)
  if mc_file2:
    ofile.write("  TFile* McF2 = new TFile(\"%s\", \"READ\");\n" %mc_file2)
  ofile.write("  TH1D *h_mc;\n")
  if mc_hist2:
    ofile.write("  TH1D *h_mc1;\n")
    ofile.write("  TH1D *h_mc2;\n")
  ofile.write("  McF -> cd(\"/\");\n")
  if mc_hist2:
    ofile.write("  TString mh1, mh2;\n")
    ofile.write("  mh1 = \"%s\";\n" % mc_hist1)
    ofile.write("  mh2 = \"%s\";\n" % mc_hist2)
    ofile.write("  h_mc1 = (TH1D*) McF -> Get(mh1.Data());\n")
    if mc_file2:
      ofile.write("  h_mc2 = (TH1D*) McF2 -> Get(mh2.Data());\n")
    else:
      ofile.write("  h_mc2 = (TH1D*) McF -> Get(mh2.Data());\n")
    ofile.write("  h_mc = new TH1D(*h_mc1);\n")
    ofile.write("  h_mc -> Divide(h_mc1, h_mc2, 1. / h_mc1 -> Integral(), 1. / h_mc2 -> Integral());\n")
  else:
    ofile.write("  TString mh;\n")
    ofile.write("  mh = \"%s\";\n" % mc_hist1)
    ofile.write("  h_mc = (TH1D*) McF -> Get(mh.Data());\n")

if data_hist1 and mc_hist1:
  hist = "h_r"
  if mcev:
    if data_hist2:
      ofile.write("  h_data1 -> Add(h_mc1, %f);\n" % (mcev / 100.))
      ofile.write("  h_data2 -> Add(h_mc2, %f);\n" % (mcev / 100.))
      ofile.write("  h_data -> Divide(h_data1, h_data2, 1. / h_data1 -> Integral(), 1. / h_data2 -> Integral());\n")
    else:
      ofile.write("  h_data -> Add(h_mc, %f);\n" % (mcev / 100.))
  ofile.write("  TH1D *h_r;\n")
  ofile.write("  h_r = new TH1D(*h_data);\n")
  if data_hist2 or mc_hist2:
    ofile.write("  h_r -> Divide(h_data, h_mc);\n")
  else:
    ofile.write("  h_r -> Divide(h_data, h_mc, 1. / h_data -> Integral(), 1. / h_mc -> Integral());\n")

ofile.write("""
  %s -> SetTitle("%s");
  %s -> GetXaxis() -> SetNdivisions(504);
  %s -> GetYaxis() -> SetRangeUser(%f, %f);

  return %s;
}
""" % (hist, title, hist, hist, y_min, y_max, hist))

ofile.write("""
void fit()
{

  SetAtlasStyle();
  gStyle -> SetOptFit(11111);

  TCanvas  *c3 = new TCanvas("c3", "c3", 50,50,600,600);
  c3 -> Update();

""")

ofile.write('  TH1D *hst;\n')
ofile.write('  hst = mk_hist();\n')

npars = c2fun.get_n_par(c2_n)
if use_eps:
  npars = npars + 1


ofile.write("""
  if (Qmax < 0.)
    Qmax = hst->GetXaxis()->GetXmax();
  TF1 *fun_c2 = new TF1("%(c2name)s", %(c2name)s, Qmin, Qmax, %(npars)d);
  fun_c2 -> SetParameters(%(def_par)s);
  fun_c2 -> SetParLimits(1, 0., 1.);
""" % {"c2name": c2fun.get_funname(c2_n),
       "npars": npars,
       "def_par": c2fun.get_parameters(c2_n),
       "params": c2fun.get_parnames(c2_n)})

ofile.write("""
  fun_c2 -> SetParNames(%(params)s);
  fun_c2 -> SetLineWidth(1);
  fun_c2 -> SetLineColor(kRed);

  fun_c2 -> FixParameter(1, 0.88);
  hst -> Fit(fun_c2, "%(fitpar)s");
  gMinuit->mnmatu(1);

  cout << "  -1  #chi2    " << fun_c2 -> GetChisquare()/fun_c2 -> GetNDF() << "    0.0" << endl;

  double cov_matrix[4][4];
  gMinuit->mnemat(cov_matrix[0], 4);

  double cov_R_lambda = cov_matrix[1][2];
  double corr_R_lambda = cov_R_lambda / (fun_c2->GetParError(1) * fun_c2->GetParError(2));

  cout << "  -2  #rho    " << corr_R_lambda << "    0.0" << endl;


  hst -> Draw();

""" % {"c2name": c2fun.get_funname(c2_n),
       "npars": npars,
       "fitpar": fitpar,
       "def_par": c2fun.get_parameters(c2_n),
       "params": c2fun.get_parnames(c2_n)})

ofile.write("""
  FILE *out_f = fopen("fit.tmp", "w");
""")
for i in xrange(npars):
  ofile.write("""  fprintf(out_f, "%d,\\"%s\\",%%f,%%f\\n", fun_c2 -> GetParameter(%d), fun_c2 -> GetParError(%d));\n""" % (i,c2fun.get_parname(c2_n, i), i, i))
ofile.write("""  fprintf(out_f, "-1,\\"#chi2\\",%f,%d\\n", fun_c2 -> GetChisquare(), fun_c2 -> GetNDF());\n""")
ofile.write('  if (corr_R_lambda < -10) {\n')
ofile.write('    corr_R_lambda = -10;\n')
ofile.write('  }\n')
ofile.write('  if (corr_R_lambda > 10) {\n')
ofile.write('    corr_R_lambda = 10;\n')
ofile.write('  }\n')
ofile.write('  fprintf(out_f, "-2,\\"#rho\\",%f,%d\\n", corr_R_lambda, 0);\n')


if rej_to > rej_from:
  ofile.write("""
  TBox *exl_box = new TBox(1000. * rej_from, %(min)f, 1000. * rej_to, %(max)f);
  exl_box -> SetFillColor(kBlack);
  exl_box -> SetFillStyle(3954);
  // exl_box -> Draw();
  double ar_y = %(min)f + 0.1 * (%(max)f - %(min)f);
  TArrow *ar1 = new TArrow(1000. * rej_from, ar_y, 1000. * rej_to, ar_y, 0.015, "<|>");
  ar1 -> Draw();

  TLatex lexc;
  lexc.SetTextSize(0.041);
  lexc.SetTextAlign(23);
  ar_y = %(min)f + 0.08 * (%(max)f - %(min)f);
  lexc.DrawLatex(500. * (rej_from + rej_to) - 100., ar_y, "excluded");
  """ % {'min': y_min, 'max': y_max})

if rej2_to > rej2_from:
  ofile.write("""
  TBox *exl2_box = new TBox(1000. * rej2_from, %(min)f, 1000. * rej2_to, %(max)f);
  exl2_box -> SetFillColor(kBlack);
  exl2_box -> SetFillStyle(3954);
  // exl_box -> Draw();
  double ar2_y = %(min)f + 0.1 * (%(max)f - %(min)f);
  TArrow *ar2 = new TArrow(1000. * rej2_from, ar2_y, 1000. * rej2_to, ar2_y, 0.015, "<|>");
  ar2 -> Draw();

  TLatex lexc2;
  lexc2.SetTextSize(0.041);
  lexc2.SetTextAlign(23);
  ar2_y = %(min)f + 0.08 * (%(max)f - %(min)f);
  lexc2.DrawLatex(500. * (rej_from + rej_to) - 100., ar2_y, "excluded");
  """ % {'min': y_min, 'max': y_max})

if Qmin > 0.:
  ofile.write("""
  TBox *und_box = new TBox(0., %f, Qmin, %f);
  und_box -> SetFillColor(kBlack);
  und_box -> SetFillStyle(3954);
  // und_box -> Draw();
""" % (y_min, y_max))

if Qmax > 0.:
  ofile.write("""
  TBox *ovr_box = new TBox(Qmax, %f, hst->GetXaxis()->GetXmax(), %f);
  ovr_box -> SetFillColor(kBlack);
  ovr_box -> SetFillStyle(3954);
  // ovr_box -> Draw();
""" % (y_min, y_max))

ofile.write("""
  ATLASLabel(0.2, 0.88, (char *) "Preliminary");

  TString picfile, stmp;
  stmp = "%s.";
  picfile = stmp;
  c3 -> Print(picfile + "png");
  c3 -> Print(picfile + "eps");
  c3 -> Print(picfile + "pdf");

""" % picname)

ofile.write("""
  TH1D *rez_hst = new TH1D(*hst);
  int i;
  double sumch2 = 0.;
  int Nch = 0.;
  for (i = 1; i <= hst -> GetNbinsX(); i++)
  {
    double Q = hst -> GetBinCenter(i);
    rez_hst -> SetBinContent(i, (fun_c2 -> Eval(Q) - hst -> GetBinContent(i)) / hst -> GetBinError(i));
    rez_hst -> SetBinError(i, 0);
    if (Q > Qmax)
      continue;
    if (Q < Qmin)
      continue;
    if (Q > rej_from * 1000. && Q < rej_to * 1000.)
      continue;
    if (Q > rej2_from * 1000. && Q < rej2_to * 1000.)
      continue;
    Nch++;
    sumch2 += pow((fun_c2 -> Eval(Q) - hst -> GetBinContent(i)) / hst -> GetBinError(i), 2.);
  }
  cout << Nch << " " << sumch2 << "\\n";
  rez_hst -> ResetStats();
  rez_hst -> SetMaximum(1.05 * rez_hst -> GetBinContent(rez_hst -> GetMaximumBin()));
  rez_hst -> SetMinimum(1.05 * rez_hst -> GetBinContent(rez_hst -> GetMinimumBin()));
  rez_hst -> GetYaxis() -> SetTitle("#epsilon");
  rez_hst -> Draw();

  double ar_y2 = rez_hst -> GetMinimum() + 0.1 * (rez_hst -> GetMaximum() - rez_hst -> GetMinimum());
  TArrow *ar2 = new TArrow(1000. * rej_from, ar_y2, 1000. * rej_to, ar_y2, 0.015, "<|>");
  ar2 -> Draw();
  TArrow *ar3 = new TArrow(1000. * rej2_from, ar_y2, 1000. * rej2_to, ar_y2, 0.015, "<|>");
  ar3 -> Draw();  
  TLatex lexc3;
  lexc3.SetTextSize(0.041);
  lexc3.SetTextAlign(23);
  ar_y2 = rez_hst -> GetMinimum() + 0.08 * (rez_hst -> GetMaximum() - rez_hst -> GetMinimum());
  if (rej2_to > rej2_from)
    lexc3.DrawLatex(500. * (rej2_from + rej_to), ar_y2, "excluded");
  else
    lexc3.DrawLatex(500. * (rej_from + rej_to), ar_y2, "excluded");
""")

if rej_to > rej_from:
  ofile.write("""
  TBox *exl_box2 = new TBox(1000. * rej_from, rez_hst -> GetMinimum(), 1000. * rej_to, rez_hst -> GetMaximum());
  exl_box2 -> SetFillColor(kBlack);
  exl_box2 -> SetFillStyle(3954);
  // exl_box2 -> Draw();
""")

if Qmin > 0.:
  ofile.write("""
  TBox *und_box2 = new TBox(0., rez_hst -> GetMinimum(), Qmin, rez_hst -> GetMaximum());
  und_box2 -> SetFillColor(kBlack);
  und_box2 -> SetFillStyle(3954);
  // und_box2 -> Draw();
""")

if Qmax > 0.:
  ofile.write("""
  TBox *ovr_box2 = new TBox(Qmax, rez_hst -> GetMinimum(), rez_hst->GetXaxis()->GetXmax(), rez_hst -> GetMaximum());
  ovr_box2 -> SetFillColor(kBlack);
  ovr_box2 -> SetFillStyle(3954);
  // ovr_box2 -> Draw();
""")

ofile.write("""
  TLine l(0., 0., hst->GetXaxis()->GetXmax(), 0.);
  l.Draw();
  TPaveText *pt = new TPaveText(0.7,0.35,0.98,0.48,"brNDC");
  TString st_ch2, st_nch;
  st_ch2.Form("sum #epsilon^{2} = %f", sumch2);
  st_nch.Form("N points = %d", Nch);
  pt -> AddText(st_ch2.Data());
  pt -> AddText(st_nch.Data());
  // pt -> Draw();
  c3 -> Print(picfile + "chi2." + "pdf");
  c3 -> Print(picfile + "chi2." + "png");
}
""")

ofile.close()

execlp("root", 'root', '-b', '-q', '-x', 'fit.C')
