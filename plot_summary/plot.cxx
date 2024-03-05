// std
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
// ROOT
#include <TGraphAsymmErrors.h>
#include <TArrow.h>
// Plot
#include "plot.h"
#include "utils.h"

//std
using std::cout;
using std::endl;


int main(int argc, char** argv) {
  help = false;
  verbose = false;
  plotName = "";
  std::string titles = "";
  std::string comments = "";
  max_binStop = 150;
  c2_type = "";
  yMin = 0.2;
  yMax = 3.0;
  xLabel = "";
  yLabel = "";
  logY = false;
  plotMult = false;
  plotKt = false;

  InputParser input(argc, argv);
  if (input.cmdOptionExists("-h") ||
      input.cmdOptionExists("--help")) {
    help = true;
  }
  if (input.cmdOptionExists("-v") ||
      input.cmdOptionExists("--verbose")) {
    verbose = true;
  }
  if (input.cmdOptionExists("--c0")) {
    paramNameVec.emplace_back("C_{0}");
    paramNameSlugVec.emplace_back("c0");
  }
  if (input.cmdOptionExists("--lambda")) {
    paramNameVec.emplace_back("#lambda");
    paramNameSlugVec.emplace_back("lambda");
  }
  if (input.cmdOptionExists("--r")) {
    paramNameVec.emplace_back("R");
    paramNameSlugVec.emplace_back("r");
  }
  if (input.cmdOptionExists("--r-out")) {
    paramNameVec.emplace_back("R_{out}");
    paramNameSlugVec.emplace_back("r_out");
  }
  if (input.cmdOptionExists("--r-side")) {
    paramNameVec.emplace_back("R_{side}");
    paramNameSlugVec.emplace_back("r_side");
  }
  if (input.cmdOptionExists("--r-long")) {
    paramNameVec.emplace_back("R_{long}");
    paramNameSlugVec.emplace_back("r_long");
  }
  if (input.cmdOptionExists("--alpha-out")) {
    paramNameVec.emplace_back("#alpha_{out}");
    paramNameSlugVec.emplace_back("alpha_out");
  }
  if (input.cmdOptionExists("--alpha-side")) {
    paramNameVec.emplace_back("#alpha_{side}");
    paramNameSlugVec.emplace_back("alpha_side");
  }
  if (input.cmdOptionExists("--alpha-long")) {
    paramNameVec.emplace_back("#alpha_{long}");
    paramNameSlugVec.emplace_back("alpha_long");
  }
  if (input.cmdOptionExists("--eps")) {
    paramNameVec.emplace_back("#epsilon");
    paramNameSlugVec.emplace_back("eps");
  }
  if (input.cmdOptionExists("--chi2")) {
    paramNameVec.emplace_back("#chi^{2}");
    paramNameSlugVec.emplace_back("chi2");
  }
  if (input.cmdOptionExists("-t")) {
    titles = input.getCmdOption("-t");
  }
  if (input.cmdOptionExists("-c")) {
    comments = input.getCmdOption("-c");
  }
  if (input.cmdOptionExists("-o")) {
    plotName = input.getCmdOption("-o");
  }
  if (input.cmdOptionExists("--y-min")) {
    yMin = std::stod(input.getCmdOption("--y-min"));
  }
  if (input.cmdOptionExists("--y-max")) {
    yMax = std::stod(input.getCmdOption("--y-max"));
  }
  if (input.cmdOptionExists("--x-label")) {
    xLabel = input.getCmdOption("--x-label");
  }
  if (input.cmdOptionExists("--y-label")) {
    yLabel = input.getCmdOption("--y-label");
  }
  if (input.cmdOptionExists("--leg-pos")) {
    std::string legendPosition = input.getCmdOption("--leg-pos");
    std::replace(legendPosition.begin(), legendPosition.end(), ',', ' ');
    std::istringstream iss(legendPosition);
    std::vector<double> positions((std::istream_iterator<double>(iss)), std::istream_iterator<double>());
    if (positions.size() == 4) {
        legX[0] = positions[0];
        legX[1] = positions[1];
        legY[0] = positions[2];
        legY[1] = positions[3];
    } else {
        std::cerr << "Invalid legend position. Expected format: x1,x2,y1,y2" << std::endl;
    }
  }
  if (input.cmdOptionExists("--max-binStop")) {
    max_binStop = std::stod(input.getCmdOption("--max-binStop"));
  }
  if (input.cmdOptionExists("--logY")) {
    logY = true;
  }
  if (input.cmdOptionExists("--mult")) {
    plotMult = true;
  }
  if (input.cmdOptionExists("--kt")) {
    plotKt = true;
  }

  if (paramNameVec.empty()) {
    paramNameVec.emplace_back("#lambda");
    paramNameSlugVec.emplace_back("lambda");
  }
  titleVec = splitString(titles, ';');
  commentVec = splitString(comments, ';');

  std::string param;
  bool fileFlag = false;
  for (int i = 0; i < argc; ++i) {
    param = argv[i];
    if (param.compare("-f") == 0) {
      fileFlag = true;
      continue;
    }
    if (fileFlag && param.find(".csv") != std::string::npos) {
      fileVec.emplace_back(param);
    }
  }
  if (fileVec.empty()) {
    help = true;
  }

  if (help) {
    cout << "INFO: Usage: " << endl;
    cout << "      ./plot -[hvto] --[c0 lambda ...] -f FILE_LIST" << endl;
    cout << "             -h                          show help" << endl;
    cout << "             -v                          verbose output" << endl;
    cout << "             --c0                        plot C0" << endl;
    cout << "             --lambda                    plot lambda (p) [default]"
         << endl;
    cout << "             --r                         plot R" << endl;
    cout << "             --r-out                     plot Rout" << endl;
    cout << "             --r-side                    plot Rside" << endl;
    cout << "             --r-long                    plot Rlong" << endl;
    cout << "             --alpha-out                 plot alpha_out" << endl;
    cout << "             --alpha-side                plot alpha_side" << endl;
    cout << "             --alpha-long                plot alpha_long" << endl;
    cout << "             --eps                       plot epsilon" << endl;
    cout << "             --chi2                      plot chi2" << endl;
    cout << "             -t \"TITLE1;TITLE2;...\"      list of titles" << endl;
    cout << "             -c \"COMMENT1;COMMENT2;...\"  list of comments"
         << endl;
    cout << "             -o PLOT_NAME                plot name" << endl;
    cout << "             --y-min VALUE               y-axis minimum" << endl;
    cout << "             --y-max VALUE               y-axis maximum" << endl;
    cout << "             --x-label LABEL             x-axis label" << endl;
    cout << "             --y-label LABEL             y-axis label" << endl;
    cout << "             --leg-pos x1,x2,y1,y2       legend position" << endl;
    cout << "             --max-binStop VALUE         maximum bin_stop" << endl;
    cout << "             --mult                      plot Multiplicity"<< endl;
    cout << "             --kt                        plot kT" << endl;
    cout << "             -f                          list of files" << endl;
    cout << endl;

    return 0;
  }

  try {
    LoadHistograms();
    Plot();
  } catch (const char* msg) {
    cout << msg << endl;

    return 1;
  }

  return 0;
}


void LoadHistograms() {
  cout << "INFO: Loading histograms..." << endl;

  size_t titleIndex = 0;
  if (yLabel.empty()) {
    for (auto paramName : paramNameVec) {
      yLabel += paramName;
      yLabel += ", ";
    }
    yLabel.pop_back();
    yLabel.pop_back();
  }

  std::vector<double> binCenterVecMult{6.09561, 14.4521, 24.6819, 34.7209,
                                       44.7141, 54.6933, 64.6605, 74.6382,
                                       84.6215, 94.56, 109.706, 134.392,
                                       162.282, 211.121};
                                     /*  6.62317, 15.4882, 25.6887, 35.7244,
                                       45.7115, 55.6937, 65.6619, 75.6394,
                                       85.6245, 95.5634, 110.685, 135.382,
                                       163.237, 211.837};*/
 /* 
  std::vector<double> binCenterVecMult{7.00125, 15.4882, 25.6887, 35.7244, 
                                       45.7115, 55.6937, 65.6619, 75.6394,
                                       85.6245, 95.5634, 110.685, 135.382,
                                       163.237, 211.837};//correct data
*/
  std::vector<double> binCenterVecKt{155.117, 255., 355., 455., 553.604,
                                     652.059, 831.163, 1194.78, 1670.};

  for (auto fileName : fileVec) {
    if (!plotMult && !plotKt) {
      if (fileName.find("mult") != std::string::npos) {
        plotMult = true;
        std::cout<< "multiplicity dependancy"<<std::endl;
      } else if (fileName.find("kt") != std::string::npos) {
        plotKt = true;
      } else {
        plotMult = true;
        std::cout<<"Mult nor kt found in file name!";
      }
      if (fileName.find("c2") != std::string::npos) {
        c2_type += "c2_";
        c2_type += fileName.substr(fileName.find("c2") + 3, 1);
      }
    }

    for (auto paramName : paramNameVec) {
      std::ifstream file(fileName);

      csvRow row;
      int paramCol = -1;
      if (file >> row) {
        for (size_t j = 0; j < row.size(); ++j) {
          if (row[j].compare(paramName) == 0) {
            paramCol = j;
          }
        }
      }

      if (paramCol < 0) {
        throw "ERROR: Parameter not found in one of the input files!";
      }
      bool noUncert = false;
      if (paramName.compare("#chi^{2}") == 0) {
        noUncert = true;
      }
      if (static_cast<size_t>(paramCol + 1) >= row.size()) {
        noUncert = true;
      }

      double binStart;
      double binStop;
      double val;
      double valErr;
      int binCenterIndex;
      std::vector<double> binCenterVec;
      if (plotMult) {
        binCenterVec = binCenterVecMult;
      }
      if (plotKt) {
        binCenterVec = binCenterVecKt;
      }
      std::vector<double> binErrVecLeft;
      std::vector<double> binErrVecRight;
      std::vector<double> valVec;
      std::vector<double> valErrVec;
      // check the highest bin_stop
      // while (file >> row) {
      //   if (std::stod(row[1]) > max_binStop)
      //     max_binStop = std::stod(row[1]);
      // }
      // file.clear(); // Clear any error flags
      // file.seekg(0, std::ios::beg); // Move to the beginning of the file
      while (file >> row) {
        binStart = std::stod(row[0]);
        binStop = std::stod(row[1]);
        if (binStop > max_binStop)
          break;
        if (binStop - binStart > 100) {
          break;
        }
        val = std::stod(row[paramCol]);
        if (noUncert) {
          valErr = 0.;
        } else {
          valErr = std::stod(row[paramCol + 1]);
        }

        binCenterIndex = GetVecIndex(binCenterVec, binStart, binStop);
        binErrVecLeft.emplace_back(binCenterVec.at(binCenterIndex) - binStart);
        binErrVecRight.emplace_back(binStop - binCenterVec.at(binCenterIndex));
        valVec.emplace_back(val);
        valErrVec.emplace_back(valErr);
      }

      if (plotMult) {
        auto rightArrow = new TArrow(
              binCenterVec.at(binCenterIndex),
              valVec.back(),
              binCenterVec.at(binCenterIndex) + binErrVecRight.back(),
              valVec.back(),
              0.02, ">");
        arrowVec.emplace_back(rightArrow);

        binCenterVec.emplace_back(binCenterVec.back() + binErrVecRight.back());
        valVec.emplace_back(-1000.);
        binErrVecLeft.emplace_back(0.);
        binErrVecRight.emplace_back(0.);
        valErrVec.emplace_back(0.);

        binErrVecRight.end()[-2] = 0.;
      }
      auto graph = new TGraphAsymmErrors(valVec.size(),
                                         &binCenterVec[0], &valVec[0],
                                         &binErrVecLeft[0], &binErrVecRight[0],
                                         &valErrVec[0], &valErrVec[0]);

      if (xLabel.empty()) {
        if (plotMult) {
          graph->GetXaxis()->SetTitle("N_{ch}");
        }
        if (plotKt) {
          graph->GetXaxis()->SetTitle("k_{T} [MeV]");
        }
      }
      graph->GetYaxis()->SetTitle(yLabel.c_str());
      if (titleIndex < titleVec.size()) {
        graph->SetTitle(titleVec.at(titleIndex).c_str());
      } else {
        graph->SetTitle(paramName.c_str());
      }

      graphVec.emplace_back(graph);
      ++titleIndex;
    }
  }
}

void Plot() {
  cout << "INFO: Plotting..." << endl;

  std::string outFilePath = "output/pdf/" + Today();
  int err = system(("mkdir -p " + outFilePath).c_str());
  if (err != 0) {
    throw "ERROR: Can't create output directory!";
  }

  if (verbose) {
    cout << "INFO: Saving output files here:" << endl;
    cout << "      " << outFilePath << endl;
  }

  if (plotName.empty()) {
    for (auto paramNameSlug : paramNameSlugVec) {
      plotName += "_";
      plotName += paramNameSlug;
    }
    plotName += "_" + c2_type;
  }

  Plotter plotter(outFilePath + "/plot" + plotName);
  for (auto graph : graphVec) {
    plotter.addGraph(graph);
  }
  for (auto arrow : arrowVec) {
    plotter.addArrow(arrow);
  }
  plotter.addNotes(commentVec);
  // plotter.legendY1 = .75;
  // plotter.legendX1 = .1;
  // plotter.legendX2 = .45;
  plotter.legendX1 = legX[0];
  plotter.legendX2 = legX[1];
  plotter.legendY1 = legY[0];
  plotter.legendY2 = legY[1];
  plotter.yMin = yMin;
  plotter.yMax = yMax;
  plotter.logY = logY;
  plotter.draw();
}
