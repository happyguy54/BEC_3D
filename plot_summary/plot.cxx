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
  max_binStop = 250;//150;
  // max_binStart = 125;
  c2_type = "";
  c2_number = "";
  ref_type = "";
  ref_num = 1;
  yMin = 0.2;
  yMax = 3.0;
  xLabel = "";
  yLabel = "";
  textVec = {};
  logY = false;
  plotMult = false;
  plotKt = false;
  plottt = false;

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
  if (input.cmdOptionExists("--tt")) {
    plottt = true;
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
  std::vector<std::string> qMaxVec;
  if (yLabel.empty()) {
    for (auto paramName : paramNameVec) {
      yLabel += paramName;
      yLabel += ", ";
    }
    yLabel.pop_back();
    yLabel.pop_back();
  }

  std::vector<double> lastgraph{};
  std::vector<double> filteredBinCenterVec{};
  std::vector<double> binCenterVecMult{6.09561, 14.4521, 24.6819, 34.7209,
                                       44.7141, 54.6933, 64.6605, 74.6382,
                                       84.6215, 94.56, 109.706, 134.392,
                                       162.282, 211.121, //older values
                                       43.48};//, 21.78}; //values added from nonunfolded multiplicity
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
                                     652.059, 831.163, 1194.78, 1670., //older values
                                     482.7}; //values added from nonunfolded multiplicity (does not matter, whole interval used)
  std::vector<double> binEdgesVectt{0.63, 0.71, 0.78, 0.85, 0.93, 1};
  std::vector<double> binCenterVectt;
  for (size_t i = 0; i < binEdgesVectt.size() - 1; ++i) {
    binCenterVectt.push_back((binEdgesVectt[i] + binEdgesVectt[i + 1]) / 2);
  }

  for (auto fileName : fileVec) {
    std::string qMax = "";
    std::string rejected = "";
    cout<<"File: "<<fileName<<endl;
    if (!plotMult && !plotKt && !plottt) {
      if (fileName.find("mult") != std::string::npos) {
        plotMult = true;
        std::cout<< "multiplicity dependancy"<<std::endl;
      } else if (fileName.find("kt") != std::string::npos) {
        plotKt = true;
        std::cout<< "kt dependancy"<<std::endl;
      } else if (fileName.find("tt") != std::string::npos) {
        plottt = true;
        std::cout<< "tt dependancy"<<std::endl;
      } else {
        plotMult = true;
        std::cout<<"Mult nor kt found in file name!";
      }
    }
    if (fileName.find("c2") != std::string::npos) {
      if (c2_type.empty()) {
        c2_type += "c2_";
        // c2_type += fileName.substr(fileName.find("c2") + 3, 1);
        size_t start = fileName.find("c2_") + 3;
        size_t end = fileName.find_first_of("_/", start);
        if (end != std::string::npos) {
            c2_type += fileName.substr(start, end - start);
            c2_number = fileName.substr(start, end - start);
        }
      } else {
        c2_type += "vs";
        size_t start = fileName.find("c2_") + 3;
        size_t end = fileName.find_first_of("_/", start);
        if (end != std::string::npos) {
            c2_type += fileName.substr(start, end - start);
            c2_number = fileName.substr(start, end - start);
        }
      }
      if (fileName.find("qMax_") != std::string::npos) {
        c2_type += "_qMax_";
        size_t start = fileName.find("qMax_") + 5;
        size_t end = fileName.find_first_of("_./", start);
        if (end != std::string::npos) {
            qMax= fileName.substr(start, end - start);
        } else {
            qMax= fileName.substr(start);
        }
        c2_type += qMax;
      }
      if (fileName.find("REJ") != std::string::npos) {
        c2_type += "rej";
        rejected = "REJ";
      }
    }
    std::string ref_mix = "";
    if (fileName.substr(fileName.find("tev") + 3, 3) == "mix") {
      // Extract the substring between "tev" and "c2"
      std::cout<<"tev mix"<<std::endl;
      size_t start = fileName.find("tev") + 3; // Start after "tev"
      size_t end = fileName.find("_c2");       // End before "c2"
      if (end != std::string::npos && start < end) {
        std::string subStr = fileName.substr(start, end - start);
        if (subStr.find("theta") != std::string::npos) {
          ref_mix += "#theta";
        }
        if (subStr.find("a_pt") != std::string::npos) {
          ref_mix += "p_{T}";
        }
        if (subStr.find("rot") != std::string::npos) {
          ref_mix += "rot";
        }
      }
    }
    if (fileName.find("fit") != std::string::npos) {
      ref_num ++;
      if (ref_type.empty()) {
        ref_type = fileName.substr(fileName.find("fit") + 4, 3);
        ref_type += ref_mix;
        std::cout<< "ref_type: " << ref_type << std::endl;
        ref_type_now.push_back(ref_type);
      }
      else if (ref_mix.empty() && ref_type_now[ref_type_now.size()-1].compare(fileName.substr(fileName.find("fit") + 4, 3)) != 0) {
        ref_type += "vs" + fileName.substr(fileName.find("fit") + 4, 3);
        ref_num = 1;
        ref_type_now.push_back(fileName.substr(fileName.find("fit") + 4, 3));
      } else if (!ref_mix.empty() && (ref_type_now[ref_type_now.size()-1].substr(3, ref_type_now[ref_type_now.size()-1].size() - 3)).compare(ref_mix) != 0) {
        std::cout<<"ref now: "<<ref_type_now[ref_type_now.size()-1].substr(3, ref_type_now[ref_type_now.size()-1].size() - 3);
        ref_type += "vs" + ref_mix;
        ref_num = 1;
        std::cout<< "ref_type: " << ref_type << std::endl;
        ref_type_now.push_back(ref_type_now[ref_type_now.size()-1].substr(0, 3) + ref_mix);
        // ref_type_now.push_back(ref_type_now[ref_type_now.size()-1]);
      }
    }
    cout << "ref number: " <<ref_num <<endl;
    cout << "c2_type: " << c2_type << endl;
    cout << "c2_number: " << c2_number << endl;

    for (auto paramName : paramNameVec) {
      std::ifstream file(fileName);

      csvRow row;
      int paramCol = -1;
      int chiCol = -1;
      std::string chi2value = "";
      if (file >> row) {
        for (size_t j = 0; j < row.size(); ++j) {
          if (row[j].compare(paramName) == 0) {
            paramCol = j;
          }
          if (row[j].compare("#chi^{2}") == 0) {
            chiCol = j;
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
      if (plottt) {
        binCenterVec = binCenterVectt;
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
        if (plotMult) {
          if (binStart == 2 && binStop == 250) {
            break;
          }
          if (binStop > max_binStop)
            break;
          if (binStart > 200) {
            continue;
          }
          // if (binStop - binStart > 100) {
          //   break;
          // }
        } else if (plotKt) {
          // if (binStop - binStart > 1000) {
          //   break;
          // }
          if (binStart < 300 || binStart > 1500) {
            if (binStop < 2000)
              continue;
          }
        }
        val = std::stod(row[paramCol]);
        if (noUncert) {
          valErr = 0.;
        } else {
          valErr = std::stod(row[paramCol + 1]);
        }
        if (chiCol >= 0) {
          chi2value = row[chiCol];
        }
        if (chi2value == "-nan" || chi2value == "nan" || chi2value.empty()) {
           std::cerr << "WARNING: Invalid chi2value detected." << std::endl;
          valErr = 10.;
          if (paramName.compare("#chi^{2}") == 0) {
            val = 10.;
          }
        }

        binCenterIndex = GetVecIndex(binCenterVec, binStart, binStop);
        if (plotMult) {
          if (binStart == 21 && binStop == 250) {
            binCenterIndex = binCenterVec.size() - 1;
          }
        } else if (plotKt) {
          if (binStart == 100 && binStop == 2000) {
            binCenterIndex = binCenterVec.size() - 1;
          }
        }
        filteredBinCenterVec.emplace_back(binCenterVec.at(binCenterIndex));
        binErrVecLeft.emplace_back(filteredBinCenterVec.back() - binStart);
        binErrVecRight.emplace_back(binStop - filteredBinCenterVec.back());
        valVec.emplace_back(val);
        valErrVec.emplace_back(valErr);
      }

      // if (plotMult) {
      {
        auto rightArrow = new TArrow(
              filteredBinCenterVec.back() - binErrVecLeft.back(),
              valVec.back(),
              filteredBinCenterVec.back() + binErrVecRight.back(),
              valVec.back(),
              0.02, ">");
        //copy rightarrow into rightarrowhead
        if (binCenterIndex == binCenterVec.size() - 1) {
          rightArrow->SetLineStyle(2);
        }
        arrowVec.emplace_back(rightArrow);

        filteredBinCenterVec.emplace_back(filteredBinCenterVec.back() + binErrVecRight.back());
        valVec.emplace_back(-1000.);
        binErrVecLeft.emplace_back(0.);
        binErrVecRight.emplace_back(0.);
        valErrVec.emplace_back(0.);

        binErrVecLeft.end()[-2] = 0.;
        binErrVecRight.end()[-2] = 0.;
      }
      auto graph = new TGraphAsymmErrors(valVec.size(),
                                         &filteredBinCenterVec[0], &valVec[0],
                                         &binErrVecLeft[0], &binErrVecRight[0],
                                         &valErrVec[0], &valErrVec[0]);
      if (xLabel.empty()) {
        if (plotMult) {
          graph->GetXaxis()->SetTitle("N_{ch}");
        }
        if (plotKt) {
          graph->GetXaxis()->SetTitle("k_{T} [MeV]");
        }
        if (plottt) {
          graph->GetXaxis()->SetTitle("t_{T}");
        }
      }
      graph->GetYaxis()->SetTitle(yLabel.c_str());
      if (titleIndex < titleVec.size()) {
        graph->SetTitle(titleVec.at(titleIndex).c_str());
      } else {
        auto parname = paramName;
        parname += "_{" + rejected + c2_number + "}";
        std::string refTypeUpper = ref_type_now[ref_type_now.size()-1].substr(0, 3);
        std::transform(refTypeUpper.begin(), refTypeUpper.end(), refTypeUpper.begin(), ::toupper);
        parname += "(" + refTypeUpper + ref_mix + ": ";
        if (chi2value == "-nan" || chi2value == "nan" || chi2value.empty()) {
          parname += chi2value;
        } else {
          float chi2valueFloat = std::stof(chi2value);
          // Convert the float back to string with 2 decimal places
          char chi2valueStr[10];
          snprintf(chi2valueStr, sizeof(chi2valueStr), "%.2f", chi2valueFloat);
          parname += chi2valueStr;
        }
        parname += ")";
        if (qMax.empty()) {
          parname += "  [Q<4.2 GeV]";
          graph->SetTitle(parname.c_str());
        } else {
          float qMaxGeV = std::stof(qMax) / 1000.0f;
          // Convert the float back to string with 2 decimal places
          char qMaxGeVStr[10];
          snprintf(qMaxGeVStr, sizeof(qMaxGeVStr), "%.2f", qMaxGeV);
          parname += "  [Q<"+ std::string(qMaxGeVStr)+"GeV]";
          graph->SetTitle(parname.c_str());
        }
        qMaxVec.emplace_back(qMax);
      }
      graphVec.emplace_back(graph);
      ++titleIndex;
    }
  }
  //check if all qMax are not the same, if yes detele them from parname title and create only textbox on righttop of graph
  if (std::adjacent_find(qMaxVec.begin(), qMaxVec.end(), std::not_equal_to<>()) == qMaxVec.end()) {
    for (auto graph : graphVec) {
      std::string parname = graph->GetTitle();
      size_t pos = parname.find("[Q<");
      //now create the textbox on righttop of graph
      auto textBox = new TLatex();
      textBox->SetTextSize(0.02);
      textBox->SetTextAlign(12);
      textBox->SetTextColor(kBlack);
      textBox->SetTextFont(42);
      textBox->SetNDC();
      textBox->SetX(0.7); // X position (near the right edge)
      textBox->SetY(0.92); // Y position (near the top edge)
      std::string qMaxtext = "[Q<" + qMaxVec[0] + " GeV]";
      textBox->SetTitle(qMaxtext.c_str());
      // textBox->DrawLatex(0.7, 0.9, qMaxtext.c_str());
      textVec.emplace_back(textBox); // Store the TLatex object globally
      if (pos != std::string::npos) {
        parname.erase(pos, parname.length());
        graph->SetTitle(parname.c_str());
      }
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
    plotName += "_" + ref_type;
    if (plotKt) {
      plotName += "_kt";
    }
  }

  Plotter plotter(outFilePath + "/plot" + plotName);
  for (auto graph : graphVec) {
    plotter.addGraph(graph, paramNameVec.size() * ref_num);
    std::cout << "paramnamesize: " << paramNameVec.size() << ", ref_num: " << ref_num << std::endl;
  }
  for (auto arrow : arrowVec) {
    plotter.addArrow(arrow, paramNameVec.size() * ref_num);
  }
  plotter.addNotes(commentVec);
  // plotter.legendY1 = .75;
  // plotter.legendX1 = .1;
  // plotter.legendX2 = .45;
  for (auto text : textVec) {
    plotter.addText(text);
  }
  plotter.legendX1 = legX[0];
  plotter.legendX2 = legX[1];
  plotter.legendY1 = legY[0];
  plotter.legendY2 = legY[1];
  plotter.yMin = yMin;
  plotter.yMax = yMax;
  plotter.logY = logY;
  plotter.draw();
}
