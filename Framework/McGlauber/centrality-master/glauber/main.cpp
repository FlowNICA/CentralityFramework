#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Fitter.h"
#include "FitterHelper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStopwatch.h"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "No argumets provided!" << std::endl;
    std::cerr << "./glauber macro.C" << std::endl;
    std::cerr << "Arguments:" << std::endl;
    std::cerr << "\tmacro.C - macro with main parameters" << std::endl;
    std::cerr << std::endl;
    return 1;
  }
  std::string macro{argv[1]};

  TStopwatch timer;
  timer.Start();

  std::cout << "glauber: Executing " + macro << std::endl;
  gROOT->Macro(macro.c_str());

  timer.Stop();
  timer.Print();
  return 0;
}