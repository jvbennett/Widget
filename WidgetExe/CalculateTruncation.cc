//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// CalculateTruncation recalculates the truncated means for a sample of tracks.
//
// For additional details, see the Widget document.
//
//********************************************************************************

#include <string>
#include "Widget/WidgetHelpers/WidgetGenerator.h"

int main( int argc, char* argv[] ) {

  std::cout << std::endl << "*** Recalculating truncation ***" << std::endl << std::endl;

  TString inputfile("");
  TString outputfile("");

  // ***************************
  // Parse command line options
  // ***************************

  for (int i = 1; i < argc; i++){

    std::string arg(argv[i]);

    if (arg == "-f"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  inputfile = argv[++i]; }
    if (arg == "-o"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outputfile = argv[++i]; }
    if (arg == "-h"){
      std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
      std::cout << "\t -f <file>\t Input file" << std::endl;
      std::cout << "\t -o <file>\t Output file" << std::endl;
      exit(1);
    }
  }

  if( inputfile == "" || outputfile == "" ){
    std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
    std::cout << "\t -c <file>\t Configuration file" << std::endl;
    std::cout << "\t -f <file>\t Input file" << std::endl;
    std::cout << "\t -o <file>\t Output file" << std::endl;
    exit(1);
  }

  // Create a widget interface and simulate the trucation
  WidgetGenerator widget;
  widget.simulateReconstruction(inputfile,outputfile);
}
