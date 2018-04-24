//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronCorrection is a program that applies electron dE/dx corrections.
//
// For additional details, see the Widget document.
//
//********************************************************************************

#include <fstream>
#include <iostream>
#include <streambuf>
#include <string>

#include "Widget/ElectronWidget/ElectronInterface.h"

int main( int argc, char* argv[] ) {

  std::cout << std::endl << "*** Running the electron calibration ***" << std::endl << std::endl;

  std::string configfile("electron.template.txt");
  int ask = 0;

  // ***************************
  // Parse command line options
  // ***************************

  for (int i = 1; i < argc; i++){

    std::string arg(argv[i]);

    if (arg == "-c"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-i"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  ask = atoi(argv[++i]); }
    if (arg == "-h"){
      std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
      std::cout << "\t -c <file>\t Configuration file" << std::endl;
      std::cout << "\t -i  <int>\t Ask before each step? (0 = no, 1 = yes)" << std::endl;
      exit(1);
    }
  }

  if( configfile == "" ){
    std::cout << "ERROR: No configuration file was given..." << std::endl;
    return 1;
  }

  // Hesitate after each if ask flag is true
  std::string ready("y");

  // Create a calibration interface
  ElectronInterface widget;


  // ---------------------------------------------------------------------------
  // 0. Set the parameters either from a config file

  if( configfile != "" )
    widget.SetupFromConfigFile(configfile);
  else{
    std::cout << "ElectronCorrection ERROR: No configuration file provided";
    return 1;
  }
  // ---------------------------------------------------------------------------
  std::cout << "Done configuring settings..." << std::endl;

  // redirect output to log file, but save a pointer to the old
  // buffer associated with std::cout and restore it for some output
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::ofstream out("widget.electron.log");


  // ---------------------------------------------------------------------------
  // 1. Apply corrections

  if( ask == 1 ){
    std::cout << "Do you want to apply corrections (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout << "Applying corrections..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.ApplyCorrections("electron.corrected.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  std::cout << std::endl << std::endl << "DONE ELECTRON CORRECTION" << std::endl;

  out.close();
}
