//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronCalibration is a program that performs the generation, simulation, 
// fitting, and correction tasks necessary for dE/dx "electron calibration".
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
  std::string ask("y");

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
      else  ask = argv[++i]; }
    if (arg == "-h"){
      std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
      std::cout << "\t -c <file>\t Configuration file" << std::endl;
      std::cout << "\t -i  <int>\t Ask before each step (y/n)" << std::endl;
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
    std::cout << "ElectronCalibration ERROR: No configuration file provided";
    return 1;
  }
  // ---------------------------------------------------------------------------
  std::cout << "Done configuring settings. Output redirected to log file." << std::endl;

  // redirect output to log file, but save a pointer to the old
  // buffer associated with std::cout and restore it for some output
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::ofstream out("widget.electron.log");


  // ---------------------------------------------------------------------------
  // 0b. Generate samples

  if( ask == "y" ){
    std::cout << "Do you want to generate the sample (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ready == "y" ){
    std::cout << "Generating the sample..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.GenerateSample("electrons.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // Get the name of the original input file (or the generated one)
  TString originfile = widget.GetFileName();



  // ---------------------------------------------------------------------------
  // Apply old constants

  if( ask == "y" ){
    std::cout << "Do you want to apply the old constants (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ready == "y" ){
    std::cout << "Applying old constants..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.ApplyCorrections("electrons.precalib.root");
    widget.SetFileName("electrons.precalib.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------



  // ---------------------------------------------------------------------------
  // 1. Run gains

  if( ask == "y" ){
    std::cout << "Do you want to determine the run gains (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ready == "y" ){
    std::cout << "Getting run gains..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.RunGains("rungains.root");

    // reapply constants with new run gains
    widget.SetFileName(originfile);
    widget.ApplyCorrections("electrons.rungains.root");
    widget.SetFileName("electrons.rungains.root");

    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------



  // ---------------------------------------------------------------------------
  // 2. Do fits of doca vs. entrance angle

  if( ask == "y" ){
    std::cout << "Do you want to do the 2d correction (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ready == "y" ){
    std::cout << "Doing 2d correction..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.TwoDCorrection("twodcorrection.root");

    // reapply constants with new 2D correction
    widget.SetFileName(originfile);
    widget.ApplyCorrections("electrons.twod.root");
    widget.SetFileName("electrons.twod.root");

    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // 3. Do fits of dE/dx vs. cos(theta): saturation correction

  if( ask == "y" ){
    std::cout << "Do you want to do the saturation correction (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ready == "y" ){
    std::cout << "Doing saturation correction..." << std::endl;

    std::cout.rdbuf( out.rdbuf() );
    widget.SaturationCorrection("saturation.root");

    // reapply constants with new saturation correction
    widget.SetFileName(originfile);
    widget.ApplyCorrections("electrons.saturation.root");
    widget.SetFileName("electrons.saturation.root");

    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  std::cout << std::endl << std::endl << "DONE ELECTRON CALIBRATION" << std::endl;

  out.close();
}
