//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// HadronCalibration is a program that performs the generation, simulation, 
// fitting, and correction tasks necessary for dE/dx calibration.
//
// For additional details, see the Widget document.
//
//********************************************************************************

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFitter.h"
#include "TLegend.h"
#include "TApplication.h"

#include "Widget/HadronWidget/HadronInterface.h"

int main( int argc, char* argv[] ) {

  std::cout << std::endl << "*** Running the calibration ***" << std::endl << std::endl;

  std::string configfile("");
  TString inputfile("");
  int part = -1;
  std::string parfile("parameters.txt");
  int ask = 1;

  // ***************************
  // Parse command line options
  // ***************************

  for (int i = 1; i < argc; i++){

    std::string arg(argv[i]);

    if (arg == "-c"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-f"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  inputfile = argv[++i]; }
    if (arg == "-p"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  parfile = argv[++i]; }
    if (arg == "-t"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  part = atoi(argv[++i]); }
    if (arg == "-i"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  ask = atoi(argv[++i]); }
    if (arg == "-h"){
      std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
      std::cout << "\t -c <file>\t Configuration file" << std::endl;
      std::cout << "\t -f <file>\t Input file" << std::endl;
      std::cout << "\t -p <file>\t Curve parameter file" << std::endl;
      std::cout << "\t -t  <int>\t Particle type (-1 = all)" << std::endl;
      std::cout << "\t -i  <int>\t Ask before each step? (0 = no, 1 = yes)" << std::endl;
      exit(1);
    }
  }

  if( configfile == "" && inputfile == "" ){
    std::cout << "ERROR: No input file was given..." << std::endl;
    return 1;
  }

  // Hesitate after each if ask flag is true
  std::string ready("y");

  // Create a calibration interface
  HadronInterface widget;

  std::vector< TString > types;
  types.push_back("pion");
  types.push_back("kaon");
  types.push_back("proton");
  types.push_back("muon");
  types.push_back("electron");

  // ---------------------------------------------------------------------------
  // 0. Set the parameters either from a config file or use the defaults

  if( configfile != "" ){
    widget.SetupFromConfigFile(configfile);
    widget.PrintParameters();
  }
  else{
    std::cout << "No configuration file provided, ";
    std::cout << "are you sure you want to use the default values (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
    if( ready != "y" ) return 1;

    int bgbins[5];
    double bgmax[5];
    double bgmin[5];

    // BESIII XYZ 2014 DATA
    bgbins[0] = 30;
    bgmax[0] = 2.12/Widget::mpion;
    bgmin[0] = 0.12/Widget::mpion;
    bgbins[1] = 10;
    bgmax[1] = 1.8/Widget::mkaon;
    bgmin[1] = 0.35/Widget::mkaon;
    bgbins[2] = 10;
    bgmax[2] = 1.7/Widget::mproton;
    bgmin[2] = 0.5/Widget::mproton;
    bgbins[3] = 15;
    bgmax[3] = 1.92/Widget::mmuon;
    bgmin[3] = 1.76/Widget::mmuon;
    bgbins[4] = 30;
    bgmax[4] = 2.3/Widget::melectron;
    bgmin[4] = 0.3/Widget::melectron;

    int cosbins[5];
    double cosmax[5];
    double cosmin[5];
    for( int i = 0; i < 4; ++i ){
      cosbins[i] = 18;
      cosmax[i] = 0.93;
      cosmin[i] = 0.0;
    }
    cosbins[4] = 18;
    cosmax[4] = 0.93;
    cosmin[4] = 0.5;
    
    // add particles
    for( int i = 0; i < 5; ++i ){
      if( inputfile != "" )
	widget.AddParticle(types[i],inputfile+types[i]+".root",0.0,
			   bgbins[i],bgmax[i],bgmin[i],
			   cosbins[i],cosmax[i],cosmin[i]);
      else
	widget.AddParticle(types[i],"widget."+types[i]+".root",0.0,
			   bgbins[i],bgmax[i],bgmin[i],
			   cosbins[i],cosmax[i],cosmin[i]);
    }
  }

  // set the file that holds the calibration constants
  widget.SetParamFile(parfile);
  // ---------------------------------------------------------------------------

  // redirect output to log file, but save a pointer to the old
  // buffer associated with std::cout and restore it for some output
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::ofstream out("widget.log");

  // ---------------------------------------------------------------------------
  // 1. Make histograms for quality validation and monitoring

  if( ask == 1 ){
    std::cout << "Do you want to make quality plots (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" )
    widget.QualityCheck( "widget.quality.root" );
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 2. Do fits of truncated means in bins of beta-gamma and cos(theta)
  //   Do not include the electron sample here

  if( ask == 1 ){
    std::cout << "Do you want to fit in bins of beta-gamma and cos(theta) (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    if( part == -1 ) // using pre-defined sample (default)
      widget.PrepareSample( "widget.uncorrected.root", false, -1 );
    else if( part == -2 ){ // no samples pre-defined so use fake samples
      for( int i = 0; i < 4; ++i ){
	widget.PrepareSample( "widget.uncorrected.root", false, i );
	std::cout.rdbuf( out.rdbuf() );
	std::cout << "Ready to move on to the next particle (y/n)? ";
	std::cin >> ready; std::cout << std::endl;
	std::cout.rdbuf( coutbuf );
	if( ready != "y" ) return 1;
      }
      return 1; // short circuit here for now
    }
    else{ // supplying individual sample for diagnostic purposes
      widget.PrepareSample( "widget.uncorrected.root", false, part );
      return 1; // short circuit here for now
    }
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 3. Determine the calibration constants for the hadron correction.
  //   SaturationCorrection does nothing for MC samples.

  if( ask == 1 ){
    std::cout << "Do you want to do to the hadron correction (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.SaturationCorrection( "widget.uncorrected.root", "sat-pars.txt" );
    widget.PrepareSample( "widget.corrected.root", true, part );
    std::cout.rdbuf( coutbuf );
  }

  // ---------------------------------------------------------------------------
  // 4. Repeat the 2D correction using the previous results as a seed
  //   SaturationCorrection does nothing for MC samples.

  if( ask == 1 ){
    std::cout << "Do you want to repeat to the hadron correction (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.SaturationCorrection( "widget.uncorrected.root", "sat-pars.fit.txt" );
    widget.PrepareSample( "widget.corrected.root", true, part );
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 5. Fit the corrected mean values in bins of beta-gamma

  if( ask == 1 ){
    std::cout << "Do you want to perform the fits in bins of beta-gamma (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.PrepareResults("widget.corrected_bgbins.root",true);
    widget.SetParamFile(parfile);
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 6. Fit the means in bins of beta-gamma and use the results to 
  //   fit the beta-gamma curve and determine the curve parameters
  if( ask == 1 ){
    std::cout << "Do you want to fit the beta-gamma curve (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.BetaGammaFits(types,"widget.corrected_bgbins.root");
    widget.SigmaFits("widget.corrected.root");
    widget.SigmaFits("widget.corrected.root");
    widget.SigmaFits("widget.corrected.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 6b. Make some plots
  if( ask == 1 ){
    std::cout << "Do you want to plot the beta-gamma curve (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.PlotBGCurve(types,"widget.corrected_bgbins.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // 7. Fit the corrected mean values in bins of beta-gamma and repeat
  //   the fits for the beta-gamma and sigma curves

  if( ask == 1 ){
    std::cout << "Do you want to repeat the fits in bins of beta-gamma (y/n)? ";
    std::cin >> ready; std::cout << std::endl;
  }
  if( ask == 0 || ready == "y" ){
    std::cout.rdbuf( out.rdbuf() );
    widget.PrepareResults("widget.corrected_bgbins.root",true);
    widget.BetaGammaFits(types,"widget.corrected_bgbins.root");
    widget.SigmaFits("widget.corrected.root");
    widget.PrepareResults("widget.corrected_bgbins.root",true);
    widget.BetaGammaFits(types,"widget.corrected_bgbins.root");
    std::cout.rdbuf( coutbuf );
  }
  // ---------------------------------------------------------------------------
}
