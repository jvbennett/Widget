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
  TString saveas("");

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
    if (arg == "-s"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  saveas = argv[++i]; }
    if (arg == "-h"){
      std::cout << std::endl << " Usage for: " << argv[0] << std::endl << std::endl;
      std::cout << "\t -c <file>\t Configuration file" << std::endl;
      std::cout << "\t -s <file>\t Save output as" << std::endl;
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
  // 1. Plot the efficiency for dE/dx

  std::cout.rdbuf( out.rdbuf() );
  widget.PlotEfficiency(saveas);
  std::cout.rdbuf( coutbuf );
  // ---------------------------------------------------------------------------
}
