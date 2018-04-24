//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// HadronInterface is a class designed to be the point of contact with the user.
// This class interfaces with the other Widget classes to perform the generation,
// simulation, fitting, and correction tasks necessary for dE/dx calibration.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef HADRONINTERFACE_H
#define HADRONINTERFACE_H

#include <string>
#include <iostream>
#include <fstream>
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

#include "WidgetHelpers/WidgetQuality.h"
#include "WidgetHelpers/WidgetGenerator.h"
#include "HadronCalibration.h"
#include "HadronPrep.h"

using std::cout;
using std::endl;

class HadronInterface{

 public:

  HadronInterface();
  virtual ~HadronInterface() {}

  // set the member variables from the given configuration file
  void SetupFromConfigFile( std::string configfile );

  // set the name of the file containing the calibration parameters
  void SetParamFile( std::string paramfile );

  // print the HadronInterface parameters
  void PrintParameters();

  // declare the number of events and beta-gamma range for each particle type
  void AddParticle( TString particle, TString filename, int nevents, int bgbins, double upperbg, double lowerbg, int cosbins, double uppercos, double lowercos );

  // generate a set of histograms for validation/monitoring
  void QualityCheck( TString outfilename );

  // generate a sample of events and fill the vectors
  void GenerateSample( TString particle );

  // simulate a sample of events and fill the vectors
  void SimulateSample( TString infilename, TString particle );

  // simulate a sample of events and fill the vectors
  void SimulateReconstruction( TString infilename, TString outfilename );

  // ---------- methods for the hadron calibration ----------

  // prepare the sample by fitting as a function of beta-gamma and cos(theta)
  void PrepareSample( TString outfilename, bool correct, int particle );

  // prepare the sample by fitting as a function of beta-gamma
  void PrepareResults( TString outfilename, bool correct );

  // perform the hadron saturation correction
  void SaturationCorrection( TString prepfilename, std::string parfile );

  // perform the beta-gamma fits
  void BetaGammaFits( std::vector< TString > particles, TString filename );

  // make some plots
  void PlotBGCurve( std::vector< TString > particles, TString filename );

  // perform the sigma fits (vs. nhit and sin(theta))
  void SigmaFits( TString filename );

  // plot efficiency
  void PlotEfficiency( TString saveas );

 private:

  TString m_filenames[5]; // names of ROOT files that contain hit information
  int m_nevents[5]; // the number of events stored in the vectors below
  int m_mcFlag; // flag for hadron correction: 0 if data, 1 if mc
  int m_type; // flag for MC type: 0 if BESIII, 1 if BELLEII

  std::string m_paramfile; // file containing calibration parameters

  std::vector< TString > m_types; // vector of particle types
  std::map< TString, int > m_index; // map particle name to vector index
  std::map< TString, double > m_mass; // map particle name to particle mass

  int m_bgbins[5]; // number of beta-gamma bins
  double m_upperbg[5]; // upper bound on beta-gamma
  double m_lowerbg[5];  // lower bound on beta-gamma

  int m_cosbins[5];  // number of cos(theta) bins
  double m_uppercos[5]; // upper bound on cos(theta)
  double m_lowercos[5];  // lower bound on cos(theta)
};
#endif
