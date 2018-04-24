//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// HadronPrep is a class that prepares a sample for the hadron calibration. It
// first bins the sample in terms of beta-gamma (p/m) and cos(theta) and fits for
// the mean of the measured dE/dx distribution and its error. This information is
// used in the HadronSaturation class. The HadronPrep class also contains methods
// to determine the predicted dE/dx mean and resolution.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef WIDGETPREP_H
#define WIDGETPREP_H

#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "HadronSaturation.h"
#include "HadronCalibration.h"

class HadronPrep{

 public:

  HadronPrep();
  virtual ~HadronPrep() {};

  // set the input file name and other variables
  HadronPrep( TString infile, int mcflag, int type, int bgbins, double upperbg, double lowerbg, int cosbins, double uppercos, double lowercos );

  // apply some formatting to TGraphs
  void FormatGraph( TGraphErrors* gr, int flag );

  // perform fits in bins of beta-gamma and cos(theta)
  void bgCosThetaFits( TString particle, TFile* outfile, bool correct, std::string paramfile );

  // perform fits in bins of beta-gamma
  void bgFits( TString particle, TFile* outfile, bool correct, std::string paramfile );

 private:

  TString m_filename; // ROOT file containing reconstructed dedx information
  int m_mcFlag; // flag for hadron correction: 0 if data, 1 if mc
  int m_type; // flag for data type: 0 if BESIII, 1 if BELLEII

  int m_bgbins; // number of beta-gamma bins
  double m_upperbg; // Upper limit on beta-gamma
  double m_lowerbg; // Lower limit on beta-gamma

  int m_cosbins;  // number of cos(theta) bins
  double m_uppercos; // Upper limit on cos(theta)
  double m_lowercos; // Lower limit on cos(theta)

  HadronSaturation m_hadsat; // object storing hadron saturation parameterization
  WidgetParameterization m_gpar; // defines beta-gamma curve and other parameterizations
};
#endif
