//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronCalibration is a class designed to perform the cosine correction
// for electrons. This entails taking the means and errors prepared by the 
// WidgetPrep class and interpolating in bins of cos(theta).
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef ELECTRONCOLLECTOR_H
#define ELECTRONCOLLECTOR_H

#include <string>
#include <iostream>
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

#include "ElectronCorrection.h"

static double PI = 3.14159265;

class Electroncollector{

 public:

  Electroncollector();
  virtual ~Electroncollector() {};

  // --------------------------------------------------
  // run gains
  // --------------------------------------------------

  // fit the dE/dx distributions
  void fitRunGains( TFile* outfile );

  // plot the mean values versus run number
  void plotRunGains( TString filename );

  // --------------------------------------------------
  // wire gains
  // --------------------------------------------------

  // --------------------------------------------------
  // electron saturation correction
  // --------------------------------------------------

  // perform the electron saturation correction
  void SaturationCorrection( TFile* oufile );

  // --------------------------------------------------
  // 2D correction
  // --------------------------------------------------

  // divide data into bins of DOCA and entrance angle
  void TwoDCorrection( TFile* outfile );

  // --------------------------------------------------
  // 1D residual cleanup
  // --------------------------------------------------

  // divide data into bins of entrance angle
  void make1DHists();

 private:

  TString m_filename; // name of ROOT file that contain hit information
  TString m_constfilename; // name of ROOT file that contains constants

  int m_mcFlag; // flag for hadron correction: 0 if data, 1 if mc
  int m_type; // flag for data type: 0 if BESIII, 1 if BELLEII
  int m_fits; // for 2D correction do truncation (0) or fits (1)

  int m_docabins; // number of DOCA bins
  double m_lowerdoca; // lower bound of DOCA
  double m_upperdoca; // upper bound of DOCA

  int m_entabins; // number of entrance angle bins
  double m_lowerenta; // lower bound of entrance angle
  double m_upperenta; // upper bound of entrance angle

  int m_costhbins; // number of cos(theta) bins
};
#endif
