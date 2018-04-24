//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// WidgetGenerator is a class designed to generate or modify a sample for use
// in the WidgetPrep class. If a pre-existing sample is provided, it is possible
// to replace the dE/dx information with a fake distribution.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef WIDGETGENERATOR_H
#define WIDGETGENERATOR_H

#include <string>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"

#include "HadronWidget/HadronSaturation.h"
#include "HadronWidget/HadronCalibration.h"

class HadronSaturation;

class WidgetGenerator{

 public:

  WidgetGenerator();
  WidgetGenerator( int nevents, double upperbg, double lowerbg );
  virtual ~WidgetGenerator() {};

  // generate a sample of fake tracks
  void generateEvents( TString filename, double genmass );

  // replace the dE/dx measurements for a sample with a fake distribution
  void simulateDedx( TString infilename, TString filename, double genmass );

  // take the hit level information and replicate the truncation
  void simulateReconstruction( TString infilename, TString outfilename );

 private:

  int m_nevents; // the number of events stored in the vectors below

  double m_upperbg; // upper bound on beta-gamma
  double m_lowerbg;  // lower bound on beta-gamma

};
#endif
