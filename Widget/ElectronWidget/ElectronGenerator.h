//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronGenerator is a class designed to generate or modify an electron sample 
// for use in the WidgetPrep class.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef ELECTRONGENERATOR_H
#define ELECTRONGENERATOR_H

#include <string>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"

#include "WidgetHelpers/WidgetConstants.h"

class ElectronGenerator{

 public:

  ElectronGenerator();
  ElectronGenerator( int nevents, double upperbg, double lowerbg );
  virtual ~ElectronGenerator() {};

  // generate a sample of fake tracks
  void generateEvents( TFile* outfile );

 private:

  int m_nevents; // the number of events stored in the vectors below

  double m_upperbg; // upper bound on beta-gamma
  double m_lowerbg;  // lower bound on beta-gamma

};
#endif
