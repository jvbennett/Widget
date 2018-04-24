//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// WidgetQuality is a class that creates some basic validation plots that are 
// useful for determining the momentum range for a sample.
//
// Plots include: cos(theta), momentum, momentum vs. cos(theta), nhit vs. pt,
//   and dE/dx vs. momentum (with lines defining mom. range used in calibration)
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef WIDGETQUALITY_H
#define WIDGETQUALITY_H

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"

#include "WidgetConstants.h"

class WidgetQuality{

 public:

  WidgetQuality();
  WidgetQuality(int type, double* bgmax, double* bgmin);
  virtual ~WidgetQuality() {};

  // make histograms for validation (no plotting)
  void makeHistograms( std::vector<TString> filenames, TString outfilename );

  // plot histograms for validation
  void plotHistograms( TString outfilename );

 private:

  int m_type; // flag for MC type: 0 if BESIII, 1 if BELLEII

  double m_mass[5]; // particle masses
  double m_upperp[5]; // upper bound on the momentum
  double m_lowerp[5];  // lower bound on the momentum
};
#endif
