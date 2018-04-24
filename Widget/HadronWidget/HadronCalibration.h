//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// HadronCalibration is a class designed to perform the hadron saturation
// correction. This entails taking the means and errors prepared by the 
// WidgetPrep class and fitting in bins of cos(theta).
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef HADRONCALIBRATION_H
#define HADRONCALIBRATION_H

#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFitResultPtr.h"

#include "HadronSaturation.h"
#include "WidgetHelpers/WidgetParameterization.h"
#include "WidgetHelpers/WidgetConstants.h"

class HadronCalibration{

 public:

  HadronCalibration();
  virtual ~HadronCalibration() {};

  // fit the beta-gamma curve
  void fitBGCurve( std::vector< TString > particles, TString filename, std::string paramfile );

  // fit the beta-gamma curve
  void plotBGCurve( std::vector< TString > particles, TString filename, std::string paramfile );

  // fit sigma vs. nhit
  void fitSigmaVsNHit( TString filename, std::string paramfile, int mcFlag, int type );

  // fit sigma vs. sin(theta)
  void fitSigmaVsSin( TString filename, std::string paramfile, int mcFlag, int type );

  // plot efficiency
  void plotEfficiency( TString filenames[5], TString saveas, std::string paramfile, int mcFlag, int type );

 private:

  HadronSaturation m_hadsat;

};
#endif
