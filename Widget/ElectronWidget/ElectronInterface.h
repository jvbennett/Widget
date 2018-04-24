//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronInterface is a class designed to be the point of contact with the user.
// This class interfaces with the other Widget classes to perform the generation,
// simulation, fitting, and correction tasks necessary for dE/dx calibration.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef ELECTRONINTERFACE_H
#define ELECTRONINTERFACE_H

#include <fstream>
#include <algorithm>

#include "WidgetHelpers/WidgetConstants.h"
#include "ElectronGenerator.h"
#include "ElectronCalibration.h"
#include "ElectronCorrection.h"

using std::cout;
using std::endl;

class ElectronInterface{

 public:

  ElectronInterface();
  virtual ~ElectronInterface() {}

  // set the member variables from the given configuration file
  void SetupFromConfigFile( std::string configfile );

  // get the input file name
  TString GetFileName(){ return m_filename; }

  // change the input file (eg. to calibrated sample)
  void SetFileName( TString newfilename );

  // generate an electron sample
  void GenerateSample( TString outfile );

  // ---------- methods for the electron calibration ----------

  // prepare and fit the saturation (dE/dx vs. cos(theta))
  void SaturationCorrection( TString outfilename );

  // prepare and fit the run gains
  void RunGains( TString outfilename );

  // prepare and fit the 2-d correction (doca vs. entrance angle)
  void TwoDCorrection( TString outfilename );

  // apply calibrations
  void ApplyCorrections( TString outfile );

 private:

  TString m_filename; // name of ROOT file that contain hit information
  TString m_constfilename; // name of file that contains constants

  double m_removeLowest;  // low end of truncation
  double m_removeHighest; // high end of truncation

  int m_mcFlag; // flag for hadron correction: 0 if data, 1 if mc
  int m_type; // flag for data type: 0 if BESIII, 1 if BELLEII
  int m_fits; // flag to indicate how to get means: 0 for truncation, 1 for fits

  int m_nevents; // the number of events stored in the vectors below

  int m_bgbins; // number of BG bins
  double m_lowerbg; // lower bound of BG
  double m_upperbg; // upper bound of BG

  int m_docabins; // number of DOCA bins
  double m_lowerdoca; // lower bound of DOCA
  double m_upperdoca; // upper bound of DOCA

  int m_entabins; // number of entrance angle bins
  double m_lowerenta; // lower bound of entrance angle
  double m_upperenta; // upper bound of entrance angle

  int m_costhbins; // number of cos(theta) bins

  ElectronCorrection ecor; // object to access corrections
};
#endif
