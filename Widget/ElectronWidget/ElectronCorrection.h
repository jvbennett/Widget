//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// ElectronCorrection is a class designed to apply calibration constants for the
// electron calibration.
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef ELECTRONCORRECTION_H
#define ELECTRONCORRECTION_H

#include "HadronWidget/HadronSaturation.h"

class ElectronCorrection{

 public:

  ElectronCorrection();
  ElectronCorrection( TString constfilename, TString infilename, double removeLowest, double removeHighest );
  virtual ~ElectronCorrection() {};

  /** Retrieve the calibration constants (placeholder for now) */
  void initializeParameters();

  /** Apply corrections to events in infile and write to outfile */
  void process(TFile* outfile);

  /** Perform a run gain correction */
  double RunGainCorrection(double& dedx) const;
  /** Perform a wire gain correction */
  double WireGainCorrection(int wireID, double& dedx) const;

  /** Perform a standard set of corrections */
  double StandardCorrection(int wireID, double costheta, double dedx) const;

  /** Perform a hadron saturation correction.
   * (Set the peak of the truncated mean for electrons to 1) */
  double HadronCorrection(double costheta,  double dedx) const;

  /** Perform the truncation */
  void calculateMeans(double& mean, double& truncatedMean, double& truncatedMeanErr, double dedx[], int size) const;
  void calculateMeans(double& mean, double& truncatedMean, double& truncatedMeanErr, std::vector<double> dedx, int size) const;

 private:

  TString m_filename;      // name of input ROOT file
  TString m_constfilename; // name of file that contains constants

  double m_removeLowest;  // low end of truncation
  double m_removeHighest; // high end of truncation

  double m_runGain;   // run gain constant

  const int c_NCDCWires = 14336; // number of wires in CDC
  double m_wireGain[14336]; // list of wire gains
  double m_valid[14336];     // whether each wire is good or not

  HadronSaturation m_hadsat; // hadron saturation object
};
#endif
