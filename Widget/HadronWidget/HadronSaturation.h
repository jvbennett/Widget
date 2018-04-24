//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
// HadronSaturation is a class designed to perform the hadron saturation
// correction. This entails taking the means and errors prepared by the 
// WidgetPrep class and fitting in bins of cos(theta).
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef HADRONCORRECTION_H
#define HADRONCORRECTION_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFitter.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

class HadronSaturation{

 public:

  HadronSaturation();
  HadronSaturation( double alpha, double gamma, double delta, double power, double ratio );
  virtual ~HadronSaturation() { clear(); };

  // set the parameters
  void setParameters( double par[] );

  // set the parameters
  void setParameters( std::string parfile );

  // fill the vectors below
  void fillSample( TString infilename );

  // print a sample of events
  void printEvents( int firstevent, int nevents );

  // perform the hadron saturation fit
  void fitSaturation();

  // clear the vectors
  void clear();

  // set the number of cosine bins
  void setCosBins( int nbins ){ m_cosbins = nbins; }

  // set the hadron saturation function flag
  void setFlag( int flag ){ m_flag = flag; }

  // some helper functions for the hadron saturation correction
  double myFunction( double alpha, double gamma, double delta, double power, double ratio );
  static void minuitFunction( int& nDim, double* gout, double& result, double* para, int flg );
  
  double D2I( double cosTheta, double D, double alpha, double gamma, double delta, double power, double ratio ) const;
  double I2D( double cosTheta, double I, double alpha, double gamma, double delta, double power, double ratio ) const;

  // define these here to be used in other classes
  double D2I( double cosTheta, double D = 1 ) const {

    double absCosTheta   = fabs(cosTheta);
    double projection    = pow(absCosTheta, m_power) + m_delta;
    double chargeDensity = D / projection;
    double numerator     = 1 + m_alpha * chargeDensity;
    double denominator   = 1 + m_gamma * chargeDensity;

    double I = D * m_ratio * numerator / denominator;

    return I;
  }

  double I2D( double cosTheta, double I = 1 ) const {

    double absCosTheta = fabs(cosTheta);
    double projection  = pow(absCosTheta, m_power) + m_delta;

    double a =  m_alpha/projection;
    double b =  1 - m_gamma/projection * (I/m_ratio);
    double c = -1.0 * I/m_ratio;

    if( b == 0 && a == 0 ){
      std::cout << "both a and b coefficiants for hadron correction are 0" << std::endl;
      return I;
    }

    double D = (a != 0) ? (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a) : -c / b;
    if( D < 0 ){
      std::cout << "D is less 0! will try another solution" << std::endl;
      D = (a != 0) ? (-b - sqrt(b * b + 4.0 * a * c)) / (2.0 * a) : -c / b;
      if( D < 0 ){
	std::cout << "D is still less 0! just return uncorrectecd value" << std::endl;
	return I;
      }
    }

    return D;
  }

 private:

  // flag for saturation function
  int m_flag;

  // the number of cosine bins
  int m_cosbins;

  // the parameters for the hadron saturation correction
  double m_alpha;
  double m_gamma;
  double m_delta;
  double m_power;
  double m_ratio;

  std::vector< double > m_dedx;      // a vector to hold dE/dx measurements
  std::vector< double > m_dedxerror; // a vector to hold dE/dx errors
  std::vector< double > m_betagamma; // a vector to hold beta-gamma values
  std::vector< double > m_costheta;  // a vector to hold cos(theta) values

};
#endif
