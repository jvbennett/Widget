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
#ifndef WIDGETPARAMETERIZATION_H
#define WIDGETPARAMETERIZATION_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "WidgetHelpers/WidgetCurve.h"
#include "WidgetSigma.h"

class WidgetParameterization{

 public:

  WidgetParameterization();
  WidgetParameterization( std::string parfile );
  virtual ~WidgetParameterization() {};

  // read in parameters from file
  void setParameters( std::string infile );

  // write out parameters to file
  void printParameters( std::string infile );

  // write out parameters to screen
  void showParameters();

  // determine the predicted mean dE/dx
  double dedxPrediction( double x );

  // determine the predicted mean dE/dx resolution
  double sigmaPrediction( double dedx, double nhit, double sin );

  // determine the predicted mean dE/dx resolution without dedx dependence
  double resPrediction( double nhit, double sin );

  // determine the predicted mean dE/dx resolution vs. sin(theta)
  double sinPrediction( double sin );

  // determine the predicted mean dE/dx resolution vs. nhit
  double nhitPrediction( double nhit );

  // set the curve parameters
  double getCurvePars( int i ){ return m_curvepars[i]; };

  // set the curve parameters
  void setCurvePars( int i, double val ){ m_curvepars[i] = val; };

  // set the dedx parameters
  double getDedxPars( int i ){ return m_dedxpars[i]; };

  // set the dedx parameters
  void setDedxPars( int i, double val ){ m_dedxpars[i] = val; };

  // set the nhit parameters
  double getNHitPars( int i ){ return m_nhitpars[i]; };

  // set the nhit parameters
  void setNHitPars( int i, double val ){ m_nhitpars[i] = val; };

  // set the sin(theta) parameters
  double getSinPars( int i ){ return m_sinpars[i]; };

  // set the sin(theta) parameters
  void setSinPars( int i, double val ){ m_sinpars[i] = val; };

 private:

  double m_curvepars[16]; // parameters for beta-gamma curve
  double m_dedxpars[3]; // parameters for sigma vs. dE/dx curve
  double m_sinpars[5]; // parameters for sigma vs. sin(theta) curve
  double m_nhitpars[5]; // parameters for sigma vs. nhit curve
};
#endif
