//********************************************************************************
// This file is part of the Widget, a package for performing dE/dx calibration.
//
// Author: Jake Bennett
// Date: July 8, 2015
//
//
// For additional details, see the Widget document.
//
//********************************************************************************
#ifndef WIDGETSIGMA_H
#define WIDGETSIGMA_H

#include <math.h>
#include <iostream>

class WidgetSigma{

 public:

  WidgetSigma() {}
  virtual ~WidgetSigma() {};

  double sigmaCurve( double *x, double *par ){

    // sometimes we want to constrain certain parameters
    double f = 0;
    if( par[0] == 1 ){
      // return dedx parameterization
      f = par[1]+par[2]*x[0];
    }
    else if( par[0] == 2 ){
      // return nhit or sin(theta) parameterization
      f = par[1]*pow(x[0],4)+par[2]*pow(x[0],3)+
	par[3]*x[0]*x[0]+par[4]*x[0]+par[5];
    }

    return f;
  }

  double operator()( double *x, double *par ){
    return sigmaCurve(x,par);
  }
};
#endif
