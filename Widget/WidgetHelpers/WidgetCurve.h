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
#ifndef WIDGETCURVE_H
#define WIDGETCURVE_H

#include <math.h>

class WidgetCurve{

 public:

  WidgetCurve() {};
  virtual ~WidgetCurve() {};

  double bgCurve( double *x, double *par ){
    
    // ************* FIX ME ******************
    // redundant parameters 3 and 5
    // ***************************************
    double f = 0;
    if( par[0] == 1 ) 
      f = par[1]*pow(sqrt(x[0]*x[0]+1),par[3])/pow(x[0],par[3]) * 
	(par[2]-par[6]*log(pow(1/x[0],par[4])) ) - par[5]+exp(par[7]+par[8]*x[0]);
    else if( par[0] == 2 ) 
      f = par[1]*pow(x[0],3)+par[2]*x[0]*x[0]+par[3]*x[0]+par[4];
    else if( par[0] == 3 ) 
      f = -1.0*par[1]*log(par[4]+pow(1/x[0],par[2]))+par[3];

    return f;
  }

  double operator()( double *x, double *par ){
    return bgCurve(x,par);
  }
};
#endif
