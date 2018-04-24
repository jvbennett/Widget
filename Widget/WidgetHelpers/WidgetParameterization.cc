#include "WidgetParameterization.h"

WidgetParameterization::WidgetParameterization(){
  setParameters("parameters.txt");
}

WidgetParameterization::WidgetParameterization( std::string infile ){
  setParameters(infile);
}

void
WidgetParameterization::setParameters( std::string infile ){

  std::cout << "Setting parameters..." << std::endl;

  std::ifstream fin;
  fin.open(infile.c_str());
  if( !fin.good() ){
    std::cout << "WARNING: CANNOT FIND parameters.txt!" << std::endl;
    exit(1);
  }

  int par;
  std::string stemp;
  fin >> par >> par;
  for( int i = 0; i <= 15; ++i ){
    fin >> par >> m_curvepars[i];
  }
  fin >> par >> par;
  for( int i = 0; i <= 2; ++i ){
    fin >> par >> m_dedxpars[i];
  }
  for( int i = 0; i <= 4; ++i ){
    fin >> par >> m_nhitpars[i];
  }
  for( int i = 0; i <= 4; ++i ){
    fin >> par >> m_sinpars[i];
  }

  fin.close();
}

void
WidgetParameterization::printParameters( std::string outfile ){

  std::cout << "WidgetParameterization: Printing parameters to file..." << std::endl;

  // write out the parameters to file
  std::ofstream fout( outfile.c_str() );
  fout << "0\t2" << std::endl << std::endl;
  for( int i = 1; i <= 16; ++i ){
    fout << i << "\t" << m_curvepars[i-1] << std::endl;
  }
  fout << std::endl << 0 << "\t" << 6 << std::endl << std::endl;
  for( int i = 1; i <= 3; ++i ){
    fout << i << "\t" << m_dedxpars[i-1] << std::endl;
  }
  fout << std::endl;
  for( int i = 4; i <= 8; ++i ){
    fout << i << "\t" << m_nhitpars[i-4] << std::endl;
  }
  fout << std::endl;
  for( int i = 9; i <= 13; ++i ){
    fout << i << "\t" << m_sinpars[i-9] << std::endl;
  }
  fout << std::endl;
  fout << "14\t1.105" << std::endl;

  fout.close();
}

void
WidgetParameterization::showParameters(){

  std::cout << "WidgetParameterization: Printing parameters to file..." << std::endl;

  std::cout << "0\t2" << std::endl << std::endl;
  for( int i = 1; i <= 16; ++i ){
    std::cout << i << "\t" << m_curvepars[i-1] << std::endl;
  }
  std::cout << std::endl << 0 << "\t" << 6 << std::endl << std::endl;
  for( int i = 1; i <= 3; ++i ){
    std::cout << i << "\t" << m_dedxpars[i-1] << std::endl;
  }
  std::cout << std::endl;
  for( int i = 4; i <= 8; ++i ){
    std::cout << i << "\t" << m_nhitpars[i-4] << std::endl;
  }
  std::cout << std::endl;
  for( int i = 9; i <= 13; ++i ){
    std::cout << i << "\t" << m_sinpars[i-9] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "14\t1.105" << std::endl;
}

double
WidgetParameterization::dedxPrediction( double bg ){
  // define the section of the curve to use
  double A = 0, B = 0, C = 0;
  if( bg < 4.5 )
    A = 1;
  else if( bg < 10 )
    B = 1;
  else
    C = 1;

  double x[1]; x[0] = bg;
  double parsA[9];
  double parsB[5];
  double parsC[5];

  parsA[0] = 1; parsB[0] = 2; parsC[0] = 3;
  for( int i = 0; i < 16; ++i ){
    if( i < 8 ) parsA[i+1] = m_curvepars[i];
    else if( i < 12 ) parsB[i%8+1] = m_curvepars[i];
    else parsC[i%12+1] = m_curvepars[i];
  }

  // calculate dE/dx from the Bethe-Bloch curve
  WidgetCurve gc;
  double partA = gc.bgCurve(x,parsA);
  double partB = gc.bgCurve(x,parsB);
  double partC = gc.bgCurve(x,parsC);

  return (A*partA+B*partB+C*partC);
}

double
WidgetParameterization::sigmaPrediction( double dedx, double nhit, double sin ){
  if( nhit < 5 ) nhit = 5;
  if( sin > 0.99 ) sin = 0.99;
  
  double x[1];
  double dedxpar[3];
  double nhitpar[6];
  double sinpar[6];

  dedxpar[0] = 1; nhitpar[0] = 2; sinpar[0] = 2;
  for( int i = 0; i < 2; ++i ){
    dedxpar[i+1] = m_dedxpars[i];
  }
  for( int i = 0; i < 5; ++i ){
    nhitpar[i+1] = m_nhitpars[i];
  }
  for( int i = 0; i < 5; ++i ){
    sinpar[i+1] = m_sinpars[i];
  }

  // determine sigma from the parameterization
  WidgetSigma gs;
  x[0] = dedx;
  double cor_dedx = gs.sigmaCurve(x,dedxpar);
  x[0] = nhit;
  double cor_nhit = gs.sigmaCurve(x,nhitpar);
  if( nhit > 35) cor_nhit = 1.0;
  x[0] = sin;
  double cor_sin = gs.sigmaCurve(x,sinpar);

  return cor_dedx * cor_sin * cor_nhit;
}

double
WidgetParameterization::resPrediction( double nhit, double sin ){
  if( nhit < 5 ) nhit = 5;
  if( sin > 0.99 ) sin = 0.99;
  
  double x[1];
  double nhitpar[6];
  double sinpar[6];

  nhitpar[0] = 2; sinpar[0] = 2;
  for( int i = 0; i < 5; ++i ){
    nhitpar[i+1] = m_nhitpars[i];
  }
  for( int i = 0; i < 5; ++i ){
    sinpar[i+1] = m_sinpars[i];
  }

  // determine sigma from the parameterization
  WidgetSigma gs;
  x[0] = nhit;
  double cor_nhit = gs.sigmaCurve(x,nhitpar);
  x[0] = sin;
  double cor_sin = gs.sigmaCurve(x,sinpar);

  if( nhit > 35) cor_nhit = 1.0;

  return cor_sin * cor_nhit;
}


double
WidgetParameterization::sinPrediction( double sin ){
  if( sin > 0.99 ) sin = 0.99;
  
  double x[1];
  double sinpar[6];

  sinpar[0] = 2;

  for( int i = 0; i < 5; ++i ){
    sinpar[i+1] = m_sinpars[i];
  }

  // determine sigma from the parameterization
  WidgetSigma gs;
  x[0] = sin;
  double cor_sin = gs.sigmaCurve(x,sinpar);

  return cor_sin;
}


double
WidgetParameterization::nhitPrediction( double nhit ){
  if( nhit < 5 ) nhit = 5;
  
  double x[1];
  double nhitpar[6];

  nhitpar[0] = 2;
  for( int i = 0; i < 5; ++i ){
    nhitpar[i+1] = m_nhitpars[i];
  }

  // determine sigma from the parameterization
  WidgetSigma gs;
  x[0] = nhit;
  double cor_nhit = gs.sigmaCurve(x,nhitpar);

  if( nhit > 35) cor_nhit = 1.0;

  return cor_nhit;
}
