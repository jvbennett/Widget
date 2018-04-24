#include "HadronSaturation.h"

static HadronSaturation *HC_obj;

HadronSaturation::HadronSaturation() :
  m_flag(0),
  m_cosbins(18),
  m_alpha(0.0266879),
  m_gamma(0.00378704),
  m_delta(0.0038601),
  m_power(1.8957),
  m_ratio(3.09834)
{
  m_dedx.clear();
  m_dedxerror.clear();
  m_betagamma.clear();
  m_costheta.clear();
  HC_obj = this;
}

HadronSaturation::HadronSaturation( double alpha, double gamma, double delta, double power, double ratio ) :
  m_flag(0),
  m_cosbins(18),
  m_alpha(alpha),
  m_gamma(gamma),
  m_delta(delta),
  m_power(power),
  m_ratio(ratio)
{
  m_dedx.clear();
  m_dedxerror.clear();
  m_betagamma.clear();
  m_costheta.clear();
}

void
HadronSaturation::setParameters( double par[] ){
  m_alpha = par[0];
  m_gamma = par[1];
  m_delta = par[2];
  m_power = par[3];
  m_ratio = par[4];
}

void
HadronSaturation::setParameters( std::string infile ){
  double hadpar[5];
  std::ifstream parfile(infile.c_str());
  if( !parfile.fail() ){
    for( int i = 0; i < 5; ++i ){
      parfile >> hadpar[i];
    }
    parfile.close();
    setParameters(hadpar);
  }
  else
    std::cout << "HadronSaturation: WARNING - parameter file " << infile << " does not exist, using defaults..." << std::endl;
}

void
HadronSaturation::fillSample( TString infilename ) {

  std::cout << "HadronSaturation: filling sample events..." << std::endl;

  // clear the vectors to be filled
  clear();
  
  std::vector< TString > types;
  types.push_back("pion");
  types.push_back("kaon");
  types.push_back("proton");
  types.push_back("muon");
  types.push_back("electron");

  // fill the containers with a previously prepared sample
  TFile* satFile = new TFile(infilename);

  for( int i = 0; i < 4; ++i ){
    TTree* satTree = (TTree*)satFile->Get(types[i]);

    double satdedx, saterror, satbg, satcosth;
    satTree->SetBranchAddress("dedx",&satdedx);
    satTree->SetBranchAddress("error",&saterror);
    satTree->SetBranchAddress("bg",&satbg);
    satTree->SetBranchAddress("costh",&satcosth);

    // fill the vectors
    for( unsigned int i = 0; i < satTree->GetEntries(); ++i ){
      satTree->GetEvent(i);

      if( saterror == 0 ) continue;

      m_dedx.push_back(satdedx);
      m_dedxerror.push_back(saterror);
      m_betagamma.push_back(satbg);
      m_costheta.push_back(satcosth);
    }
  }

  satFile->Close();
}

void
HadronSaturation::printEvents( int firstevent = 0, int nevents = 10 ){

  if( firstevent < 0 || (nevents+firstevent) > m_dedx.size() ){
    std::cout << "HadronSaturation: trying to print events out of range (" << m_dedx.size() << ")" << std::endl;
    exit(1);
  }

  std::cout << std::endl << std::endl;
  for( int i = firstevent; i < nevents; ++i ){
    std::cout << "Event " << i << ":\t bg = " << m_betagamma[i] << "\t cos(theta) = " << m_costheta[i] << "\t dE/dx mean = " << m_dedx[i] << "\t dE/dx error = " << m_dedxerror[i] << std::endl;
  }
  std::cout << std::endl;
}

double 
HadronSaturation::D2I( double cosTheta, double D, double alpha, double gamma, double delta, double power, double ratio ) const {

  double absCosTheta   = fabs(cosTheta);
  double projection    = pow(absCosTheta, power) + delta;
  double chargeDensity = D / projection;
  double numerator     = 1 + alpha * chargeDensity;
  double denominator   = 1 + gamma * chargeDensity;

  double I = D * ratio * numerator / denominator;

  return I;
}

double 
HadronSaturation::I2D( double cosTheta, double I, double alpha, double gamma, double delta, double power, double ratio ) const {

  double absCosTheta = fabs(cosTheta);
  double projection  = pow(absCosTheta, power) + delta;

  double a =  alpha/projection;
  double b =  1 - gamma/projection * (I/ratio);
  double c = -1.0 * I/ratio;

  if( b == 0 && a == 0 ){
    std::cout << "both a and b coefficiants for hadron correction are 0" << std::endl;
    return I;
  }

  double D = (a != 0) ? (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a) : -c / b;
  if( D < 0 ){
    //    std::cout << "D is less 0! will try another solution" << std::endl;
    D = (a != 0) ? (-b - sqrt(b * b + 4.0 * a * c)) / (2.0 * a) : -c / b;
    if( D < 0 ){
      //      std::cout << "D is still less 0! just return uncorrectecd value" << std::endl;
      return I;
    }
  }

  return D;
}

double
HadronSaturation::myFunction( double alpha, double gamma, double delta, double power, double ratio ){

  unsigned int nevt = m_dedx.size();
  double chisq = 0, dedxsum = 0;
  std::vector< double > vdedxavg;

  // Compute the average value (across cos(theta)) for each bin of beta-gamma.
  // NOTE: the correction is not constrained to a certain value (1), it
  // changes as a function of beta gamma...
  double dedxcor = 0;
  for( unsigned int i = 0; i < nevt; i++ ){
    if( m_flag == 0 )
      dedxcor = D2I(m_costheta[i],I2D(m_costheta[i],1.0,alpha,gamma,delta,power,ratio)*m_dedx[i],alpha,gamma,delta,power,ratio);
    else
      dedxcor = D2I(m_costheta[i],m_dedx[i],alpha,gamma,delta,power,ratio);
    dedxsum += dedxcor;

    if( (i+1)%m_cosbins == 0 ){
      vdedxavg.push_back(dedxsum/m_cosbins);
      dedxsum = 0;
    }
  }

  // Construct a chi^2 value for the difference between the average in a cosine bin
  // to the actual values
  for( unsigned int i = 0; i < nevt; i++ ){
    if( m_flag == 0 )
      dedxcor = D2I(m_costheta[i],I2D(m_costheta[i],1.0,alpha,gamma,delta,power,ratio)*m_dedx[i],alpha,gamma,delta,power,ratio);
    else
      dedxcor = D2I(m_costheta[i],m_dedx[i],alpha,gamma,delta,power,ratio);

    int j = (int)i/m_cosbins;
    chisq += pow((dedxcor-vdedxavg[j])/m_dedxerror[i],2);
  }

  return chisq;
}

void
HadronSaturation::minuitFunction( int& nDim, double* gout, double& result, double* par, int flag ){
  result = HC_obj->myFunction(par[0],par[1],par[2],par[3],par[4]);
}

void
HadronSaturation::fitSaturation() {

  std::cout << "Performing the hadron saturation fit..." << std::endl;

  // Construct the fitter
  TFitter* minimizer = new TFitter(5);
  minimizer->SetFCN(HadronSaturation::minuitFunction);

  TRandom* rand = new TRandom();
  
  minimizer->SetParameter(0,"alpha",m_alpha,rand->Rndm()*0.001,0.0,0.5);
  minimizer->SetParameter(1,"gamma",m_gamma,rand->Rndm()*0.001,-2,2);
  minimizer->SetParameter(2,"delta",m_delta,rand->Rndm()*0.001,-2,2);
  minimizer->SetParameter(3,"power",m_power,rand->Rndm()*0.01,-5,5);
  minimizer->SetParameter(4,"ratio",m_ratio,rand->Rndm()*0.01,-100,100);

  // Set minuit fitting strategy
  double arg[10];
  arg[0] = 2;
  minimizer->ExecuteCommand("SET STR",arg,1);

  // Suppress much of the output of MINUIT
  arg[0] = -1;
  minimizer->ExecuteCommand("SET PRINT",arg,1);

  // Minimize with MIGRAD
  arg[0] = 10000;

  double fitpar[5], fiterr[5];
  for( int i = 0; i < 20; ++i ){

    minimizer->SetParameter(0,"alpha",m_alpha,rand->Rndm()*0.001,0.0,0.5);
    minimizer->SetParameter(1,"gamma",m_gamma,rand->Rndm()*0.001,-2,2);
    minimizer->SetParameter(2,"delta",m_delta,rand->Rndm()*0.001,-2,2);
    minimizer->SetParameter(3,"power",m_power,rand->Rndm()*0.01,-5,5);
    minimizer->SetParameter(4,"ratio",m_ratio,rand->Rndm()*0.01,-100,100);

    int status = minimizer->ExecuteCommand("MIGRAD",arg,1);
    status = minimizer->ExecuteCommand("MIGRAD",arg,1);
    minimizer->PrintResults(1,0);
    std::cout << "Fit status: " << status << std::endl;

    int counter = 0;
    while( status != 0 && counter < 10 ){
      counter++;

      minimizer->SetParameter(0,"alpha",m_alpha,rand->Rndm()*0.001,0.0,0.5);
      minimizer->SetParameter(1,"gamma",m_gamma,rand->Rndm()*0.001,-2,2);
      minimizer->SetParameter(2,"delta",m_delta,rand->Rndm()*0.001,-2,2);
      minimizer->SetParameter(3,"power",m_power,rand->Rndm()*0.01,-5,5);
      minimizer->SetParameter(4,"ratio",m_ratio,rand->Rndm()*0.01,-100,100);

      status = minimizer->ExecuteCommand("MIGRAD",arg,1);
      std::cout << "Fit status: " << status << std::endl;
    }
    if( status != 0 ){
      std::cout << "HadronSaturation::ERROR - BAD FIT!" << std::endl;
      return;
    }
  }

  for( int par = 0; par < 5; ++par ){
    fitpar[par] = minimizer->GetParameter(par);
    fiterr[par] = minimizer->GetParError(par);
  }

  std::cout << "HadronSaturation: fit results" << std::endl;
  std::cout << "\t alpha (" << m_alpha << "): " << fitpar[0] << " +- " << fiterr[0] << std::endl;
  std::cout << "\t gamma (" << m_gamma << "): " << fitpar[1] << " +- " << fiterr[1] << std::endl;
  std::cout << "\t delta (" << m_delta << "): " << fitpar[2] << " +- " << fiterr[2] << std::endl;
  std::cout << "\t power (" << m_power << "): " << fitpar[3] << " +- " << fiterr[3] << std::endl;
  std::cout << "\t ratio (" << m_ratio << "): " << fitpar[4] << " +- " << fiterr[4] << std::endl;

  // function value, estimated distance to minimum, errdef
  // variable parameters, number of parameters
  double chi2, edm, errdef;
  int nvpar, nparx;
  minimizer->GetStats(chi2,edm,errdef,nvpar,nparx);
  std::cout << "\t\t Fit chi^2: " << chi2 << std::endl;
  
  std::ofstream parfile;
  parfile.open("sat-pars.fit.txt");

  parfile << fitpar[0] << std::endl;
  parfile << fitpar[1] << std::endl;
  parfile << fitpar[2] << std::endl;
  parfile << fitpar[3] << std::endl;
  parfile << fitpar[4] << std::endl;
  parfile.close();

  delete minimizer;
}

void
HadronSaturation::clear() {

  m_dedx.clear();
  m_dedxerror.clear();
  m_betagamma.clear();
  m_costheta.clear();
}
