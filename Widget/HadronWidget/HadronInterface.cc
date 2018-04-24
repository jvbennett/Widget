#include "HadronInterface.h"

HadronInterface::HadronInterface(){
  // set up a vector of particle names
  m_types.push_back("pion");
  m_types.push_back("kaon");
  m_types.push_back("proton");
  m_types.push_back("muon");
  m_types.push_back("electron");

  // if no files are specified, assume sample is BESIII MC
  m_mcFlag = 1;
  m_type = 0;

  // set up the values in the map of particle names
  for( int i = 0; i < 5; ++i ){
    m_index[m_types[i]] = i;
  }

  // and masses
  m_mass["pion"] = Widget::mpion;
  m_mass["kaon"] = Widget::mkaon;
  m_mass["proton"] = Widget::mproton;
  m_mass["muon"] = Widget::mmuon;
  m_mass["electron"] = Widget::melectron;

  // initialize various parameters
  for( int i = 0; i < 5; ++i ){
    m_filenames[i] = "";
    m_nevents[i] = 0;
    m_bgbins[i] = 15;
    m_upperbg[i] = 5.0;
    m_lowerbg[i] = 0.1;
    m_cosbins[i] = 18;    
    m_uppercos[i] = 0.93;
    m_lowercos[i] = 0.0;
  }

  m_paramfile = "parameters.txt";
}

void
HadronInterface::SetupFromConfigFile( std::string configfile ){

  cout << "HadronInterface: setting up from " << configfile << endl;

  std::ifstream input(configfile.c_str());

  int type;
  std::string varname, vartype, line;
  while( std::getline(input, line) ){

    // skip comment lines
    line.erase( std::find( line.begin(), line.end(), '#' ), line.end() );
    if( line.empty() ) continue;

    // parse the lines in the configuration file
    type = -1;
    std::istringstream iss(line);
    iss >> vartype;

    // check for a mcFlag
    if( vartype == "mcFlag" ){
      iss >> m_mcFlag;
      continue;
    }
    else if( vartype == "type" ){
      iss >> m_type;
      continue;
    }
    else if( vartype == "pion") // the first word in the line should be the particle type
      type = 0;
    else if( vartype == "kaon")
      type = 1;
    else if( vartype == "proton")
      type = 2;
    else if( vartype == "muon")
      type = 3;
    else if( vartype == "electron")
      type = 4;
    else{
      cout << "HadronInterface: ERROR: unknown particle type " << vartype << " - this line will not be read" << endl;
      continue;
    }

    // the next word should be the variable type
    iss >> varname;
    if( varname == "filename")
      iss >> m_filenames[type];
    else if( varname == "nevents")
      iss >> m_nevents[type];
    else if( varname == "bgbins")
      iss >> m_bgbins[type];
    else if( varname == "upperbg")
      iss >> m_upperbg[type];
    else if( varname == "lowerbg")
      iss >> m_lowerbg[type];
    else if( varname == "upperp"){
      iss >> m_upperbg[type];
      m_upperbg[type] = m_upperbg[type]/m_mass[m_types[type]];
    }
    else if( varname == "lowerp"){
      iss >> m_lowerbg[type];
      m_lowerbg[type] = m_lowerbg[type]/m_mass[m_types[type]];
    }
    else if( varname == "cosbins")
      iss >> m_cosbins[type];
    else if( varname == "uppercos")
      iss >> m_uppercos[type];
    else if( varname == "lowercos")
      iss >> m_lowercos[type];
    else{
      cout << "HadronInterface: ERROR: unknown variable type " << vartype << " - this line will not be read" << endl;
      continue;
    }
  }

  cout << "HadronInterface: Done parsing configuration file..." << endl;

  m_paramfile = "parameters.txt";
}

void
HadronInterface::PrintParameters(){

  cout << endl << endl << "Printing parameters" << endl << endl;

  for( int i = 0; i < 5; ++i ){
    cout << m_types[i] << endl;
    cout << "\tFilename: " << m_filenames[i] << endl;
    cout << "\tBeta-gamma range: " << 
      m_bgbins[i] << " bins, from " << m_upperbg[i] << " to " << m_lowerbg[i] << endl;
    cout << "\tCos(theta) range: " << 
      m_cosbins[i] << " bins, from " << m_uppercos[i] << " to " << m_lowercos[i] << endl;
  }
}

void
HadronInterface::SetParamFile( std::string paramfile ){
  m_paramfile = paramfile;
}

void
HadronInterface::AddParticle( TString particle, TString filename, int nevents, int bgbins, double upperbg, double lowerbg, int cosbins, double uppercos, double lowercos ){
  int index = m_index[particle];
  
  m_filenames[index] = filename;
  m_nevents[index] = nevents;
  m_bgbins[index] = bgbins;
  m_upperbg[index] = upperbg;
  m_lowerbg[index] = lowerbg;
  m_cosbins[index] = cosbins;
  m_uppercos[index] = uppercos;
  m_lowercos[index] = lowercos;
}

void
HadronInterface::QualityCheck( TString outfilename ) {
  
  cout << "HadronInterface: making quality plots..." << endl;

  std::vector<TString> filenames;
  for( int i = 0; i < 5; ++i ){
    filenames.push_back(m_filenames[i]);
  }

  WidgetQuality gq(m_type,m_upperbg,m_lowerbg);
  gq.makeHistograms(filenames, outfilename);
  gq.plotHistograms(outfilename);
}

void
HadronInterface::GenerateSample( TString particle = "pion" ) {
  
  cout << "HadronInterface: generating sample events..." << endl;

  int index = m_index[particle];

  // generate a sample of events
  WidgetGenerator generator(m_nevents[index], m_upperbg[index], m_lowerbg[index]);
  generator.generateEvents(m_filenames[index],m_mass[particle]);

  cout << "\t generated " << m_nevents[index] << "..." << endl;
}

void
HadronInterface::SimulateSample( TString infilename, TString particle ) {
  
  cout << "HadronInterface: simulating the dE/dx response for sample events..." << endl;

  int index = m_index[particle];

  // simulate the dE/dx response for a sample of events
  WidgetGenerator generator(m_nevents[index], m_upperbg[index], m_lowerbg[index]);
  generator.simulateDedx(infilename,m_filenames[index],m_mass[particle]);

  cout << "\t simulated " << m_nevents[index] << "..." << endl;
}

void
HadronInterface::SimulateReconstruction( TString infilename, TString outfilename ) {
  
  cout << "HadronInterface: simulating the dE/dx reconstruction..." << endl;

  // simulate the dE/dx response for a sample of events
  WidgetGenerator generator;
  generator.simulateReconstruction(infilename,outfilename);

  cout << "\t simulated reconstruction finished..." << endl;
}

// -------------------- HADRON CORRECTIONS --------------------


/* Fit the truncated mean distributions in bins of \beta\gamma and cos(theta) and 
   save the means and widths in the output file. */

void
HadronInterface::PrepareSample( TString outfilename, bool correct, int particle ) {

  cout << "HadronInterface: fitting in bins of beta-gamma and cos(theta)..." << endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  if( particle == -1 ){
    // make sure the samples exist
    for( int i = 0; i < 5; ++i ){ // <---- adding electrons...
      TFile* testfile = new TFile(m_filenames[i]);
      if( testfile->GetListOfKeys()->IsEmpty() ){
	cout << "\t No hadron sample for " << m_types[i] << "... " << endl;
	continue;
      }
      testfile->Close();

      // prepare the sample by fitting in bins of beta-gamma and cos(theta)
      HadronPrep prep(m_filenames[i], m_mcFlag, m_type, m_bgbins[i], m_upperbg[i], m_lowerbg[i], m_cosbins[i], m_uppercos[i], m_lowercos[i]);
      prep.bgCosThetaFits(m_types[i],outfile,correct,m_paramfile);
    }
  }
  else{
    // make sure the samples exist
    TFile* testfile = new TFile(m_filenames[particle]);
    if( testfile->GetListOfKeys()->IsEmpty() ){
      cout << "\t No hadron sample for " << m_types[particle] << "... " << endl;
    }
    testfile->Close();

    // prepare the sample by fitting in bins of beta-gamma and cos(theta)
    HadronPrep prep(m_filenames[particle], m_mcFlag, m_type, m_bgbins[particle], m_upperbg[particle], m_lowerbg[particle], m_cosbins[particle], m_uppercos[particle], m_lowercos[particle]);
    prep.bgCosThetaFits(m_types[particle],outfile,correct,m_paramfile);
  }

  outfile->Close();
}


/* Fit the truncated mean distributions in bins of \beta\gamma and save the means
   and widths in the output file. */

void
HadronInterface::PrepareResults( TString outfilename, bool correct ) {

  cout << "HadronInterface: getting the universal beta-gamma fits..." << endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  // make sure the samples exist
  for( int i = 0; i < 5; ++i ){
    TFile* testfile = new TFile(m_filenames[i]);
    if( testfile->GetListOfKeys()->IsEmpty() ){
      cout << "\t No hadron sample for " << m_types[i] << "... " << endl;
      continue;
    }
    testfile->Close();

    // prepare the sample by fitting in bins of beta-gamma and cos(theta)
    HadronPrep prep(m_filenames[i], m_mcFlag, m_type, m_bgbins[i], m_upperbg[i], m_lowerbg[i], m_cosbins[i], m_uppercos[i], m_lowercos[i]);
    prep.bgFits(m_types[i],outfile,correct,m_paramfile);
  }
  outfile->Close();
}

void
HadronInterface::SaturationCorrection( TString prepfilename, std::string parfile ) {

  // do not do the correction on mc samples...
  if( m_mcFlag ) return;

  cout << "HadronInterface: performing the hadron saturation correction..." << endl;

  // make sure the file exists
  TFile* testfile = new TFile(prepfilename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t CAUTION: Fitted file does not exist... " << endl;
    exit(1);
  }
  testfile->Close();

  // create the HadronSaturation object and fill the prepared samples
  HadronSaturation hadsat;
  hadsat.setParameters(parfile);
  hadsat.fillSample(prepfilename);
  hadsat.printEvents(0,10);

  // perform the hadron saturation correction
  hadsat.fitSaturation();
}

void
HadronInterface::BetaGammaFits( std::vector< TString > particles, TString filename ) {

  cout << "HadronInterface: performing the beta-gamma curve fits..." << endl;

  // make sure the file exists
  TFile* testfile = new TFile(filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t CAUTION: Fitted file does not exist... " << endl;
    exit(1);
  }
  testfile->Close();

  // create the HadronCalibration object and fit the prepared samples
  HadronCalibration hadcal;
  hadcal.fitBGCurve(particles,filename,m_paramfile);
  SetParamFile("parameters.bgcurve.fit");
}

void
HadronInterface::PlotBGCurve( std::vector< TString > particles, TString filename ) {

  cout << "HadronInterface: performing the beta-gamma curve fits..." << endl;

  // make sure the file exists
  TFile* testfile = new TFile(filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t CAUTION: Fitted file does not exist... " << endl;
    exit(1);
  }
  testfile->Close();

  // create the HadronCalibration object and fit the prepared samples
  HadronCalibration hadcal;
  hadcal.plotBGCurve(particles,filename,m_paramfile);
}

void
HadronInterface::SigmaFits( TString filename ) {

  cout << "HadronInterface: performing the sigma curve fits..." << endl;

  // make sure the file exists
  TFile* testfile = new TFile(filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t CAUTION: Fitted file does not exist... " << endl;
    exit(1);
  }
  testfile->Close();

  // create the HadronCalibration object and fit the prepared samples
  // only the pion samples are used for the sigma fits
  HadronCalibration hadcal;
  hadcal.fitSigmaVsNHit(m_filenames[0],m_paramfile,m_mcFlag,m_type);
  SetParamFile("parameters.sigmanhit.fit");
  hadcal.fitSigmaVsSin(m_filenames[0],m_paramfile,m_mcFlag,m_type);
  SetParamFile("parameters.sigmasin.fit");
}


void
HadronInterface::PlotEfficiency(TString saveas) {

  cout << "HadronInterface: plotting the efficiency..." << endl;

  // create the HadronCalibration object and fill histograms
  HadronCalibration hadcal;
  SetParamFile("parameters.bgcurve.fit");

  TString flavor("pion");
  hadcal.plotEfficiency(m_filenames,saveas,m_paramfile,m_mcFlag,m_type);
}
