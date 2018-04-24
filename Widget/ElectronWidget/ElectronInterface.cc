#include "ElectronInterface.h"

ElectronInterface::ElectronInterface(){

  // if no files are specified, assume sample is MC (1) and not data (0)
  m_mcFlag = 1;

  // if no files are specified, assume sample is BESIII (0) and not Belle II (1)
  m_type = 0;

  // initialize various parameters
  m_filename = "";
  m_constfilename = "";
  m_nevents = 0;

  // default truncation
  m_removeLowest = 0.05;
  m_removeHighest = 0.25;

  // specify whether to take truncated means (0) or perform fits (1)
  m_fits = 0;
}

void
ElectronInterface::SetupFromConfigFile( std::string configfile ){

  cout << "ElectronInterface: setting up from " << configfile << endl;

  std::ifstream input(configfile.c_str());

  std::string varname, vartype, line;
  while( std::getline(input, line) ){

    // skip comment lines
    line.erase( std::find( line.begin(), line.end(), '#' ), line.end() );
    if( line.empty() ) continue;

    // parse the lines in the configuration file
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
    else if( vartype == "fit" ){
      iss >> m_fits;
      continue;
    }
    else if( vartype == "constfile" ){
      iss >> m_constfilename;
      continue;
    }
    else if( vartype == "lowest" ){
      iss >> m_removeLowest;
      continue;
    }
    else if( vartype == "highest" ){
      iss >> m_removeHighest;
      continue;
    }
    else if( vartype != "electron" ){ // the first word in the line should be the particle type
      std::cout << "ElectronInterface: ERROR: unknown particle type " << vartype << " - this line will not be read" << endl;
      continue;
    }

    // the next word should be the variable type
    iss >> varname;
    if( varname == "filename")
      iss >> m_filename;
    else if( varname == "nevents")
      iss >> m_nevents;
    else if( varname == "docabins")
      iss >> m_docabins;
    else if( varname == "upperdoca")
      iss >> m_upperdoca;
    else if( varname == "lowerdoca")
      iss >> m_lowerdoca;
    else if( varname == "bgbins")
      iss >> m_bgbins;
    else if( varname == "upperp"){
      iss >> m_upperbg;
      m_upperbg = m_upperbg/Widget::melectron;
    }
    else if( varname == "lowerp"){
      iss >> m_lowerbg;
      m_lowerbg = m_lowerbg/Widget::melectron;
    }
    else if( varname == "entabins")
      iss >> m_entabins;
    else if( varname == "upperenta")
      iss >> m_upperenta;
    else if( varname == "lowerenta")
      iss >> m_lowerenta;
    else if( varname == "costhbins")
      iss >> m_costhbins;
    else{
      std::cout << "ElectronInterface: ERROR: unknown variable type " << vartype << " - this line will not be read" << endl;
      continue;
    }
  }

  std::cout << "ElectronInterface: Done parsing configuration file..." << endl;
}


void
ElectronInterface::SetFileName( TString newfilename ){
  m_filename = newfilename;
}


void
ElectronInterface::GenerateSample( TString outfilename ){

  std::cout << "ElectronInterface: generating electron sample..." << std::endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");

  // bin the data in run numbers and fit for the mean in each bin
  ElectronGenerator egen(100000,2.0,0.2);
  egen.generateEvents(outfile);

  // if we are generating a sample, use it for calibration
  m_filename = outfilename;
}



// -------------------- ELECTRON CORRECTIONS --------------------

void
ElectronInterface::RunGains( TString outfilename ) {

  std::cout << "ElectronInterface: getting run gains..." << std::endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  // make sure the samples exist
  TFile* testfile = new TFile(m_filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t No electron sample" << endl;
    exit;
  }
  testfile->Close();

  // bin the data in run numbers and fit for the mean in each bin
  ElectronCalibration ecalib(m_filename, m_constfilename, m_mcFlag, m_type, m_fits, m_docabins, m_upperdoca, m_lowerdoca, m_entabins, m_upperenta, m_lowerenta, m_costhbins);
  ecalib.fitRunGains(outfile);

  // plot the run gains
  ecalib.plotRunGains(outfilename);
}

void
ElectronInterface::TwoDCorrection( TString outfilename ) {

  std::cout << "WidgetInterface: getting two dimensional correction..." << std::endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  // make sure the samples exist
  TFile* testfile = new TFile(m_filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t No electron sample" << endl;
    exit;
  }
  testfile->Close();

  // bin the data in doca and entrance angle and fit for the mean in each bin
  ElectronCalibration ecalib(m_filename, m_constfilename, m_mcFlag, m_type, m_fits, m_docabins, m_upperdoca, m_lowerdoca, m_entabins, m_upperenta, m_lowerenta, m_costhbins);
  ecalib.TwoDCorrection(outfile);

  outfile->Close();
}

void 
ElectronInterface::SaturationCorrection( TString outfilename ) {

  std::cout << "WidgetInterface: doing electron saturation correction..." << std::endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  // make sure the samples exist
  TFile* testfile = new TFile(m_filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t No electron sample" << endl;
    exit;
  }
  testfile->Close();

  // bin the data in doca and entrance angle and fit for the mean in each bin
  ElectronCalibration ecalib(m_filename, m_constfilename, m_mcFlag, m_type, m_fits, m_docabins, m_upperdoca, m_lowerdoca, m_entabins, m_upperenta, m_lowerenta, m_costhbins);
  ecalib.SaturationCorrection(outfile);

  outfile->Close();
}

void
ElectronInterface::ApplyCorrections( TString outfilename ) {

  std::cout << "WidgetInterface: applying electron corrections..." << std::endl;

  TFile* outfile = new TFile(outfilename,"RECREATE");
  // make sure the samples exist
  TFile* testfile = new TFile(m_filename);
  if( testfile->GetListOfKeys()->IsEmpty() ){
    cout << "\t No electron sample" << endl;
    exit;
  }
  testfile->Close();

  // bin the data in doca and entrance angle and fit for the mean in each bin
  ElectronCorrection ecor(m_constfilename,m_filename,m_removeLowest,m_removeHighest);
  ecor.process(outfile);

  outfile->Close();
}
