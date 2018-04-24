#include "ElectronCorrection.h"

ElectronCorrection::ElectronCorrection() : 
  m_constfilename("dEdxConstants.root"),
  m_filename(""),
  m_removeLowest(0.05),
  m_removeHighest(0.25)
{
  initializeParameters();
}


ElectronCorrection::ElectronCorrection( TString constfilename, TString infilename, double removeLowest, double removeHighest ) : 
  m_constfilename(constfilename),
  m_filename(infilename),
  m_removeLowest(removeLowest),
  m_removeHighest(removeHighest)
{
  initializeParameters();
}


void ElectronCorrection::process(TFile* outFile)
{

  TFile* dedxInput = new TFile(m_filename);
  TTree* track = (TTree*)dedxInput->Get("track");

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx truncated mean with electron saturation correction
  double dedxerr; // dE/dx error with electron saturation correction
  double mean;    // dE/dx mean with electron saturation correction
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  int nhits;      // number of hits on this track

  const int maxhits = 100;
  int wire[maxhits];    // wire number of this hit
  double dedxhit[maxhits]; // dE/dx for this hit

  track->SetBranchAddress("mean", &mean);
  track->SetBranchAddress("dedx", &dedxpub);
  track->SetBranchAddress("dedxsat", &dedx);
  track->SetBranchAddress("dedxerr", &dedxerr);
  track->SetBranchAddress("pF", &p);
  track->SetBranchAddress("costh", &costh);
  track->SetBranchAddress("numGoodLayerHits", &nhits);
  track->SetBranchAddress("wire", wire);
  track->SetBranchAddress("dedxhit", dedxhit);

  outFile->cd();
  TTree* newtrack = track->CloneTree(0);

  // --------------------------------------------------
  //  LOOP OVER EACH DEDX MEASUREMENT (TRACK LEVEL)
  // --------------------------------------------------

  for( unsigned int i = 0; i < track->GetEntries(); ++i ){
    track->GetEvent(i);

    if( nhits <= 0 ){
      std::cout << "No good hits on this track...";
      continue;
    }


    // --------------------------------------------------
    //  LOOP OVER EACH DEDX MEASUREMENT (HIT LEVEL)
    // --------------------------------------------------

    for (int j = 0; j < nhits; ++j) {
      double newdedx = StandardCorrection(wire[j], costh, dedxhit[j]);
      dedxhit[j] = newdedx;
    } // end loop over hits

    // recalculate means -> mean, dedx, and dedxerr will be replaced
    std::cout << "MEANS = " << dedx << "\t";
    calculateMeans(mean, dedx, dedxerr, dedxhit, nhits);
    std::cout << dedx << std::endl << std::endl;

    newtrack->Fill();
  } // end loop over tracks

  outFile->cd();
  newtrack->AutoSave();
  outFile->Close();
  dedxInput->Close();
}


void ElectronCorrection::initializeParameters()
{

  std::cout << "ElectronCorrection: initializing calibration constants..." << std::endl;

  // For now just initialize the parameters to an arbitrary values for
  // debugging. Eventually, this should get the constants from the
  // calibration database.
  m_runGain = 1.0;
  for (int i = 0; i < c_NCDCWires; ++i) {
    m_wireGain[i] = 1.0;
    m_valid[i] = 1.0;
  }

  // get the hadron saturation parameters
  // if the parameters do not exist, use the values in the default constructor
  double hadpar[5];
  std::ifstream parfile("sat-pars.txt");
  if( !parfile.fail() ){
    for( int i = 0; i < 5; ++i ){
      parfile >> hadpar[i];
    }
    parfile.close();
    m_hadsat.setParameters(hadpar);
  }
}


double ElectronCorrection::RunGainCorrection(double& dedx) const
{

  if (m_runGain != 0) {
    double newDedx = dedx / m_runGain;
    return newDedx;
  } else
    return dedx;
}

double ElectronCorrection::WireGainCorrection(int wireID, double& dedx) const
{

  if (m_valid[wireID] && m_wireGain[wireID] != 0) {
    double newDedx = dedx / m_wireGain[wireID];
    return newDedx;
  } else
    return dedx;
}


double ElectronCorrection::HadronCorrection(double costheta, double dedx) const
{

  double newDedx = m_hadsat.D2I(costheta,m_hadsat.I2D(costheta,1.0)*dedx);
  return newDedx;
}


double ElectronCorrection::StandardCorrection(int wireID, double costheta, double dedx) const
{

  double temp = dedx;

  temp = RunGainCorrection(temp);

  temp = WireGainCorrection(wireID, temp);

  //  temp = HadronCorrection(costheta, temp);

  return temp;
}


void ElectronCorrection::calculateMeans(double &mean, double &truncatedMean, double &truncatedMeanErr,
					double dedx[], int size) const
{
  // Calculate the truncated average by skipping the lowest & highest
  // events in the array of dE/dx values
  std::vector<double> sortedDedx;
  for( int i = 0; i < size; ++i )
    sortedDedx.push_back(dedx[i]);
  std::sort(sortedDedx.begin(), sortedDedx.end());

  double truncatedMeanTmp = 0.0;
  double meanTmp = 0.0;
  double sum_of_squares = 0.0;
  int numValuesTrunc = 0;
  const int numDedx = sortedDedx.size();
  const int lowEdgeTrunc = int(numDedx * m_removeLowest);
  const int highEdgeTrunc = int(numDedx * (1 - m_removeHighest));
  for (int i = 0; i < numDedx; i++) {
    std::cout << sortedDedx[i] << "+\t";
    meanTmp += sortedDedx[i];
    if (i >= lowEdgeTrunc and i < highEdgeTrunc) {
      truncatedMeanTmp += sortedDedx[i];
      sum_of_squares += sortedDedx[i] * sortedDedx[i];
      numValuesTrunc++;
    }
  }
  std::cout << std::endl;

  if (numDedx != 0) {
    meanTmp /= numDedx;
  }
  if (numValuesTrunc != 0) {
    truncatedMeanTmp /= numValuesTrunc;
  } else {
    truncatedMeanTmp = meanTmp;
  }

  mean = meanTmp;
  truncatedMean = truncatedMeanTmp;

  if (numValuesTrunc > 1) {
    truncatedMeanErr = sqrt(sum_of_squares / double(numValuesTrunc) - truncatedMeanTmp * truncatedMeanTmp) / double(numValuesTrunc - 1);
  } else {
    truncatedMeanErr = 0;
  }
}


void ElectronCorrection::calculateMeans(double &mean, double &truncatedMean, double &truncatedMeanErr,
					std::vector<double> dedx, int size) const
{
  // Calculate the truncated average by skipping the lowest & highest
  // events in the array of dE/dx values
  std::vector<double> sortedDedx = dedx;
  std::sort(sortedDedx.begin(), sortedDedx.end());

  double truncatedMeanTmp = 0.0;
  double meanTmp = 0.0;
  double sum_of_squares = 0.0;
  int numValuesTrunc = 0;
  const int numDedx = sortedDedx.size();
  const int lowEdgeTrunc = int(numDedx * m_removeLowest);
  const int highEdgeTrunc = int(numDedx * (1 - m_removeHighest));
  for (int i = 0; i < numDedx; i++) {
    std::cout << sortedDedx[i] << "+\t";
    meanTmp += sortedDedx[i];
    if (i >= lowEdgeTrunc and i < highEdgeTrunc) {
      truncatedMeanTmp += sortedDedx[i];
      sum_of_squares += sortedDedx[i] * sortedDedx[i];
      numValuesTrunc++;
    }
  }
  std::cout << std::endl;

  if (numDedx != 0) {
    meanTmp /= numDedx;
  }
  if (numValuesTrunc != 0) {
    truncatedMeanTmp /= numValuesTrunc;
  } else {
    truncatedMeanTmp = meanTmp;
  }

  mean = meanTmp;
  truncatedMean = truncatedMeanTmp;

  if (numValuesTrunc > 1) {
    truncatedMeanErr = sqrt(sum_of_squares / double(numValuesTrunc) - truncatedMeanTmp * truncatedMeanTmp) / double(numValuesTrunc - 1);
  } else {
    truncatedMeanErr = 0;
  }
}
