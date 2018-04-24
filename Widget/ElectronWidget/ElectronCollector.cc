#include "ElectronCalibration.h"

ElectronCalibration::ElectronCalibration() :
  m_filename(""),
  m_constfilename(""),
  m_mcFlag(1),
  m_type(0),
  m_fits(0),
  m_docabins(20),
  m_lowerdoca(-1.0),
  m_upperdoca(1.0),
  m_entabins(20),
  m_lowerenta(-1.0),
  m_upperenta(1.0),
  m_costhbins(20)
{}

ElectronCalibration::ElectronCalibration( TString infile, TString constfile, int mcFlag, int type, int fits, int docabins, double upperdoca, double lowerdoca, int entabins, double upperenta, double lowerenta, int costhbins ) :
  m_filename(infile),
  m_constfilename(constfile),
  m_mcFlag(mcFlag),
  m_type(type),
  m_fits(fits),
  m_docabins(docabins),
  m_upperdoca(upperdoca),
  m_lowerdoca(lowerdoca),
  m_entabins(entabins),
  m_upperenta(upperenta),
  m_lowerenta(lowerenta),
  m_costhbins(costhbins)
{}

void
ElectronCalibration::fitRunGains( TFile* outfile ){

  // supress TCanvas::Print messages
  //  gErrorIgnoreLevel = kWarning;

  TFile* infile = new TFile(m_filename);
  TTree* electron = (TTree*)infile->Get("track");
  
  std::cout << "ElectronCalibration: getting run gains for file " << m_filename << std::endl;

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  int run;        // run number per event
  int nhits;      // number of (all) hits on this track
  double dedx;    // dE/dx with electron saturation correction
  double costh;   // cosine of track polar angle

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  electron->SetBranchAddress("run", &run);
  electron->SetBranchAddress("dedxsat", &dedx);
  electron->SetBranchAddress("costh", &costh);

  if( m_type == 0 )
    electron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( m_type == 1 )
    electron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill histograms and fit with Gaussian functions
  // to extract the means and errors
  TTree* satTree = new TTree("electron","dE/dx means and errors");

  int fitrun;      // run number for this bin
  double fitmean;     // mean dE/dx value for this bin
  double fitmeanerr;  // error on mean dE/dx value for this bin
  double fitwidth;    // dE/dx width for this bin
  double fitwidtherr; // error on dE/dx width for this bin

  satTree->Branch("run",&fitrun,"run/I");
  satTree->Branch("dedx",&fitmean,"dedx/D");
  satTree->Branch("dedxerr",&fitmeanerr,"dedxerr/D");
  satTree->Branch("width",&fitwidth,"width/D");
  satTree->Branch("widtherr",&fitwidtherr,"widtherr/D");

  // Create the histograms to be fit (before correction)
  TH1F* hdedx_run = new TH1F("dedx_run","dedx_run",200,0,2);

  // Container to store the constants
  std::map<int, double> m_rungains;

  // Fill the histograms to be fitted
  double lastrun = -1;
  int nentries = electron->GetEntries();
  for( unsigned int index = 0; index < nentries; ++index ){
    electron->GetEvent(index);

    int nhit = 0;
    if( m_type == 0 )
      nhit = std::floor(b3nhit);
    else if( m_type == 1 )
      nhit = b2nhit;

    // clean up bad events and restrict the momentum range
    if( nhit < 0 || nhit > 100 || dedx <= 0 || costh != costh )
      continue;

    if( lastrun == -1 ) lastrun = run;

    if( run == lastrun && index < (nentries-1) ) hdedx_run->Fill(dedx);
    else{
      if( hdedx_run->Integral() < 50 ){
	lastrun = run;
	continue;
      }

      // fill some details for this bin
      fitrun = lastrun;
      // fit the dE/dx distribution in bins of run number
      hdedx_run->Fit("gaus","ql");
      fitmean = hdedx_run->GetFunction("gaus")->GetParameter(1);
      fitmeanerr = hdedx_run->GetFunction("gaus")->GetParError(1);
      fitwidth = hdedx_run->GetFunction("gaus")->GetParameter(2);
      fitwidtherr = hdedx_run->GetFunction("gaus")->GetParError(2);
      satTree->Fill();

      m_rungains[fitrun] = fitmean;

      hdedx_run->Reset();
      hdedx_run->Fill(dedx);

      lastrun = run;
    }
  }

  satTree->Print();
  /*
  outfile->cd();
  TH2F* rungains = new TH2F("rungains","",m_rungains.size(),0,m_rungains.size());
  rungains->SetStats(0);
  satTree->Project("runsigma","width:run");
  TH1F* runsigma = (TH1F*)outfile->Get("runsigma");
  runsigma->SetStats(0);

  TCanvas* can = new TCanvas("electrons","",1200,600);
  can->Divide(2,1);
  can->cd(1);
  rungains->SetTitle(";Run number;dE/dx");
  rungains->Draw("P");

  can->cd(2);
  runsigma->SetTitle(";Run number;dE/dx #sigma");
  runsigma->Draw("P");
  can->SaveAs("plots/rungains.eps");

  delete hdedx_run;
  delete can;
  */

  // write out the data to file
  outfile->cd();
  //  rungains->Write();
  //  runsigma->Write();
  satTree->Write();
  outfile->Close();
  infile->Close();
}

void
ElectronCalibration::plotRunGains( TString filename ){

  std::cout << "ElectronCalibration: plotting run gains for file " << m_filename << std::endl;

  TFile* infile = new TFile(filename);
  TTree* intree = (TTree*)infile->Get("electron");

  intree->Project("rungains","dedx:run");
  TH1F* rungains = (TH1F*)infile->Get("rungains");
  intree->Project("runsigma","dedx:width");
  TH1F* runsigma = (TH1F*)infile->Get("runsigma");

  TCanvas* can = new TCanvas("electrons","",1200,600);
  can->Divide(2,1);
  can->cd(1);
  rungains->SetTitle(";Run number;dE/dx");
  rungains->Draw("P");

  can->cd(2);
  runsigma->SetTitle(";Run number;dE/dx #sigma");
  runsigma->Draw("P");
  can->SaveAs("plots/rungains.eps");

  delete can;

  infile->Close();
}

void
ElectronCalibration::TwoDCorrection( TFile* outfile ){
  
  std::cout << "ElectronCalibration: reading in file " << m_filename << std::endl;

  TFile* infile = new TFile(m_filename);
  TTree* track = (TTree*)infile->Get("track");

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  int kMaxEntries = 100;
  double doca[kMaxEntries];    // distance of closest approach for each hit
  double enta[kMaxEntries];    // CDC cell entrance angle for each hit
  double dedxhit[kMaxEntries]; // dE/dx for each hit
  int nhits;                   // number of hits on this track

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  track->SetBranchAddress("doca", doca);
  track->SetBranchAddress("enta", enta);
  track->SetBranchAddress("dedxhit", dedxhit);

  if( m_type == 0 )
    track->SetBranchAddress("numGoodHits", &b3nhit);
  else if( m_type == 1 )
    track->SetBranchAddress("numGoodLayerHits", &b2nhit);

  double docastep = (m_upperdoca-m_lowerdoca)/m_docabins;
  double entastep[m_entabins];

  // --------------------------------------------------
  // SET THE ENTA BIN SIZES
  // --------------------------------------------------

  // Fill the histograms to be fitted
  TH1F* totalenta = new TH1F("totalenta","",314,0,0.785);
  for( unsigned int index = 0; index < track->GetEntries(); ++index ){
    track->GetEvent(index);

    if( m_type == 0 )
      nhits = std::floor(b3nhit);
    else if( m_type == 1 )
      nhits = b2nhit;
    for( int hit = 0; hit < nhits; ++hit ){

      // use approx. rotational symmetry to map to [-pi/4,pi/4]
      if( enta[hit] > PI/4.0 ) enta[hit] -= PI/2.0;
      else if( enta[hit] < -1.0*PI/4.0 ) enta[hit] += PI/2.0;
      totalenta->Fill(std::abs(enta[hit]));
    }
  }
  entastep[0] = -0.785;
  entastep[39] = 0.785;
  int bin = 20; int binmin = 1;
  double total = totalenta->Integral()*2.0/(m_entabins+1);
  for( int i = 2; i <= 314; ++i ){
    if( totalenta->Integral(binmin,i) >= total ){
      entastep[bin] = 0.0025*i;
      bin++;
      binmin = i+1;
    }
  }
  for( int i = 0; i < m_entabins; ++i ){
    if( i > 0 && i < 20 ) entastep[i] = -1.0*entastep[39-i];
  }

  // --------------------------------------------------
  // OUTPUT TTREE
  // --------------------------------------------------

  TTree* outTree = new TTree("electrons","2D means and errors");

  double satdoca;  // DOCA for this bin
  double satenta;  // entrance angle for this bin
  double satdedx;  // mean dE/dx value for this bin
  double saterror; // error on ^

  outTree->Branch("doca",&satdoca,"doca/D");
  outTree->Branch("enta",&satenta,"enta/D");
  outTree->Branch("dedx",&satdedx,"dedx/D");
  outTree->Branch("error",&saterror,"error/D");

  TH2F* twod = new TH2F("twod","",m_docabins,m_lowerdoca,m_upperdoca,m_entabins,m_lowerenta,m_upperenta);

  // --------------------------------------------------
  // IF DOING TRUNCATION, FILL CONTAINERS AND TRUNCATE
  // --------------------------------------------------

  if( m_fits == 0 ){

    std::vector<double> docaent[m_docabins][m_entabins];
    int docaenthits[m_docabins][m_entabins];
    for( int i = 0; i < m_docabins; ++i ){
      for( int j = 0; j < m_docabins; ++j ){
	docaent[i][j].clear();
	docaenthits[i][j] = 0;
      }
    }

    for( unsigned int index = 0; index < track->GetEntries(); ++index ){
      track->GetEvent(index);

      if( m_type == 0 )
	nhits = std::floor(b3nhit);
      else if( m_type == 1 )
	nhits = b2nhit;
      for( int hit = 0; hit < nhits; ++hit ){

	// use approx. rotational symmetry to map to [-pi/4,pi/4]
	if( enta[hit] > PI/4.0 ) enta[hit] -= PI/2.0;
	else if( enta[hit] < -1.0*PI/4.0 ) enta[hit] += PI/2.0;

	if( doca[hit] > m_upperdoca || doca[hit] < m_lowerdoca ) continue;
	if( enta[hit] > m_upperenta || enta[hit] < m_lowerenta ) continue;

	int docaBin = (int)((doca[hit]-m_lowerdoca)/(m_upperdoca-m_lowerdoca) * m_docabins);
	int entaBin = (int)((enta[hit]-m_lowerenta)/(m_upperenta-m_lowerenta) * m_entabins);

	docaent[docaBin][entaBin].push_back(dedxhit[hit]);
	docaenthits[docaBin][entaBin]++;
      } // end of hit loop
    } // end of event loop

    double tmean;
    ElectronCorrection ecor;
    for( int i = 0; i < m_docabins; ++i ){
      for( int j = 0; j < m_entabins; ++j ){

	// fill some details for this bin
	satdoca = m_lowerdoca+0.5*docastep+i*docastep;
	satenta = m_lowerenta+0.5*entastep[j]+j*entastep[j];

	// calculate the truncated mean and error
	ecor.calculateMeans(tmean,satdedx,saterror,docaent[i][j],docaenthits[i][j]);
	twod->SetBinContent(i,j,satdedx);

	std::cout << i << "\t" << j << "\t" << satdedx << "\t" << saterror << std::endl;

	// fill the output TTree
	outTree->Fill();
      }
    }
  }  

  // --------------------------------------------------
  // IF DOING FITS, LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  if( m_fits == 1 ){

    // Create the histograms to be fit (before correction)
    TH1F* hdoca_enta[m_docabins][m_entabins];

    // initialize the histograms
    for( int i = 0; i < m_docabins; ++i ){
      for( int j = 0; j < m_entabins; ++j ){
	char histname[100];
	sprintf(histname,"doca_%d_enta_%d",i,j);

	hdoca_enta[i][j] = new TH1F(histname,histname,200,0,200);
      }
    }

    // Fill the histograms to be fitted
    for( unsigned int index = 0; index < track->GetEntries(); ++index ){
      track->GetEvent(index);

      if( m_type == 0 )
	nhits = std::floor(b3nhit);
      else if( m_type == 1 )
	nhits = b2nhit;
      for( int hit = 0; hit < nhits; ++hit ){

	// use approx. rotational symmetry to map to [-pi/4,pi/4]
	if( enta[hit] > PI/4.0 ) enta[hit] -= PI/2.0;
	else if( enta[hit] < -1.0*PI/4.0 ) enta[hit] += PI/2.0;

	if( doca[hit] > m_upperdoca || doca[hit] < m_lowerdoca ) continue;
	if( enta[hit] > m_upperenta || enta[hit] < m_lowerenta ) continue;

	int docaBin = (int)((doca[hit]-m_lowerdoca)/(m_upperdoca-m_lowerdoca) * m_docabins);
	int entaBin = (int)((enta[hit]-m_lowerenta)/(m_upperenta-m_lowerenta) * m_entabins);

	hdoca_enta[docaBin][entaBin]->Fill(dedxhit[hit]);
      } // end of hit loop
    } // end of event loop

  
    // --------------------------------------------------
    // FIT IN BINS OF DOCA AND ENTRANCE ANGLE
    // --------------------------------------------------

    // fit the histograms
    int fs1; // make sure none of the fits fail...
    for( int i = 0; i < m_docabins; ++i ){
      for( int j = 0; j < m_entabins; ++j ){

	// fill some details for this bin
	satdoca = m_lowerdoca+0.5*docastep+i*docastep;
	satenta = m_lowerenta+0.5*entastep[j]+j*entastep[j];

	int fs; // check if the fit succeeded
	fs = hdoca_enta[i][j]->Fit("landau","ql");
	double mean = hdoca_enta[i][j]->GetFunction("landau")->GetParameter(1);
	double width = hdoca_enta[i][j]->GetFunction("landau")->GetParameter(2);
	//      fs = hdoca_enta[i][j]->Fit("landau","ql","",mean-2.5*width,mean+2.5*width);
	fs = hdoca_enta[i][j]->Fit("landau","ql");
	if( fs != 0 ) fs1++;

	satdedx = hdoca_enta[i][j]->GetFunction("landau")->GetParameter(1);
	saterror = hdoca_enta[i][j]->GetFunction("landau")->GetParError(1);

	if( fs == 0 ) twod->SetBinContent(i,j,satdedx);

	// fill the tree for this bin
	outTree->Fill();
      }
    }
    if( fs1 != 0 ) std::cout << "\t\t" <<"MEAN FITS FAILED: " << fs1 << std::endl;


    // Print the histograms for quality control
    TCanvas* ctmp = new TCanvas("tmp","tmp",900,900);
    ctmp->Divide(3,3);
    std::stringstream psname; psname << "plots/twoDFits.ps[";
    ctmp->Print(psname.str().c_str());
    psname.str(""); psname << "plots/twoDFits.ps";
    for( int i = 0 ; i < m_docabins; ++i ){
      for( int j = 0; j < m_entabins; ++j ){
	ctmp->cd(j%9+1);
	hdoca_enta[i][j]->Draw();
	if((j+1)%9==0)
	  ctmp->Print(psname.str().c_str());
      }
    }
    psname.str(""); psname << "plots/twoDFits.ps]";
    ctmp->Print(psname.str().c_str());
    delete ctmp;

    for( int i = 0; i < m_docabins; ++i ){
      for( int j = 0; j < m_entabins; ++j ){
	delete hdoca_enta[i][j];
      }
    }
    delete totalenta;
  }

  // write out the data to file
  outfile->cd();
  outTree->Write();
  twod->Write();
  infile->Close();
}

void
ElectronCalibration::SaturationCorrection( TFile* outfile ){

  TFile* infile = new TFile(m_filename);
  TTree* electron = (TTree*)infile->Get("track");
  
  std::cout << "ElectronCalibration: electron saturation correction " << m_filename << std::endl;

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  int run;     // run number per event
  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double costh;   // cosine of track polar angle

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  electron->SetBranchAddress("run", &run);
  electron->SetBranchAddress("dedx", &dedxpub);
  electron->SetBranchAddress("dedxsat", &dedx);
  electron->SetBranchAddress("costh", &costh);

  if( m_type == 0 )
    electron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( m_type == 1 )
    electron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill histograms and fit with Gaussian functions
  // to extract the means and errors
  TTree* satTree = new TTree("electron","saturation correction");

  double costhbin;    // center of cos(theta) bin
  double fitmean;     // mean dE/dx value for this bin
  double fitmeanerr;  // error on mean dE/dx value for this bin
  double fitwidth;    // dE/dx width for this bin
  double fitwidtherr; // error on dE/dx width for this bin

  satTree->Branch("costh",&costhbin,"costh/D");
  satTree->Branch("dedx",&fitmean,"dedx/D");
  satTree->Branch("dedxerr",&fitmeanerr,"dedxerr/D");
  satTree->Branch("width",&fitwidth,"width/D");
  satTree->Branch("widtherr",&fitwidtherr,"widtherr/D");

  // Create the histograms to be fit (before correction)
  TH1F* hdedx_costh[m_costhbins];
  for( int i = 0; i < m_costhbins; ++i ){
    hdedx_costh[i] = new TH1F(TString::Format("dedx_costh%i",i),";dE/dx;cos(#theta)",200,0,2);
  }

  // Fill the histograms to be fitted
  for( unsigned int index = 0; index < electron->GetEntries(); ++index ){
    electron->GetEvent(index);

    int nhit = 0;
    if( m_type == 0 )
      nhit = std::floor(b3nhit);
    else if( m_type == 1 )
      nhit = b2nhit;

    // clean up bad events and restrict the momentum range
    if( nhit < 0 || nhit > 100 || dedx <= 0 || costh != costh )
      continue;

    int bin = floor((costh+1.0)/2*m_costhbins);
    hdedx_costh[bin]->Fill(dedxpub);
  }

  // Fit the histograms
  for( unsigned int i = 0; i < m_costhbins; ++i ){
    costhbin = (2.0*i+1)/m_costhbins-1;
    hdedx_costh[i]->Fit("gaus","ql");
    fitmean = hdedx_costh[i]->GetFunction("gaus")->GetParameter(1);
    fitmeanerr = hdedx_costh[i]->GetFunction("gaus")->GetParError(1);
    fitwidth = hdedx_costh[i]->GetFunction("gaus")->GetParameter(2);
    fitwidtherr = hdedx_costh[i]->GetFunction("gaus")->GetParError(2);
    satTree->Fill();
  }

  outfile->cd();
  satTree->Project("satmean","dedx:costh");
  TH1F* satmean = (TH1F*)outfile->Get("satmean");
  satmean->SetStats(0);
  satTree->Project("satsigma","width:costh");
  TH1F* satsigma = (TH1F*)outfile->Get("satsigma");
  satsigma->SetStats(0);

  TCanvas* can = new TCanvas("electrons","",1200,600);
  can->Divide(2,1);
  can->cd(1);
  satmean->SetTitle(";cos(#theta);dE/dx");
  satmean->Draw("P");

  can->cd(2);
  satsigma->SetTitle(";cos(#theta);dE/dx sigma");
  satsigma->Draw("P");
  can->SaveAs("plots/satmean.eps");

  for( int i = 0; i < m_costhbins; ++i ){
    delete hdedx_costh[i];
  }
  delete can;

  // write out the data to file
  outfile->cd();
  satmean->Write();
  satsigma->Write();
  satTree->Write();
  outfile->Close();

  infile->Close();
}
