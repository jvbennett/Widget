#include "WidgetGenerator.h"

WidgetGenerator::WidgetGenerator() : 
  m_nevents(1000000),
  m_upperbg(5.0),
  m_lowerbg(0.1)
{}

WidgetGenerator::WidgetGenerator( int nevents, double upperbg, double lowerbg ) : 
  m_nevents(nevents),
  m_upperbg(upperbg),
  m_lowerbg(lowerbg)
{}

void
WidgetGenerator::generateEvents( TString filename, double genmass ){

  std::cout << "WidgetGenerator: generating " << m_nevents << " particles of mass " << genmass << std::endl;
  std::cout << m_upperbg << "\t" << m_lowerbg << std::endl;


  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  double db;      // the nearest distance to the IP in the xy plane
  double dz;      // the nearest distance to the IP in the z plane
  double chiPi;   // PID chi value
  double nhit;    // number of hits on this track
  double eop;     // energy over momentum in the calorimeter

  double dedx_cur, dedx_res;

  // do not overwrite a file if it already exists
  TFile* genfile = new TFile(filename,"CREATE");
  TTree* gentree = new TTree("track","Fake track sample");

  gentree->Branch("dedx",&dedxpub,"dedx/D");
  gentree->Branch("dedxcur",&dedx_cur,"dedxcur/D");
  gentree->Branch("dedxres",&dedx_res,"dedxres/D");
  gentree->Branch("dedxsat",&dedx,"dedxsat/D");
  gentree->Branch("pF",&p,"pF/D");
  gentree->Branch("bg",&bg,"bg/D");
  gentree->Branch("costh",&costh,"costh/D");
  gentree->Branch("db",&db,"db/D");
  gentree->Branch("dz",&dz,"dz/D");
  gentree->Branch("chiPi",&chiPi,"chiPi/D");
  gentree->Branch("numGoodHits",&nhit,"numGoodHits/D");
  gentree->Branch("eopst",&eop,"eopst/D");

  TH2F* dedx_costh = new TH2F("dedx_costh","dE/dx vs. cos(#theta);cos(#theta);dE/dx",
			      100,-1.0,1.0,100,0.0,10.0);
  TH2F* dedx_bg    = new TH2F("dedx_bg","dE/dx vs. #beta#gamma;#beta#gamma;dE/dx",
			      100,m_lowerbg,m_upperbg,100,0.0,10.0);

  WidgetParameterization gpar;
  HadronSaturation hadsat;

  TRandom *r0 = new TRandom();
  // set the seed to the current machine clock
  r0->SetSeed(0);

  if( genmass == 0 ){
    std::cout << "WidgetGenerator::ERROR - please input non-zero particle mass" << std::endl;
    return;
  }

  double mass = genmass;
  double massvec[5];
  massvec[0] = Widget::mpion;
  massvec[1] = Widget::mkaon;
  massvec[2] = Widget::mproton;
  massvec[3] = Widget::mmuon;
  massvec[4] = Widget::melectron;

  // generate the events
  for( int i = 0; i < m_nevents; ++i ){

    bg = r0->Rndm(i)*(m_upperbg+m_lowerbg)+m_lowerbg;
    while( bg < m_lowerbg || bg > m_upperbg )
      bg = r0->Rndm(i)*(m_upperbg+m_lowerbg)+m_lowerbg;
    p = bg*mass;

    costh = 2*(r0->Rndm(i))-1;
    while( costh < -1.0 || costh > 1.0 )
      costh = 2*(r0->Rndm(i))-1;

    nhit = r0->Gaus(25);
    eop = 0.5;
    chiPi = r0->Gaus(1.0);

    dedx_cur = gpar.dedxPrediction(bg);
    dedx_res = gpar.sigmaPrediction(bg,nhit,sqrt(1-costh*costh));
    if( mass == Widget::melectron || mass == Widget::mmuon ) dedx_res = 0.05;

    // fake the hadron saturation
    dedx = r0->Gaus(hadsat.I2D(costh,dedx_cur)/hadsat.I2D(costh,1.0),dedx_res/3.0);
    //    if( dedx > 2 ) continue;

    /*    while( dedx > 2 )
	  dedx = r0->Gaus(hadsat.I2D(costh,dedx_cur)/hadsat.I2D(costh,1.0),dedx_res);*/

    dedxpub = dedx;

    dedx_costh->Fill(costh,dedx);
    dedx_bg->Fill(bg,dedx);

    gentree->Fill();
  }
  delete r0;

  genfile->cd();
  gentree->Write();
  dedx_costh->Write();
  dedx_bg->Write();
  genfile->Close();
}

void
WidgetGenerator::simulateDedx( TString infilename, TString outfilename, double genmass ){

  TFile* infile = new TFile(infilename);
  TTree* intree = (TTree*)infile->Get("track");

  double dedx, bg, costh, p;
  intree->SetBranchAddress("p",&p);
  intree->SetBranchAddress("costh",&costh);

  // do not overwrite a file that already exists
  TFile* outfile = new TFile(outfilename,"CREATE");
  TTree* tracks = new TTree("track","Fake track sample");
  tracks->Branch("dedxsat",&dedx,"dedxsat/D");
  tracks->Branch("p",&p,"p/D");
  tracks->Branch("costh",&costh,"costh/D");
    
  TH2F* dedx_costh = new TH2F("dedx_costh","dE/dx vs. cos(#theta);cos(#theta);dE/dx",
			      100,-1.0,1.0,100,0.0,10.0);
  TH2F* dedx_bg    = new TH2F("dedx_bg","dE/dx vs. #beta#gamma;#beta#gamma;dE/dx",
			      100,m_lowerbg,m_upperbg,100,0.0,10.0);

  HadronSaturation hadsat;

  TRandom *r0 = new TRandom();
  // set the seed to the current machine clock
  r0->SetSeed(0);

  double mass = genmass;
  double massvec[5];
  massvec[0] = Widget::mpion;
  massvec[1] = Widget::mkaon;
  massvec[2] = Widget::mproton;
  massvec[3] = Widget::mmuon;
  massvec[4] = Widget::melectron;

  // use the real momentum and angle, but replace the dE/dx measurement
  for( int i = 0; i < intree->GetEntries(); ++i ){
    intree->GetEvent(i);

    if( i < 10 )
      std::cout << "Event " << i << ": p = " << p << ", cos(theta) = " << costh << std::endl;

    if( genmass == 0 ) mass = massvec[i%5];
    bg = p/mass;

    // fake the hadron saturation
    dedx = r0->Gaus(hadsat.I2D(costh,0.8)/hadsat.I2D(costh,1.0),0.05);
    while( dedx > 2 )
      dedx = r0->Gaus(hadsat.I2D(costh,0.8)/hadsat.I2D(costh,1.0),0.05);

    dedx_costh->Fill(costh,dedx);
    dedx_bg->Fill(bg,dedx);
    
    tracks->Fill();
  }
  delete r0;

  outfile->cd();
  tracks->Write();
  dedx_costh->Write();
  dedx_bg->Write();
  outfile->Close();
}

void
WidgetGenerator::simulateReconstruction( TString infilename, TString outfilename ){

  TFile* infile = new TFile(infilename);
  TTree* intree = (TTree*)infile->Get("track");

  int nhits;
  double dedxhit[100], adcraw[100], path[100];
  intree->SetBranchAddress("dedxhit",dedxhit);
  intree->SetBranchAddress("adcraw",adcraw);
  intree->SetBranchAddress("path",path);

  double layer[100], layerdedx[100], layerdx[100];
  intree->SetBranchAddress("layer",layer);
  intree->SetBranchAddress("layerdx",layerdx);
  intree->SetBranchAddress("layerdedx",layerdedx);

  double mean, dedx, dedxerr, costh;
  intree->SetBranchAddress("numLayerHits",&nhits);
  intree->SetBranchAddress("mean",&mean);
  intree->SetBranchAddress("dedx",&dedx);
  intree->SetBranchAddress("dedxerr",&dedxerr);
  intree->SetBranchAddress("costh",&costh);

  // do not overwrite a file that already exists
  TFile* outfile = new TFile(outfilename,"CREATE");
  TTree* tracks = intree->CloneTree(0);

  double newdedxhit[100];
  double newmean, newdedx, newerr;
  tracks->Branch("newmean",&newmean,"newmean/D");
  tracks->Branch("newdedx",&newdedx,"newdedx/D");
  tracks->Branch("newerr",&newerr,"newerr/D");

  for( int i = 0; i < intree->GetEntries(); ++i ){
    intree->GetEntry(i);

    // clear the simulated measurements
    newmean = 0; newdedx = 0; newerr = 0;

    // sort the dedx hits for truncation
    std::vector<double> sortedDedx;
    for( int j = 0; j < nhits; ++j ){
      sortedDedx.push_back(layerdedx[j]);
    }
    std::sort(sortedDedx.begin(), sortedDedx.end());

    // recalculate the mean values
    int lb = int(nhits * 0.05 + 0.5);
    int ub = int(nhits * 0.75 + 0.5);
    double ntrunc = 0;
    for( int j = 0; j < nhits; ++j ){
      newmean += sortedDedx[j];
      if( j >= lb && j < ub ){
	newdedx += sortedDedx[j];
	newerr += sortedDedx[j]*sortedDedx[j];
	ntrunc++;
      }
    }

    if( nhits != 0 ) newmean /= nhits;
    if( ntrunc != 0 ) newdedx /= ntrunc;
    else newdedx = newmean;

    if( ntrunc > 1 ) newerr = sqrt(newerr/ntrunc - newdedx*newdedx)/(ntrunc-1);
    else newerr = 0;

    tracks->Fill();
  }

  tracks->AutoSave();
  outfile->Close();
  infile->Close();
}
