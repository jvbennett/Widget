#include "ElectronGenerator.h"

ElectronGenerator::ElectronGenerator() : 
  m_nevents(1000000),
  m_upperbg(2.0),
  m_lowerbg(0.2)
{}

ElectronGenerator::ElectronGenerator( int nevents, double upperbg, double lowerbg ) : 
  m_nevents(nevents),
  m_upperbg(upperbg),
  m_lowerbg(lowerbg)
{}

void
ElectronGenerator::generateEvents( TFile* genfile ){

  std::cout << "ElectronGenerator: generating " << m_nevents << " electrons" << std::endl;
  std::cout << "\t\t" << m_upperbg << "\t" << m_lowerbg << std::endl;

  int run = 1;
  int nhit;    // number of hits on this track
  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double q;       // track charge
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  double db;      // the nearest distance to the IP in the xy plane
  double dz;      // the nearest distance to the IP in the z plane
  double chiPi;   // PID chi value
  double eop;     // energy over momentum in the calorimeter

  double dedx_cur, dedx_res;

  // track level variables
  TTree* gentrack = new TTree("track","Fake track sample");

  gentrack->Branch("run",&run,"run/I");
  gentrack->Branch("numGoodHits",&nhit,"numGoodHits/I");
  gentrack->Branch("dedx",&dedxpub,"dedx/D");
  gentrack->Branch("dedxcur",&dedx_cur,"dedxcur/D");
  gentrack->Branch("dedxres",&dedx_res,"dedxres/D");
  gentrack->Branch("dedxsat",&dedx,"dedxsat/D");
  gentrack->Branch("q",&q,"q/D");
  gentrack->Branch("pF",&p,"pF/D");
  gentrack->Branch("bg",&bg,"bg/D");
  gentrack->Branch("costh",&costh,"costh/D");
  gentrack->Branch("db",&db,"db/D");
  gentrack->Branch("dz",&dz,"dz/D");
  gentrack->Branch("chiPi",&chiPi,"chiPi/D");
  gentrack->Branch("eopst",&eop,"eopst/D");

  // hit level variables
  TTree* genhit = new TTree("hit","Fake hit sample");

  const int kMaxHits = 100;
  double doca[kMaxHits];
  double enta[kMaxHits];
  double dedxhit[kMaxHits];

  genhit->Branch("numGoodHits",&nhit,"numGoodHits/I");
  genhit->Branch("doca",doca,"doca[numGoodHits]/D");
  genhit->Branch("enta",enta,"enta[numGoodHits]/D");
  genhit->Branch("dedxhit",dedxhit,"dedxhit[numGoodHits]/D");
  genhit->Branch("dedx",&dedxpub,"dedx/D");
  genhit->Branch("dedxsat",&dedx,"dedxsat/D");
  genhit->Branch("pF",&p,"pF/D");
  genhit->Branch("costh",&costh,"costh/D");
  genhit->Branch("db",&db,"db/D");
  genhit->Branch("dz",&dz,"dz/D");
  genhit->Branch("chiPi",&chiPi,"chiPi/D");

  TH2F* dedx_costh = new TH2F("dedx_costh","dE/dx vs. cos(#theta);cos(#theta);dE/dx",
			      100,-1.0,1.0,100,0.0,10.0);
  TH2F* dedx_pq    = new TH2F("dedx_pq","dE/dx vs. pq;pq  [GeV];dE/dx",
			      100,-1.0*m_upperbg,m_upperbg,100,0.8,1.2);

  TRandom *r0 = new TRandom();
  // set the seed to the current machine clock
  r0->SetSeed(0);

  double mass = Widget::melectron;

  TF1 empcor("empcor","1-0.15*TMath::Exp(-5*x)",0.2,2);
  TF1 docashape("docashape","-78.3*x^2+75.5",-1,1);
  TF1 enta1("enta1","1500*TMath::Landau(x,0.125,0.05)",0,3);
  TF1 enta2("enta2","2000*TMath::Landau(TMath::Abs(x),0.125,0.05)",-3,0);
  TF1 entashape("entashape","enta1+enta2",-3,3);
  TF1 satshape("satshape","1.0-2.0*TMath::Landau(TMath::Abs(x),0.01,0.05)",-1,1);

  // generate the hit level events
  for( int i = 0; i < m_nevents; ++i ){
    if( i % 100 == 0 ){
      std::cout << i << " events" << std::endl;
      run++;
    }

    double q = (r0->Rndm(i) < 0.5) ? -1.0 : 1.0;

    p = r0->Rndm(i)*(m_upperbg-m_lowerbg)+m_lowerbg;
    while( p < m_lowerbg || p > m_upperbg )
      p = r0->Rndm(i)*(m_upperbg-m_lowerbg)+m_lowerbg;
    bg = p/mass;

    costh = 2*(r0->Rndm(i))-1;
    while( costh < -1.0 || costh > 1.0 )
      costh = 2*(r0->Rndm(i))-1;

    dedx = (q < 0) ? (2 - empcor.Eval(p)) : empcor.Eval(p);
    dedx = dedx*satshape.Eval(costh);
    dedx_cur = 1.0;
    dedx_res = 1.05-empcor.Eval(p);

    nhit = r0->Gaus(25);
    while( nhit > 100 || nhit < 0 )
      nhit = r0->Gaus(25);
    eop = 0.5;
    chiPi = r0->Gaus(1.0);

    for( int j = 0; j < nhit; ++j ){
      enta[j] = entashape.GetRandom();
      doca[j] = docashape.GetRandom();
      dedxhit[j] = r0->Landau(dedx,dedx_res);
    }

    dedxpub = dedx;

    dedx_costh->Fill(costh,dedx);
    dedx_pq->Fill(p*q,dedx);

    gentrack->Fill();
    genhit->Fill();
  }
  delete r0;

  genfile->cd();
  gentrack->Write();
  genhit->Write();
  dedx_costh->Write();
  dedx_pq->Write();
  genfile->Close();
}
