#include "HadronCalibration.h"

HadronCalibration::HadronCalibration() {}

void
HadronCalibration::fitBGCurve( std::vector< TString > particles, TString filename, std::string paramfile ){

  // read in a file that contains fit results for bg bins
  TFile* infile = new TFile(filename);

  // multigraphs to hold the curve and residual results
  TMultiGraph* gr = new TMultiGraph("gr",";#beta#gamma;dE/dx");
  TMultiGraph* gr_res = new TMultiGraph("gr_res",";dE/dx;#sigma");
  TMultiGraph* grp = new TMultiGraph("grp",";Momentum;dE/dx");
  TMultiGraph* grp_res = new TMultiGraph("grp_res",";Momentum;#sigma");

  // multigraphs to hold the chi and sigma distributions
  TMultiGraph* grchi = new TMultiGraph("grchi",";#beta#gamma;#chi mean");
  TMultiGraph* grsigma = new TMultiGraph("grsigma",";#beta#gamma;#sigma");
  TMultiGraph* grchip = new TMultiGraph("grchip",";Momentum;#chi mean");
  TMultiGraph* grsigmap = new TMultiGraph("grsigmap",";Momentum;#sigma");

  const int npart = 5;
  TGraphErrors particle_dedx[npart];
  TGraph particle_res[npart];
  TGraphErrors particle_dedxp[npart];
  TGraph particle_resp[npart];
  TGraph newchi[npart];
  TGraph newsigma[npart];
  TGraph newchip[npart];
  TGraph newsigmap[npart];

  // --------------------------------------------------
  // FILL BG CURVE VALUES
  // --------------------------------------------------

  double mass = 0.0;
  for( int i = 0; i < particles.size(); ++i ){
    TString particle = particles[i];

    if( particle == "pion" ) mass = Widget::mpion;
    else if( particle == "kaon" ) mass = Widget::mkaon;
    else if( particle == "proton" ) mass = Widget::mproton;
    else if( particle == "muon" ) mass = Widget::mmuon;
    else if( particle == "electron" ) mass = Widget::melectron;
    if( mass == 0.0 ) exit(1);

    TTree* hadron = (TTree*)infile->Get(particle);
  
    std::cout << "HadronCalibration: reading " << particle << " in file " << filename << std::endl;

    double dedx_pub;  // dE/dx without electron saturation correction
    double dedx;      // dE/dx with electron saturation correction
    double bg;    // track beta-gamma

    double dedxerr;   // 
    double dedx_res;  // 
    //    double nhit_avg;  // 
    double sin_avg;   // 
    double chi;       // 
    double chi_pub;   // 
    double sigma;     // 
    double sigma_pub; // 

    hadron->SetBranchAddress("dedx",&dedx);
    hadron->SetBranchAddress("dedx_pub",&dedx_pub);
    hadron->SetBranchAddress("bg_avg",&bg);
    hadron->SetBranchAddress("error",&dedxerr);
    hadron->SetBranchAddress("dedxres_avg",&dedx_res);
    //    hadron->SetBranchAddress("nhit_avg",&nhit_avg);
    hadron->SetBranchAddress("sinth_avg",&sin_avg);
    hadron->SetBranchAddress("chi",&chi);
    hadron->SetBranchAddress("chi_pub",&chi_pub);
    hadron->SetBranchAddress("sigma",&sigma);
    hadron->SetBranchAddress("sigma_pub",&sigma_pub);

    for( int j = 0; j < hadron->GetEntries(); ++j ){
      hadron->GetEvent(j);

      particle_dedx[i].SetPoint(j,bg,dedx);
      particle_dedx[i].SetPointError(j,0,dedxerr);
      particle_res[i].SetPoint(j,dedx,dedx_res);

      particle_dedxp[i].SetPoint(j,bg*mass,dedx);
      particle_dedxp[i].SetPointError(j,0,dedxerr);
      particle_resp[i].SetPoint(j,bg*mass,dedx_res);

      newchi[i].SetPoint(j,bg,chi);
      newsigma[i].SetPoint(j,bg,sigma);
      newchip[i].SetPoint(j,bg*mass,chi);
      newsigmap[i].SetPoint(j,bg*mass,sigma);
    }

    newchi[i].SetMarkerStyle(4);
    newchi[i].SetMarkerColor(i+1);
    if( i == 4 ) newchi[i].SetMarkerColor(i+2);
    newchip[i].SetMarkerStyle(4);
    newchip[i].SetMarkerColor(i+1);
    if( i == 4 ) newchip[i].SetMarkerColor(i+2);

    newsigma[i].SetMarkerStyle(4);
    newsigma[i].SetMarkerColor(i+1);
    if( i == 4 ) newsigma[i].SetMarkerColor(i+2);
    newsigmap[i].SetMarkerStyle(4);
    newsigmap[i].SetMarkerColor(i+1);
    if( i == 4 ) newsigmap[i].SetMarkerColor(i+2);

    particle_dedx[i].SetMarkerStyle(4);
    particle_dedx[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_dedx[i].SetMarkerColor(i+2);
    particle_dedxp[i].SetMarkerStyle(4);
    particle_dedxp[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_dedxp[i].SetMarkerColor(i+2);

    particle_res[i].SetMarkerStyle(4);
    particle_res[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_res[i].SetMarkerColor(i+2);
    particle_resp[i].SetMarkerStyle(4);
    particle_resp[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_resp[i].SetMarkerColor(i+2);

    gr->Add(&particle_dedx[i]);
    gr_res->Add(&particle_res[i]);
    grp->Add(&particle_dedxp[i]);
    grp_res->Add(&particle_resp[i]);
    grchi->Add(&newchi[i]);
    grsigma->Add(&newsigma[i]);
    grchip->Add(&newchip[i]);
    grsigmap->Add(&newsigmap[i]);
  }


  // --------------------------------------------------
  // FIT BG CURVE
  // --------------------------------------------------

  WidgetParameterization gpar(paramfile);
  WidgetCurve* gc = new WidgetCurve();

  double bgmin1 = 0.1, bgmax1 = 5;
  double bgmin2 = 4, bgmax2 = 12;
  double bgmin3 = 8, bgmax3 = 2000;

 TF1* fdedx1 = new TF1("fdedx1",gc,bgmin1,bgmax1,9,"WidgetCurve");
  fdedx1->SetParameter(0,1);
  for( int i = 1; i < 9; ++i ){
    fdedx1->SetParameter(i,gpar.getCurvePars(i-1));
  }
  // fix one of the redundant parameters
  fdedx1->FixParameter(4,1);

  TF1* fdedx2 = new TF1("fdedx2",gc,bgmin2,bgmax2,5,"WidgetCurve");
  fdedx2->SetParameter(0,2);
  for( int i = 1; i < 5; ++i ){
    fdedx2->SetParameter(i,gpar.getCurvePars(7+i));
  }

  TF1* fdedx3 = new TF1("fdedx3",gc,bgmin3,bgmax3,5,"WidgetCurve");
  fdedx3->SetParameter(0,3);
  for( int i = 1; i < 5; ++i ){
    fdedx3->SetParameter(i,gpar.getCurvePars(11+i));
  }

  TF1* fdedx4 = new TF1("fdedx4","1",100,100000);

  TCanvas* bandcan = new TCanvas("bandcan","dE/dx",820,750);
  grp->Draw("APE");
  bandcan->SaveAs("plots/dedxbands.eps");
  delete bandcan;

  TCanvas* logbandcan = new TCanvas("logbandcan","dE/dx",820,750);
  logbandcan->cd()->SetLogy();
  grp->Draw("APE");
  logbandcan->SaveAs("plots/dedxbands_log.eps");
  delete logbandcan;

  TCanvas* bgcurvecan = new TCanvas("bgcurvecan","dE/dx",820,750);
  bgcurvecan->cd()->SetLogy();
  bgcurvecan->cd()->SetLogx();
  gr->Draw("APE");

  int stat1 = gr->Fit("fdedx1","","",bgmin1,bgmax1);
  for( int i = 0; i < 50; ++i ){
    if( stat1 == 0 ) break;
    stat1 = gr->Fit("fdedx1","","",bgmin1,bgmax1);
  }
  int stat2 = gr->Fit("fdedx2","","",bgmin2,bgmax2);
  for( int i = 0; i < 50; ++i ){
    if( stat2 == 0 ) break;
    stat2 = gr->Fit("fdedx2","","",bgmin2,bgmax2);
  }
  int stat3 = gr->Fit("fdedx3","","",bgmin3,bgmax3);
  for( int i = 0; i < 50; ++i ){
    if( stat3 == 0 ) break;
    stat3 = gr->Fit("fdedx3","","",bgmin3,bgmax3);
  }

  // if the fit was successful, write out the updated parameters
  if( stat1 != 0 ) std::cout << "WARNING: BG FIT 1 FAILED..." << std::endl;
  for( int i = 1; i < 9; ++i ){
    gpar.setCurvePars(i-1,fdedx1->GetParameter(i));
  }
  if( stat2 != 0 ) std::cout << "WARNING: BG FIT 2 FAILED..." << std::endl;
  for( int i = 1; i < 5; ++i ){
    gpar.setCurvePars(7+i,fdedx2->GetParameter(i));
  }
  if( stat3 != 0 ) std::cout << "WARNING: BG FIT 3 FAILED..." << std::endl;
  for( int i = 1; i < 5; ++i ){
    gpar.setCurvePars(11+i,fdedx3->GetParameter(i));
  }
  
  fdedx1->SetLineColor(kBlack);
  fdedx1->Draw("same");
  fdedx2->SetLineColor(kBlue);
  fdedx2->Draw("same");
  fdedx3->SetLineColor(kRed);
  fdedx3->Draw("same");
  fdedx4->SetLineColor(kGray);
  fdedx4->Draw("same");

  TLegend* tleg = new TLegend(0.4,0.6,0.6,0.8);
  for( int i = 0; i < 5; ++i ){
    tleg->AddEntry(&particle_dedx[i],particles[i],"p");
  }
  tleg->Draw("same");
  
  bgcurvecan->SaveAs("plots/bgcurve.eps");
  delete bgcurvecan;


  // --------------------------------------------------
  // GET RESIDUALS AND CHIS
  // --------------------------------------------------

  double A = 4.5, B = 10;

  TMultiGraph* fit_res = new TMultiGraph("fit_res",";#beta#gamma;Residual");
  TGraph particle_residual[npart];
  int respoint = 1;
  double rmin = 1.0, rmax = 1.0;
  for( int i = 0; i < npart; ++i ){
    for( int j = 0; j < particle_dedx[i].GetN(); ++j ){
      double x, y, fit;
      particle_dedx[i].GetPoint(j,x,y);
      if( y == 0 ) continue;
      if( x < A ) 
	fit = fdedx1->Eval(x);
      else if( x < B )
	fit = fdedx2->Eval(x);
      else
	fit = fdedx3->Eval(x);

      // the curve is just 1 for electrons...
      if( npart == 4 ) fit = 1.0;

      particle_residual[i].SetPoint(respoint++,x,fit/y);
      if( fit/y < rmin ) rmin = fit/y;
      else if( fit/y > rmax ) rmax = fit/y;
    }

    particle_residual[i].SetMarkerStyle(4);
    particle_residual[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_residual[i].SetMarkerColor(i+2);
    fit_res->Add(&particle_residual[i]);
  }
  fit_res->SetMinimum(rmin*0.98);
  fit_res->SetMaximum(rmax*1.02);

  TCanvas* bgrescan = new TCanvas("bgrescan","dE/dx",820,750);
  bgrescan->cd()->SetLogx();
  fit_res->Draw("AP");
  tleg->Draw("same");

  bgrescan->SaveAs("plots/bgresidual.eps");
  delete bgrescan;


  // --------------------------------------------------
  // FIT SIGMA VS DEDX
  // --------------------------------------------------

  TF1* sigvsdedx = new TF1("sigvsdedx","[0]+[1]*x",0.0,10.0);
  sigvsdedx->SetParameter(0,gpar.getDedxPars(0));
  sigvsdedx->SetParameter(1,gpar.getDedxPars(1));

  TCanvas* sigcan = new TCanvas("sigcan","dE/dx",820,750);
  sigcan->cd();
  gr_res->Draw("APE");
  tleg->Draw("same");

  int status = gr_res->Fit("sigvsdedx","qm","",0.0,4.0);
  for( int i = 0; i < 10; ++i ){
    if( status == 0 ) break;
    status = gr_res->Fit("sigvsdedx","q","",0.0,4.0);
  }
  sigcan->SaveAs("plots/sigmavsdedx.eps");  
  delete sigcan;

  if( status != 0 ) std::cout << "WARNING: SIGMA VS DEDX FIT FAILED..." << std::endl;
  gpar.setDedxPars(0,sigvsdedx->GetParameter(0));
  gpar.setDedxPars(1,sigvsdedx->GetParameter(1));

  // write out the (possibly) updated parameters to file
  gpar.printParameters("parameters.bgcurve.fit");

  sigvsdedx->Draw("same");
  tleg->Draw("same");


  // --------------------------------------------------
  // PLOT CHI MEAN VS BETA-GAMMA
  // --------------------------------------------------

  TLine* line0 = new TLine(0,0,10000,0);
  line0->SetLineStyle(kDashed);
  line0->SetLineColor(kRed);
  TLine* line1 = new TLine(0,1,10000,1);
  line1->SetLineStyle(kDashed);
  line1->SetLineColor(kRed);

  TCanvas* cbgcan = new TCanvas("cbgcan","Mean #chi vs. #beta#gamma");
  grchi->SetMinimum(-0.5);
  grchi->SetMaximum(+0.5);
  grchi->Draw("AP");
  line0->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/chimeanVSbetagamma.eps");

  grchi->GetXaxis()->SetLimits(0,25);
  grchi->Draw("AP");
  line0->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/chimeanVSbetagamma_zoom.eps");

  //  grsigma->GetXaxis()->SetTitle("#beta#gamma");
  //  grsigma->GetYaxis()->SetTitle("#chi mean");
  grsigma->SetMinimum(0.5);
  grsigma->SetMaximum(1.5);
  grsigma->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/sigmaVSbetagamma.eps");

  grsigma->GetXaxis()->SetLimits(0,25);
  grsigma->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/sigmaVSbetagamma_zoom.eps");

  TCanvas* cpcan = new TCanvas("cpcan","Mean #chi vs. momentum");
  grchip->SetMinimum(-0.5);
  grchip->SetMaximum(+0.5);
  grchip->Draw("AP");
  line0->Draw("same");
  tleg->Draw("same");
  cpcan->SaveAs("plots/chimeanVSmomentum.eps");

  grsigmap->SetMinimum(0.5);
  grsigmap->SetMaximum(1.5);
  grsigmap->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cpcan->SaveAs("plots/sigmaVSmomentum.eps");

  delete cbgcan;
  delete cpcan;
  delete line0;
  delete line1;
  delete gr;
  delete gr_res;
  delete grp;
  delete grp_res;
  delete grchi;
  delete grsigma;
  delete grchip;
  delete grsigmap;
}

void
HadronCalibration::plotBGCurve( std::vector< TString > particles, TString filename, std::string paramfile ){

  // read in a file that contains fit results for bg bins
  TFile* infile = new TFile(filename);

  // multigraphs to hold the curve and residual results
  TMultiGraph* gr = new TMultiGraph("gr",";#beta#gamma;dE/dx");
  TMultiGraph* gr_res = new TMultiGraph("gr_res",";dE/dx;#sigma");
  TMultiGraph* grp = new TMultiGraph("grp",";Momentum;dE/dx");
  TMultiGraph* grp_res = new TMultiGraph("grp_res",";Momentum;#sigma");

  // multigraphs to hold the chi and sigma distributions
  TMultiGraph* grchi = new TMultiGraph("grchi",";#beta#gamma;#chi mean");
  TMultiGraph* grsigma = new TMultiGraph("grsigma",";#beta#gamma;#sigma");
  TMultiGraph* grchip = new TMultiGraph("grchip",";Momentum;#chi mean");
  TMultiGraph* grsigmap = new TMultiGraph("grsigmap",";Momentum;#sigma");

  const int npart = 5;
  TGraphErrors particle_dedx[npart];
  TGraph particle_res[npart];
  TGraphErrors particle_dedxp[npart];
  TGraph particle_resp[npart];
  TGraph newchi[npart];
  TGraph newsigma[npart];
  TGraph newchip[npart];
  TGraph newsigmap[npart];

  // --------------------------------------------------
  // FILL BG CURVE VALUES
  // --------------------------------------------------

  double mass = 0.0;
  for( int i = 0; i < particles.size(); ++i ){
    TString particle = particles[i];

    if( particle == "pion" ) mass = Widget::mpion;
    else if( particle == "kaon" ) mass = Widget::mkaon;
    else if( particle == "proton" ) mass = Widget::mproton;
    else if( particle == "muon" ) mass = Widget::mmuon;
    else if( particle == "electron" ) mass = Widget::melectron;
    if( mass == 0.0 ) exit(1);

    TTree* hadron = (TTree*)infile->Get(particle);
  
    std::cout << "HadronCalibration: reading " << particle << " in file " << filename << std::endl;

    double dedx_pub;  // dE/dx without electron saturation correction
    double dedx;      // dE/dx with electron saturation correction
    double bg;    // track beta-gamma

    double dedxerr;   // 
    double dedx_res;  // 
    //    double nhit_avg;  // 
    double sin_avg;   // 
    double chi;       // 
    double chi_pub;   // 
    double sigma;     // 
    double sigma_pub; // 

    hadron->SetBranchAddress("dedx",&dedx);
    hadron->SetBranchAddress("dedx_pub",&dedx_pub);
    hadron->SetBranchAddress("bg_avg",&bg);
    hadron->SetBranchAddress("error",&dedxerr);
    hadron->SetBranchAddress("dedxres_avg",&dedx_res);
    //    hadron->SetBranchAddress("nhit_avg",&nhit_avg);
    hadron->SetBranchAddress("sinth_avg",&sin_avg);
    hadron->SetBranchAddress("chi",&chi);
    hadron->SetBranchAddress("chi_pub",&chi_pub);
    hadron->SetBranchAddress("sigma",&sigma);
    hadron->SetBranchAddress("sigma_pub",&sigma_pub);

    for( int j = 0; j < hadron->GetEntries(); ++j ){
      hadron->GetEvent(j);

      particle_dedx[i].SetPoint(j,bg,dedx);
      particle_dedx[i].SetPointError(j,0,dedxerr);
      particle_res[i].SetPoint(j,dedx,dedx_res);

      particle_dedxp[i].SetPoint(j,bg*mass,dedx);
      particle_dedxp[i].SetPointError(j,0,dedxerr);
      particle_resp[i].SetPoint(j,bg*mass,dedx_res);

      newchi[i].SetPoint(j,bg,chi);
      newsigma[i].SetPoint(j,bg,sigma);
      newchip[i].SetPoint(j,bg*mass,chi);
      newsigmap[i].SetPoint(j,bg*mass,sigma);
    }

    newchi[i].SetMarkerStyle(4);
    newchi[i].SetMarkerColor(i+1);
    if( i == 4 ) newchi[i].SetMarkerColor(i+2);
    newchip[i].SetMarkerStyle(4);
    newchip[i].SetMarkerColor(i+1);
    if( i == 4 ) newchip[i].SetMarkerColor(i+2);

    newsigma[i].SetMarkerStyle(4);
    newsigma[i].SetMarkerColor(i+1);
    if( i == 4 ) newsigma[i].SetMarkerColor(i+2);
    newsigmap[i].SetMarkerStyle(4);
    newsigmap[i].SetMarkerColor(i+1);
    if( i == 4 ) newsigmap[i].SetMarkerColor(i+2);

    particle_dedx[i].SetMarkerStyle(4);
    particle_dedx[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_dedx[i].SetMarkerColor(i+2);
    particle_dedxp[i].SetMarkerStyle(4);
    particle_dedxp[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_dedxp[i].SetMarkerColor(i+2);

    particle_res[i].SetMarkerStyle(4);
    particle_res[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_res[i].SetMarkerColor(i+2);
    particle_resp[i].SetMarkerStyle(4);
    particle_resp[i].SetMarkerColor(i+1);
    if( i == 4 ) particle_resp[i].SetMarkerColor(i+2);

    gr->Add(&particle_dedx[i]);
    gr_res->Add(&particle_res[i]);
    grp->Add(&particle_dedxp[i]);
    grp_res->Add(&particle_resp[i]);
    grchi->Add(&newchi[i]);
    grsigma->Add(&newsigma[i]);
    grchip->Add(&newchip[i]);
    grsigmap->Add(&newsigmap[i]);
  }

  TCanvas* bandcan = new TCanvas("bandcan","dE/dx",820,750);
  grp->Draw("APE");
  bandcan->SaveAs("plots/dedxbands.eps");
  delete bandcan;

  TCanvas* logbandcan = new TCanvas("logbandcan","dE/dx",820,750);
  logbandcan->cd()->SetLogy();
  grp->Draw("APE");
  logbandcan->SaveAs("plots/dedxbands_log.eps");
  delete logbandcan;

  TCanvas* bgcurvecan = new TCanvas("bgcurvecan","dE/dx",820,750);
  bgcurvecan->cd()->SetLogy();
  bgcurvecan->cd()->SetLogx();
  gr->Draw("APE");

  // --------------------------------------------------
  // PLOT CHI MEAN VS BETA-GAMMA
  // --------------------------------------------------

  TLine* line0 = new TLine(0,0,10000,0);
  line0->SetLineStyle(kDashed);
  line0->SetLineColor(kRed);
  TLine* line1 = new TLine(0,1,10000,1);
  line1->SetLineStyle(kDashed);
  line1->SetLineColor(kRed);

  TCanvas* cbgcan = new TCanvas("cbgcan","Mean #chi vs. #beta#gamma");
  grchi->SetMinimum(-0.5);
  grchi->SetMaximum(+0.5);
  grchi->Draw("AP");
  line0->Draw("same");
  TLegend* tleg = new TLegend(0.4,0.6,0.6,0.8);
  for( int i = 0; i < 5; ++i ){
    tleg->AddEntry(&particle_dedx[i],particles[i],"p");
  }
  tleg->Draw("same");
  cbgcan->SaveAs("plots/chimeanVSbetagamma.eps");

  grchi->GetXaxis()->SetLimits(0,25);
  grchi->Draw("AP");
  line0->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/chimeanVSbetagamma_zoom.eps");

  //  grsigma->GetXaxis()->SetTitle("#beta#gamma");
  //  grsigma->GetYaxis()->SetTitle("#chi mean");
  grsigma->SetMinimum(0.5);
  grsigma->SetMaximum(1.5);
  grsigma->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/sigmaVSbetagamma.eps");

  grsigma->GetXaxis()->SetLimits(0,25);
  grsigma->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cbgcan->SaveAs("plots/sigmaVSbetagamma_zoom.eps");

  TCanvas* cpcan = new TCanvas("cpcan","Mean #chi vs. momentum");
  grchip->SetMinimum(-0.5);
  grchip->SetMaximum(+0.5);
  grchip->Draw("AP");
  line0->Draw("same");
  tleg->Draw("same");
  cpcan->SaveAs("plots/chimeanVSmomentum.eps");

  grsigmap->SetMinimum(0.5);
  grsigmap->SetMaximum(1.5);
  grsigmap->Draw("AP");
  line1->Draw("same");
  tleg->Draw("same");
  cpcan->SaveAs("plots/sigmaVSmomentum.eps");

  TFile* outfile = new TFile("dedxbands.root","RECREATE");
  outfile->cd();
  gr->Write();
  gr_res->Write();
  grp->Write();
  grp_res->Write();
  grchi->Write();
  grsigma->Write();
  grchip->Write();
  grsigmap->Write();
  outfile->Close();

  delete cbgcan;
  delete cpcan;
  delete line0;
  delete line1;
  delete gr;
  delete gr_res;
  delete grp;
  delete grp_res;
  delete grchi;
  delete grsigma;
  delete grchip;
  delete grsigmap;
}

void
HadronCalibration::fitSigmaVsNHit( TString filename, std::string paramfile, int mcflag = 0, int type = 0 ){

  TFile* infile = new TFile(filename);
  TTree* hadron = (TTree*)infile->Get("track");
  
  std::cout << "HadronCalibration: reading in file " << filename << std::endl;

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  double db;      // the nearest distance to the IP in the xy plane
  double dz;      // the nearest distance to the IP in the z plane
  double chiPi;   // PID chi value
  double eop;     // energy over momentum in the calorimeter

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  hadron->SetBranchAddress("dedx", &dedxpub);
  hadron->SetBranchAddress("dedxsat", &dedx);
  hadron->SetBranchAddress("pF", &p);
  hadron->SetBranchAddress("costh", &costh);
  hadron->SetBranchAddress("db", &db);
  hadron->SetBranchAddress("dz", &dz);
  hadron->SetBranchAddress("chiPi", &chiPi);
  hadron->SetBranchAddress("eopst", &eop);

  if( type == 0 )
    hadron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( type == 1 )
    hadron->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    hadron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  const int m_nhitbins = 20;
  const double m_uppernhit = 32.5, m_lowernhit = 5.5;

  double nhitstep = (m_uppernhit-m_lowernhit)/m_nhitbins;

  // Create the histograms to be fit
  TH1F* hdedx_nhit[m_nhitbins];
  TH1F* hdedx_nhit_pub[m_nhitbins];

  // initialize the histograms
  for( int i = 0; i < m_nhitbins; ++i ){
    char histname[100],histname2[100],histname3[100],histname4[100],histname5[100],histname6[100];
    sprintf(histname, "dedx_nhit_%d",i);
    sprintf(histname4, "dedxpub_nhit_%d",i);

    hdedx_nhit[i] = new TH1F(histname,histname,200,-1,1);
    hdedx_nhit_pub[i] = new TH1F(histname4,histname4,200,-1,1);
  }

  // Create some containers to calculate averages
  double sumnhit[m_nhitbins];
  double sumbg[m_nhitbins];
  double sumsin[m_nhitbins];
  int sumsize[m_nhitbins];
  for( int i = 0; i < m_nhitbins; ++i ){
    sumnhit[i] = 0;
    sumbg[i] = 0;
    sumsin[i] = 0;
    sumsize[i] = 0;
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
  WidgetParameterization gpar(paramfile);

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill the histograms to be fitted
  for( unsigned int index = 0; index < hadron->GetEntries(); ++index ){
    hadron->GetEvent(index);

    bg = fabs(p)/Widget::mpion;

    int nhit;
    if( type == 0 )
      nhit = std::floor(b3nhit);
    else if( type == 1 )
      nhit = b2nhit;

    if( dedx <= 0 || nhit <= m_lowernhit || nhit >= m_uppernhit )
      continue;

    int nhitBin = (int)((nhit-m_lowernhit)/(m_uppernhit-m_lowernhit) * m_nhitbins);

    double dedx_new = dedx;
    if( !mcflag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_cur = gpar.dedxPrediction(bg);
    double dedx_res = gpar.sigmaPrediction(dedx,nhit,sqrt(1-costh*costh));
    double chi_new  = (dedx_new-dedx_cur)/dedx_res;

    double res_cor = gpar.sinPrediction(sqrt(1-costh*costh));
    int i = (int) ( (fabs(nhit)-m_lowernhit)/nhitstep );
    hdedx_nhit[i]->Fill((dedx_new-dedx_cur)/res_cor);
    hdedx_nhit_pub[i]->Fill((dedx_new-dedx_cur)/res_cor);

    if( fabs(nhit-(m_lowernhit+0.5*nhitstep+nhitBin*nhitstep)) < 0.5*nhitstep ){
      sumnhit[nhitBin] += nhit;
      sumbg[nhitBin] += bg;
      sumsin[nhitBin] += sqrt(1-costh*costh);
      sumsize[nhitBin] += 1;
    }
  }// end of event loop


  // --------------------------------------------------
  // FIT IN BINS OF NHIT
  // --------------------------------------------------

  // fit the histograms with Gaussian functions
  // and extract the means and errors
  TTree* nhitTree = new TTree("sigmavsnhit","dE/dx means and errors");

  double nhitdedx;          // dE/dx mean value for this bin
  double nhitsigma;         // dE/dx width value for this bin
  double nhitdedxpub;       // dE/dx mean value for this bin
  double nhitsigmapub;      // dE/dx width value for this bin
  double nhitavg;           // average nhit value for this bin
  double nhitbgavg;         // average bg value for this bin
  double nhitsinavg;        // average sin(theta) value for this bin

  nhitTree->Branch("dedx",&nhitdedx,"dedx/D");
  nhitTree->Branch("width",&nhitsigma,"width/D");
  nhitTree->Branch("dedx_pub",&nhitdedxpub,"dedx_pub/D");
  nhitTree->Branch("sigma_pub",&nhitsigmapub,"sigma_pub/D");
  nhitTree->Branch("nhit_avg",&nhitavg,"nhit_avg/D");
  nhitTree->Branch("bg_avg",&nhitbgavg,"bg_avg/D");
  nhitTree->Branch("sin_avg",&nhitsinavg,"sin_avg/D");

  double nhits[m_nhitbins], nhitserr[m_nhitbins], nhitres[m_nhitbins], nhitreserr[m_nhitbins];

  // Fit the histograms
  double avg_sigma = 0.0;
  for( int i = 0; i < m_nhitbins; ++i ){

    // fill some details for this bin
    nhitavg = sumnhit[i]/sumsize[i];
    nhitbgavg = sumbg[i]/sumsize[i];
    nhitsinavg = sumsin[i]/sumsize[i];

    // fit the dE/dx distribution in bins of nhit
    int fs = hdedx_nhit[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
    double mean = hdedx_nhit[i]->GetFunction("gaus")->GetParameter(1);
    double width = hdedx_nhit[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_nhit[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;

    nhitdedx = hdedx_nhit[i]->GetFunction("gaus")->GetParameter(1);
    nhitsigma = hdedx_nhit[i]->GetFunction("gaus")->GetParameter(2);


    // fit the dE/dx distribution in bins of nhit
    fs = hdedx_nhit_pub[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
    mean = hdedx_nhit_pub[i]->GetFunction("gaus")->GetParameter(1);
    width = hdedx_nhit_pub[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_nhit_pub[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;

    nhitdedxpub = hdedx_nhit_pub[i]->GetFunction("gaus")->GetParameter(1);
    nhitsigmapub = hdedx_nhit_pub[i]->GetFunction("gaus")->GetParameter(2);


    // same some information for fitting sigma vs nhit
    nhits[i] = nhitavg;
    nhitserr[i] = 0.0;
    nhitres[i] = nhitsigma;
    nhitreserr[i] = hdedx_nhit[i]->GetFunction("gaus")->GetParError(2);

    avg_sigma += nhitres[i];

    // fill the tree for this bin
    nhitTree->Fill();
  }

  avg_sigma = avg_sigma/m_nhitbins;
  for( int i = 0; i < m_nhitbins; ++i ){
    nhitres[i] = nhitres[i]/avg_sigma;
    nhitreserr[i] = nhitreserr[i]/avg_sigma;
  }

  // Print the histograms for quality control
  TCanvas* ctmp = new TCanvas("tmp","tmp",900,900);
  ctmp->Divide(3,3);
  std::stringstream psname; psname << "plots/sigmavsnhit_fits.ps[";
  ctmp->Print(psname.str().c_str());
  psname.str(""); psname << "plots/sigmavsnhit_fits.ps";
  for( int i = 0 ; i < m_nhitbins; ++i ){
    ctmp->cd(i%9+1);
    hdedx_nhit[i]->Draw();
    if((i+1)%9==0)
      ctmp->Print(psname.str().c_str());
  }
  psname.str(""); psname << "plots/sigmavsnhit_fits.ps]";
  ctmp->Print(psname.str().c_str());
  delete ctmp;


  // --------------------------------------------------
  // FIT SIGMA VS NHIT CURVE
  // --------------------------------------------------


  TCanvas* sigvsnhitcan = new TCanvas("sigvsnhitcan","#sigma vs. nHit",600,600);

  TGraphErrors* gr = new TGraphErrors(m_nhitbins,nhits,nhitres,nhitserr,nhitreserr);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.3);
  gr->SetTitle(";nHit;#sigma");
  /*
  // Set stat options
  gStyle->SetOptStat(1111111);
  // Set y-position (fraction of pad size)
  gStyle->SetStatY(0.65);                
  // Set x-position (fraction of pad size)
  gStyle->SetStatX(0.9);                
  // Set width of stat-box (fraction of pad size)
  gStyle->SetStatW(0.15);                
  // Set height of stat-box (fraction of pad size)
  gStyle->SetStatH(0.15);                
  */
  gr->Draw("AP");

  WidgetSigma* gs = new WidgetSigma();

  double nhitmin = 2, nhitmax = 50;
  TF1* fsigma = new TF1("fsigma",gs,nhitmin,nhitmax,6,"WidgetSigma");
  fsigma->SetParameter(0,2);
  for( int i = 1; i < 6; ++i ){
    fsigma->SetParameter(i,gpar.getNHitPars(i-1));
  }

  // if the fit succeeds, write out the new parameters
  int status = gr->Fit("fsigma","qm","",nhitmin,nhitmax);
  for( int i = 0; i < 10; ++i ){
    if( status == 0 ) break;
    status = gr->Fit("fsigma","q","",nhitmin,nhitmax);
  }

  if( status != 0 ) std::cout << "WARNING: SIGMA VS NHIT FIT FAILED..." << std::endl;
  for( int j = 1; j < 6; ++j ){
    gpar.setNHitPars(j-1,fsigma->GetParameter(j));
  }
  
  fsigma->Draw("same");
  sigvsnhitcan->SaveAs("plots/sigmavsnhit.eps");
  delete sigvsnhitcan;

  // write out the (possibly) updated parameters to file
  gpar.printParameters("parameters.sigmanhit.fit");
}

void
HadronCalibration::fitSigmaVsSin( TString filename, std::string paramfile, int mcflag = 0, int type = 0 ){

  TFile* infile = new TFile(filename);
  TTree* hadron = (TTree*)infile->Get("track");
  
  std::cout << "HadronCalibration: reading in file " << filename << std::endl;

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  double db;      // the nearest distance to the IP in the xy plane
  double dz;      // the nearest distance to the IP in the z plane
  double chiPi;   // PID chi value
  double eop;     // energy over momentum in the calorimeter

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  hadron->SetBranchAddress("dedx", &dedxpub);
  hadron->SetBranchAddress("dedxsat", &dedx);
  hadron->SetBranchAddress("pF", &p);
  hadron->SetBranchAddress("costh", &costh);
  hadron->SetBranchAddress("db", &db);
  hadron->SetBranchAddress("dz", &dz);
  hadron->SetBranchAddress("chiPi", &chiPi);
  hadron->SetBranchAddress("eopst", &eop);

  if( type == 0 )
    hadron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( type == 1 )
    hadron->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    hadron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  const int m_sinbins = 20;
  const double m_uppersin = 1.0, m_lowersin = 0.36;

  double sinstep = (m_uppersin-m_lowersin)/m_sinbins;

  // Create the histograms to be fit
  TH1F* hdedx_sin[m_sinbins];
  TH1F* hdedx_sin_pub[m_sinbins];

  // initialize the histograms
  for( int i = 0; i < m_sinbins; ++i ){
    char histname[100],histname2[100],histname3[100],histname4[100],histname5[100],histname6[100];
    sprintf(histname, "dedx_sin_%d",i);
    sprintf(histname4, "dedxpub_sin_%d",i);

    hdedx_sin[i] = new TH1F(histname,histname,200,-1,1);
    hdedx_sin_pub[i] = new TH1F(histname4,histname4,200,-1,1);
  }

  // Create some containers to calculate averages
  double sumnhit[m_sinbins];
  double sumbg[m_sinbins];
  double sumsin[m_sinbins];
  int sumsize[m_sinbins];
  for( int i = 0; i < m_sinbins; ++i ){
    sumnhit[i] = 0;
    sumbg[i] = 0;
    sumsin[i] = 0;
    sumsize[i] = 0;
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
  WidgetParameterization gpar(paramfile);

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill the histograms to be fitted
  for( unsigned int index = 0; index < hadron->GetEntries(); ++index ){
    hadron->GetEvent(index);

    if( costh != costh ) continue;

    bg = fabs(p)/Widget::mpion;
    double sin = sqrt(1-costh*costh);

    int nhit;
    if( type == 0 )
      nhit = std::floor(b3nhit);
    else if( type == 1 )
      nhit = b2nhit;

    if( dedx <= 0 || nhit <= 0 || nhit >= 100 || sin > m_uppersin || sin < m_lowersin )
      continue;

    int sinBin = (int)((sin-m_lowersin)/(m_uppersin-m_lowersin) * m_sinbins);

    double dedx_new = dedx;
    if( !mcflag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_cur = gpar.dedxPrediction(bg);
    double dedx_res = gpar.sigmaPrediction(dedx,nhit,sqrt(1-costh*costh));
    double chi_new  = (dedx_new-dedx_cur)/dedx_res;

    double res_cor = gpar.nhitPrediction(nhit);
    int i = (int) ( (sin-m_lowersin)/sinstep );
    hdedx_sin[i]->Fill((dedx_new-dedx_cur)/res_cor);
    hdedx_sin_pub[i]->Fill((dedx_new-dedx_cur)/res_cor);

    if( fabs(sin-(m_lowersin+0.5*sinstep+sinBin*sinstep)) < 0.5*sinstep ){
      sumnhit[sinBin] += nhit;
      sumbg[sinBin] += bg;
      sumsin[sinBin] += sqrt(1-costh*costh);
      sumsize[sinBin] += 1;
    }
  }// end of event loop


  // --------------------------------------------------
  // FIT IN BINS OF SIN
  // --------------------------------------------------

  // fit the histograms with Gaussian functions
  // and extract the means and errors
  TTree* sinTree = new TTree("sigmavssin","dE/dx means and errors");

  double sindedx;          // dE/dx mean value for this bin
  double sinsigma;         // dE/dx width value for this bin
  double sindedxpub;       // dE/dx mean value for this bin
  double sinsigmapub;      // dE/dx width value for this bin
  double sinnhitavg;       // average nhit value for this bin
  double sinbgavg;         // average bg value for this bin
  double sinavg;           // average sin(theta) value for this bin

  sinTree->Branch("dedx",&sindedx,"dedx/D");
  sinTree->Branch("width",&sinsigma,"width/D");
  sinTree->Branch("dedx_pub",&sindedxpub,"dedx_pub/D");
  sinTree->Branch("sigma_pub",&sinsigmapub,"sigma_pub/D");
  sinTree->Branch("nhit_avg",&sinnhitavg,"nhit_avg/D");
  sinTree->Branch("bg_avg",&sinbgavg,"bg_avg/D");
  sinTree->Branch("sin_avg",&sinavg,"sin_avg/D");

  double sins[m_sinbins], sinserr[m_sinbins], sinres[m_sinbins], sinreserr[m_sinbins];

  // Fit the histograms
  double avg_sigma = 0.0;
  for( int i = 0; i < m_sinbins; ++i ){

    // fill some details for this bin
    sinnhitavg = sumnhit[i]/sumsize[i];
    sinbgavg = sumbg[i]/sumsize[i];
    sinavg = sumsin[i]/sumsize[i];


    // fit the dE/dx distribution in bins of sin
    int fs = hdedx_sin[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
    double mean = hdedx_sin[i]->GetFunction("gaus")->GetParameter(1);
    double width = hdedx_sin[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_sin[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;

    sindedx = hdedx_sin[i]->GetFunction("gaus")->GetParameter(1);
    sinsigma = hdedx_sin[i]->GetFunction("gaus")->GetParameter(2);


    // fit the dE/dx distribution in bins of sin
    fs = hdedx_sin_pub[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
    mean = hdedx_sin_pub[i]->GetFunction("gaus")->GetParameter(1);
    width = hdedx_sin_pub[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_sin_pub[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;

    sindedxpub = hdedx_sin_pub[i]->GetFunction("gaus")->GetParameter(1);
    sinsigmapub = hdedx_sin_pub[i]->GetFunction("gaus")->GetParameter(2);


    // same some information for fitting sigma vs sin
    sins[i] = sinavg;
    sinserr[i] = 0.0;
    sinres[i] = sinsigma;
    sinreserr[i] = hdedx_sin[i]->GetFunction("gaus")->GetParError(2);

    avg_sigma += sinres[i];

    // fill the tree for this bin
    sinTree->Fill();
  }

  avg_sigma = avg_sigma/m_sinbins;
  for( int i = 0; i < m_sinbins; ++i ){
    sinres[i] = sinres[i]/avg_sigma;
    sinreserr[i] = sinreserr[i]/avg_sigma;
  }

  // Print the histograms for quality control
  TCanvas* ctmp = new TCanvas("tmp","tmp",900,900);
  ctmp->Divide(3,3);
  std::stringstream psname; psname << "plots/sigmavssin_fits.ps[";
  ctmp->Print(psname.str().c_str());
  psname.str(""); psname << "plots/sigmavssin_fits.ps";
  for( int i = 0 ; i < m_sinbins; ++i ){
    ctmp->cd(i%9+1);
    hdedx_sin[i]->Draw();
    if((i+1)%9==0)
      ctmp->Print(psname.str().c_str());
  }
  psname.str(""); psname << "plots/sigmavssin_fits.ps]";
  ctmp->Print(psname.str().c_str());
  delete ctmp;


  // --------------------------------------------------
  // FIT SIGMA VS SIN CURVE
  // --------------------------------------------------


  TCanvas* sigvssincan = new TCanvas("sigvssincan","#sigma vs. sin(#theta)",600,600);

  TGraphErrors* gr = new TGraphErrors(m_sinbins,sins,sinres,sinserr,sinreserr);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.3);
  gr->SetTitle(";sin(#theta);#sigma");
  gr->Draw("AP");

  WidgetSigma* gs = new WidgetSigma();

  TF1* fsigma = new TF1("fsigma",gs,m_lowersin,m_uppersin,6,"WidgetSigma");
  fsigma->SetParameter(0,2);
  for( int i = 1; i < 6; ++i ){
    fsigma->SetParameter(i,gpar.getSinPars(i-1));
  }

  // if the fit succeeds, write out the new parameters
  int status = gr->Fit("fsigma","q","",m_lowersin,m_uppersin);
  for( int i = 0; i < 10; ++i ){
    if( status == 0 ) break;
    status = gr->Fit("fsigma","q","",m_lowersin,m_uppersin);
  }

  if( status != 0 ) std::cout << "WARNING: SIGMA VS SIN FIT FAILED..." << std::endl;
  for( int j = 1; j < 6; ++j ){
    gpar.setSinPars(j-1,fsigma->GetParameter(j));
  }

  fsigma->Draw("same");
  sigvssincan->SaveAs("plots/sigmavssin.eps");
  delete sigvssincan;

  // write out the (possibly) updated parameters to file
  gpar.printParameters("parameters.sigmasin.fit");
}




void
HadronCalibration::plotEfficiency( TString filenames[5], TString saveas, std::string paramfile, int mcflag = 0, int type = 0 ){

  TFile* pifile = new TFile(filenames[0]);
  TFile* kfile = new TFile(filenames[1]);
  TTree* pion = (TTree*)pifile->Get("track");
  TTree* kaon = (TTree*)kfile->Get("track");
  
  std::cout << "HadronCalibration: reading in file " << filenames[0] << std::endl;

  // --------------------------------------------------
  // INITIALIZE CONTAINERS
  // --------------------------------------------------

  double dedxpub; // dE/dx without electron saturation correction
  double dedx;    // dE/dx with electron saturation correction
  double p;       // track momentum
  double bg;      // track beta-gamma
  double costh;   // cosine of track polar angle
  double db;      // the nearest distance to the IP in the xy plane
  double dz;      // the nearest distance to the IP in the z plane
  double chiPi;   // PID chi value
  double eop;     // energy over momentum in the calorimeter

  // Belle II variables
  int b2nhit;       // number of hits on this track

  // BES III variables
  double b3nhit;    // number of hits on this track

  pion->SetBranchAddress("dedx", &dedxpub);
  pion->SetBranchAddress("dedxsat", &dedx);
  pion->SetBranchAddress("pF", &p);
  pion->SetBranchAddress("costh", &costh);
  pion->SetBranchAddress("db", &db);
  pion->SetBranchAddress("dz", &dz);
  pion->SetBranchAddress("chiPi", &chiPi);
  pion->SetBranchAddress("eopst", &eop);

  if( type == 0 )
    pion->SetBranchAddress("numGoodHits", &b3nhit);
  else if( type == 1 )
    pion->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    pion->SetBranchAddress("numGoodLayerHits", &b2nhit);


  kaon->SetBranchAddress("dedx", &dedxpub);
  kaon->SetBranchAddress("dedxsat", &dedx);
  kaon->SetBranchAddress("pF", &p);
  kaon->SetBranchAddress("costh", &costh);
  kaon->SetBranchAddress("db", &db);
  kaon->SetBranchAddress("dz", &dz);
  kaon->SetBranchAddress("chiPi", &chiPi);
  kaon->SetBranchAddress("eopst", &eop);

  if( type == 0 )
    kaon->SetBranchAddress("numGoodHits", &b3nhit);
  else if( type == 1 )
    kaon->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    kaon->SetBranchAddress("numGoodLayerHits", &b2nhit);

  const int m_nhitbins = 20;
  const double m_uppernhit = 32.5, m_lowernhit = 5.5;

  double nhitstep = (m_uppernhit-m_lowernhit)/m_nhitbins;


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
  else
    std::cout << "WARNING: NO SATUARTION PARAMETERS!!!" << std::endl;

  WidgetParameterization gpar(paramfile);

  const int nbins = 20;
  TFile* outfile = new TFile(saveas,"RECREATE");

  TH1F* pinumer = new TH1F("pinumer","",nbins,0,2);
  TH1F* pidenom = new TH1F("pidenom","",nbins,0,2);

  TH2F* pihdedx = new TH2F("pihdedx","",100,0,2,100,0,4);
  TH2F* pipredpi = new TH2F("pipredpi","",100,0,2,100,0,4);
  TH2F* pipredk = new TH2F("pipredk","",100,0,2,100,0,4);

  TH2F* pichivp = new TH2F("pichivp","",100,0,2,100,-10,10);

  TH1F* knumer = new TH1F("knumer","",nbins,0,2);
  TH1F* kdenom = new TH1F("kdenom","",nbins,0,2);

  TH2F* khdedx = new TH2F("khdedx","",100,0,2,100,0,4);
  TH2F* kpredpi = new TH2F("kpredpi","",100,0,2,100,0,4);
  TH2F* kpredk = new TH2F("kpredk","",100,0,2,100,0,4);

  TH2F* kchivp = new TH2F("kchivp","",100,0,2,100,-10,10);

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill the pions
  for( unsigned int index = 0; index < pion->GetEntries(); ++index ){
    pion->GetEvent(index);

    int nhit;
    if( type == 0 )
      nhit = std::floor(b3nhit);
    else if( type == 1 )
      nhit = b2nhit;

    if( dedx <= 0 || nhit <= m_lowernhit || nhit >= m_uppernhit )
      continue;

    double dedx_new = dedx;
    if( !mcflag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_res = gpar.sigmaPrediction(dedx,nhit,sqrt(1-costh*costh));

    pihdedx->Fill(fabs(p),dedx_new);

    double pibg = fabs(p)/Widget::mpion;
    double dedx_cur_pi = gpar.dedxPrediction(pibg);
    double chi_new_pi  = (dedx_new-dedx_cur_pi)/dedx_res;

    pipredpi->Fill(fabs(p),dedx_cur_pi);

    double kbg = fabs(p)/Widget::mkaon;
    double dedx_cur_k = gpar.dedxPrediction(kbg);
    double chi_new_k  = (dedx_new-dedx_cur_k)/dedx_res;

    pipredk->Fill(fabs(p),dedx_cur_k);

    pichivp->Fill(fabs(p),(chi_new_k*chi_new_k-chi_new_pi*chi_new_pi));
    
    pidenom->Fill(fabs(p));
    if( chi_new_k*chi_new_k > chi_new_pi*chi_new_pi ) pinumer->Fill(fabs(p));
  }


  // Fill the kaons
  for( unsigned int index = 0; index < kaon->GetEntries(); ++index ){
    kaon->GetEvent(index);

    int nhit;
    if( type == 0 )
      nhit = std::floor(b3nhit);
    else if( type == 1 )
      nhit = b2nhit;

    if( dedx <= 0 || nhit <= m_lowernhit || nhit >= m_uppernhit )
      continue;
    if( fabs(p) < 0.5 and dedx < 1 ) continue;

    double dedx_new = dedx;
    if( !mcflag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_res = gpar.sigmaPrediction(dedx,nhit,sqrt(1-costh*costh));

    khdedx->Fill(fabs(p),dedx_new);

    double pibg = fabs(p)/Widget::mpion;
    double dedx_cur_pi = gpar.dedxPrediction(pibg);
    double chi_new_pi  = (dedx_new-dedx_cur_pi)/dedx_res;

    kpredpi->Fill(fabs(p),dedx_cur_pi);

    double kbg = fabs(p)/Widget::mkaon;
    double dedx_cur_k = gpar.dedxPrediction(kbg);
    double chi_new_k  = (dedx_new-dedx_cur_k)/dedx_res;

    kpredk->Fill(fabs(p),dedx_cur_k);

    kchivp->Fill(fabs(p),(chi_new_k*chi_new_k-chi_new_pi*chi_new_pi));
    
    kdenom->Fill(fabs(p));
    if( chi_new_k*chi_new_k < chi_new_pi*chi_new_pi ) knumer->Fill(fabs(p));
  }

  TEfficiency* piteff = new TEfficiency(*pinumer,*pidenom);
  TEfficiency* kteff = new TEfficiency(*knumer,*kdenom);

  outfile->cd();
  pinumer->Write();
  pidenom->Write();
  piteff->Write();
  pihdedx->Write();
  pipredpi->Write();
  pipredk->Write();
  pichivp->Write();

  knumer->Write();
  kdenom->Write();
  kteff->Write();
  khdedx->Write();
  kpredpi->Write();
  kpredk->Write();
  kchivp->Write();
  outfile->Close();
}
