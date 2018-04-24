#include "HadronPrep.h"

HadronPrep::HadronPrep() :
  m_filename("widget.gen.root"),
  m_mcFlag(1),
  m_type(0),
  m_bgbins(15),
  m_upperbg(5.0),
  m_lowerbg(0.1),
  m_cosbins(18),
  m_uppercos(0.93),
  m_lowercos(0.0)
{}

HadronPrep::HadronPrep( TString infile, int mcFlag, int type, int bgbins, double upperbg, double lowerbg, int cosbins, double uppercos, double lowercos ){
  m_filename = infile;
  m_mcFlag = mcFlag;
  m_type = type;
  m_bgbins = bgbins;
  m_upperbg = upperbg;
  m_lowerbg = lowerbg;
  m_cosbins = cosbins;
  m_uppercos = uppercos;
  m_lowercos = lowercos;
}

void
HadronPrep::FormatGraph( TGraphErrors* gr, int flag ){
  if( flag == 0 ){
    gr->SetTitle(";cos(#theta);#chi_{mean}");
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(.7);
    gr->SetMaximum(1);
    gr->SetMinimum(-1);
    gr->GetXaxis()->SetRangeUser(-1,1);
  }
  else if( flag == 1 ){
    gr->SetMarkerStyle(22);
    gr->SetMarkerColor(9);
    gr->SetMarkerSize(.7);
  }
  else if( flag == 2 ){
    gr->SetTitle(";cos(#theta);#sigma");
    gr->SetMarkerStyle(24);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(.7);
    gr->SetMaximum(2);
    gr->SetMinimum(0);
    gr->GetXaxis()->SetRangeUser(-1,1);
  }
}

void
HadronPrep::bgCosThetaFits( TString particle, TFile* outfile, bool correct, std::string paramfile ){

  // supress TCanvas::Print messages
  gErrorIgnoreLevel = kWarning;

  m_gpar.setParameters(paramfile);

  double mass = 0.0;
  if( particle == "pion" ) mass = Widget::mpion;
  else if( particle == "kaon" ) mass = Widget::mkaon;
  else if( particle == "proton" ) mass = Widget::mproton;
  else if( particle == "muon" ) mass = Widget::mmuon;
  else if( particle == "electron" ) mass = Widget::melectron;
  if( mass == 0.0 ) exit(1);

  TFile* infile = new TFile(m_filename);
  TTree* hadron = (TTree*)infile->Get("track");
  
  std::cout << "HadronPrep: reading in file " << m_filename << std::endl;
  std::cout << "\t beta-gamma range: " << m_lowerbg << " to " << m_upperbg << std::endl;
  std::cout << "\t cos(theta) range: " << m_lowercos << " to " << m_uppercos << std::endl;

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

  if( m_type == 0 )
    hadron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( m_type == 1 )
    hadron->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    hadron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  double bgstep = (m_upperbg-m_lowerbg)/m_bgbins;
  double cosstep = (m_uppercos-m_lowercos)/m_cosbins;

  // Create the histograms to be fit (before correction)
  TH1F* hdedx_costh[m_bgbins][m_cosbins];
  TH1F* hdedx_costh_pub[m_bgbins][m_cosbins];
  TH1F* hchi_costh[m_bgbins][m_cosbins];
  TH1F* hchi_costh_pub[m_bgbins][m_cosbins];

  // initialize the histograms
  for( int i = 0; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      char histname[100],histname2[100],histname3[100],histname4[100],histname5[100];
      sprintf(histname, "dedx_costh_p_%d_cos_%d",i,j);
      sprintf(histname2, "chi_costh_p_%d_cos_%d",i,j);
      sprintf(histname3, "chipub_costh_p_%d_cos_%d",i,j);
      sprintf(histname4, "dedxpub_costh_p_%d_cos_%d",i,j);
      sprintf(histname5, "dedxnew_costh_p_%d_cos_%d",i,j);


      if( particle == "electron" ){
	hdedx_costh[i][j] = new TH1F(histname,histname,200,0,2);
	hchi_costh[i][j] = new TH1F(histname2,histname2,200,-2,2);
	hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,-5,5);
	hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,2);
      }
      else if( particle == "proton" ){
	if( i < 5 ){
	  hdedx_costh[i][j] = new TH1F(histname,histname,200,0,6);
	  hchi_costh[i][j] = new TH1F(histname2,histname2,200,-30,30);
	  hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,0,40);
	  hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,6);
	}
	else if( i < 9 ){
	  hdedx_costh[i][j] = new TH1F(histname,histname,200,0,4);
	  hchi_costh[i][j] = new TH1F(histname2,histname2,200,-5,5);
	  hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,-20,20);
	  hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,4);
	}
	else{
	  hdedx_costh[i][j] = new TH1F(histname,histname,200,0,2);
	  hchi_costh[i][j] = new TH1F(histname2,histname2,200,-5,5);
	  hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,-5,5);
	  hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,2);
	}
      }
      else if( i < 2 ){
	hdedx_costh[i][j] = new TH1F(histname,histname,200,0,4);
	hchi_costh[i][j] = new TH1F(histname2,histname2,200,-5,5);
	hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,-5,5);
	hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,4);
      }
      else{
	hdedx_costh[i][j] = new TH1F(histname,histname,200,0,2);
	hchi_costh[i][j] = new TH1F(histname2,histname2,200,-5,5);
	hchi_costh_pub[i][j] = new TH1F(histname3,histname3,200,-5,5);
	hdedx_costh_pub[i][j] = new TH1F(histname4,histname4,200,0,2);
      }
    }
  }

  // Create some containers to calculate averages
  double sumcos[m_bgbins][m_cosbins];
  double sumsin[m_bgbins][m_cosbins];
  double sumbg[m_bgbins][m_cosbins];
  double sumressq[m_bgbins][m_cosbins];
  int sumsize[m_bgbins][m_cosbins];
  for( int i = 0; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      sumcos[i][j] = 0;
      sumsin[i][j] = 0;
      sumbg[i][j] = 0;
      sumressq[i][j] = 0;
      sumsize[i][j] = 0;
    }
  }

  // Create some containers for means and errors
  double means[m_bgbins][m_cosbins];
  double errors[m_bgbins][m_cosbins];
  double widths[m_bgbins][m_cosbins];
  for( int i = 0; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      means[i][j] = 0;
      errors[i][j] = 0;
      widths[i][j] = 0;
    }
  }


  // get the hadron saturation parameters
  // if the parameters do not exist, use the values in the default constructor
  m_hadsat.setParameters("sat-pars.fit.txt");

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill the histograms to be fitted
  for( unsigned int index = 0; index < hadron->GetEntries(); ++index ){
    hadron->GetEvent(index);

    bg = fabs(p)/mass;
    costh = fabs(costh);

    int nhit;
    if( m_type == 0 )
      nhit = std::floor(b3nhit);
    else if( m_type == 1 )
      nhit = b2nhit;

    // clean up bad events and restrict the momentum range
    if( nhit < 0 || nhit > 100 || dedx <= 0 || costh != costh ||
	bg <= m_lowerbg || bg >= m_upperbg || costh <= m_lowercos || costh >= m_uppercos )
      continue;

    // apply an E/p cut for pions
    if( particle == "pion" && eop > 0.75 ) continue;

    // use loose dE/dx restrictions to remove contaminations
    if( particle == "proton" && (dedxpub < 0 || dedxpub > 100 || dedxpub < 0.85/fabs(p)) )
      continue;
    if( particle == "kaon" && dedxpub < 0.4/fabs(p) )
      continue;

    int bgBin = (int)((bg-m_lowerbg)/(m_upperbg-m_lowerbg) * m_bgbins);
    int cosBin = (int)((costh-m_lowercos)/(m_uppercos-m_lowercos) * m_cosbins);

    double dedx_pre = dedx;
    double dedx_new = dedx_pre;
    if( !m_mcFlag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_cur = m_gpar.dedxPrediction(bg);
    double dedx_res = m_gpar.sigmaPrediction(dedx_cur,nhit,sqrt(1-costh*costh));
    if( dedx_res == 0 ){
      std::cout << "RESOLUTION IS ZERO!!!" << std::endl;
      continue;
    }

    double chi_new  = (dedx_new-dedx_cur)/dedx_res;
    if( particle == "electron" ) chi_new = (dedx_new-1)/dedx_res;
    if( correct ) hdedx_costh[bgBin][cosBin]->Fill(dedx_new);
    else hdedx_costh[bgBin][cosBin]->Fill(dedx_pre);

    if( correct ) hchi_costh[bgBin][cosBin]->Fill(chi_new);
    else hchi_costh[bgBin][cosBin]->Fill(chiPi);

    hdedx_costh_pub[bgBin][cosBin]->Fill(dedxpub);
    hchi_costh_pub[bgBin][cosBin]->Fill(chiPi);

    //    if( fabs(bg-(m_lowerbg+0.5*bgstep+bgBin*bgstep)) < 0.5*bgstep && 
    //	fabs(costh-(m_lowercos+0.5*cosstep+cosBin*cosstep)) < 0.5*cosstep ){
      sumcos[bgBin][cosBin] += costh;
      sumsin[bgBin][cosBin] += sqrt(1-costh*costh);
      sumbg[bgBin][cosBin] += bg;
      sumressq[bgBin][cosBin] += pow(dedx_res,2);
      sumsize[bgBin][cosBin] += 1;
      //    }

  }// end of event loop

  // --------------------------------------------------
  // FIT IN BINS OF BETA-GAMMA AND COS(THETA)
  // --------------------------------------------------

  // fit the histograms with Gaussian functions
  // and extract the means and errors
  TTree* satTree = new TTree(particle,"dE/dx means and errors");

  double satbg;          // beta-gamma value for this bin
  double satcosth;       // cos(theta) value for this bin
  double satdedx;        // mean dE/dx value for this bin
  double saterror;       // error on ^
  double satbg_avg;      // average beta-gamma value for this sample
  double satcosth_avg;   // average cos(theta) value for this sample
  double satsinth_avg;   // average sin(theta) value for this sample
  double satdedxres_avg; // average dE/dx error squared for this sample
  double satpubdedx;     // mean "public" dE/dx value for this bin
  double satpuberror;    // error on ^
  double satpubwidth;    // width of ^ distribution
  double satchi;         // mean chi value for this bin
  double satchierr;      // error on ^
  double satchiwidth;    // width of ^ distribution
  double satchipub;      // mean "public" chi value for this bin
  double satchipuberr;   // error on ^
  double satchipubwidth; // width of ^ distribution
  double ratio;          // ratio of the predicted mean to that of the average

  satTree->Branch("bg",&satbg,"bg/D");
  satTree->Branch("costh",&satcosth,"costh/D");
  satTree->Branch("dedx",&satdedx,"dedx/D");
  satTree->Branch("error",&saterror,"error/D");
  satTree->Branch("bg_avg",&satbg_avg,"bg_avg/D");
  satTree->Branch("costh_avg",&satcosth_avg,"costh_avg/D");
  satTree->Branch("sinth_avg",&satsinth_avg,"sinth_avg/D");
  satTree->Branch("dedxres_avg",&satdedxres_avg,"dedxres_avg/D");
  satTree->Branch("dedx_pub",&satpubdedx,"dedx_pub/D");
  satTree->Branch("dedxerr_pub",&satpuberror,"dedxerr_pub/D");
  satTree->Branch("chi",&satchi,"chi/D");
  satTree->Branch("chi_pub",&satchipub,"chi_pub/D");
  satTree->Branch("sigma",&satchiwidth,"sigma/D");
  satTree->Branch("sigma_pub",&satchipubwidth,"sigma_pub/D");
  satTree->Branch("ratio",&ratio,"ratio/D");

  // Fit the histograms
  int fs1=0, fs2=0, fs3=0, fs4=0;
  double dedxmax = 1.0, dedxmin = 1.0;
  for( int i = 0; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){

      // fill some details for this bin
      satbg = m_lowerbg+0.5*bgstep+i*bgstep;
      satcosth = m_lowercos+0.5*cosstep+j*cosstep;

      satbg_avg = sumbg[i][j]/sumsize[i][j];
      satcosth_avg = sumcos[i][j]/sumsize[i][j];
      satsinth_avg = sumsin[i][j]/sumsize[i][j];
      satdedxres_avg = sumressq[i][j]/sumsize[i][j];

      int fs;
      // fit the dE/dx distribution in bins of beta-gamma and cosine
      fs = hdedx_costh[i][j]->Fit("gaus","ql");
      double mean = hdedx_costh[i][j]->GetFunction("gaus")->GetParameter(1);
      double width = hdedx_costh[i][j]->GetFunction("gaus")->GetParameter(2);
      fs = hdedx_costh[i][j]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
      if( fs != 0 ){
	// std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
	fs1++;
      }

      satdedx = hdedx_costh[i][j]->GetFunction("gaus")->GetParameter(1);
      saterror = hdedx_costh[i][j]->GetFunction("gaus")->GetParError(1);
      means[i][j] = satdedx;
      errors[i][j] = saterror;
      widths[i][j] = hdedx_costh[i][j]->GetFunction("gaus")->GetParameter(2);

      if( satdedx > dedxmax ) dedxmax = satdedx;
      if( satdedx < dedxmin ) dedxmin = satdedx;

      // fit the dE/dx distribution without correction in bins of beta-gamma and cosine
      fs = hdedx_costh_pub[i][j]->Fit("gaus","ql");
      mean = hdedx_costh_pub[i][j]->GetFunction("gaus")->GetParameter(1);
      width = hdedx_costh_pub[i][j]->GetFunction("gaus")->GetParameter(2);
      fs = hdedx_costh_pub[i][j]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
      if( fs != 0 ){
	// std::cout << "\t\t" <<"MEAN PUB FIT STATUS " << fs << std::endl;
	fs2++;
      }

      satpubdedx = hdedx_costh_pub[i][j]->GetFunction("gaus")->GetParameter(1);
      satpuberror = hdedx_costh_pub[i][j]->GetFunction("gaus")->GetParError(1);
      satpubwidth = hdedx_costh_pub[i][j]->GetFunction("gaus")->GetParameter(2);


      // fit the chi distribution
      fs = hchi_costh[i][j]->Fit("gaus","ql");
      mean = hchi_costh[i][j]->GetFunction("gaus")->GetParameter(1);
      width = hchi_costh[i][j]->GetFunction("gaus")->GetParameter(2);
      fs = hchi_costh[i][j]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
      if( fs != 0 ){
	// std::cout << "\t\t" <<"CHI FIT STATUS " << fs << std::endl;
	fs3++;
      }

      satchi = hchi_costh[i][j]->GetFunction("gaus")->GetParameter(1);
      satchierr = hchi_costh[i][j]->GetFunction("gaus")->GetParError(1);
      satchiwidth = hchi_costh[i][j]->GetFunction("gaus")->GetParameter(2);


      // fit the chi distribution without correction
      /*
      fs = hchi_costh_pub[i][j]->Fit("gaus","ql");	    
      mean = hchi_costh_pub[i][j]->GetFunction("gaus")->GetParameter(1);
      width = hchi_costh_pub[i][j]->GetFunction("gaus")->GetParameter(2);
      fs = hchi_costh_pub[i][j]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
      if( fs != 0 ){
	// std::cout << "\t\t" <<"CHI PUB FIT STATUS " << fs << std::endl;
	fs4++;
      }

      satchipub = hchi_costh_pub[i][j]->GetFunction("gaus")->GetParameter(1);
      satchipuberr = hchi_costh_pub[i][j]->GetFunction("gaus")->GetParError(1);
      satchipubwidth = hchi_costh_pub[i][j]->GetFunction("gaus")->GetParameter(2);
      */

      // determine the ratio of the predicted mean at a given bg to that of the average
      ratio = m_gpar.dedxPrediction(satbg_avg)/m_gpar.dedxPrediction(satbg);

      // fill the tree for this bin
      satTree->Fill();
    }
  }

  if( fs1+fs2+fs3+fs4 != 0 ) std::cout << "FIT RESULTS: " << particle << std::endl;
  if( fs1 != 0 ) std::cout << "\t\t" <<"MEAN FIT FAILS " << fs1 << std::endl;
  if( fs2 != 0 ) std::cout << "\t\t" <<"MEAN PUB FIT FAILS " << fs2 << std::endl;
  if( fs3 != 0 ) std::cout << "\t\t" <<"CHI FIT FAILS " << fs3 << std::endl;
  if( fs4 != 0 ) std::cout << "\t\t" <<"CHI PUB FIT FAILS " << fs4 << std::endl;

  std::string corname("uncorrected");
  if( correct == true ) corname = "corrected";

  // Print the histograms for quality control
  TCanvas* ctmp = new TCanvas("tmp","tmp",900,900);
  ctmp->Divide(3,3);
  std::stringstream psname; psname << "plots/dedx_" << particle << "_" << corname << ".ps[";
  ctmp->Print(psname.str().c_str());
  psname.str(""); psname << "plots/dedx_" << particle << "_" << corname << ".ps";
  for( int i = 0 ; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      ctmp->cd(j%9+1);
      hdedx_costh[i][j]->Draw();
      if((j+1)%9==0)
        ctmp->Print(psname.str().c_str());
    }
  }
  psname.str(""); psname << "plots/dedx_" << particle << "_" << corname << ".ps]";
  ctmp->Print(psname.str().c_str());
  delete ctmp;


  // Print the histograms for quality control
  TCanvas* ctmp2 = new TCanvas("tmp2","tmp2",900,900);
  ctmp2->Divide(3,3);
  psname.str(""); psname << "plots/dedx_" << particle << "_pub_" << corname << ".ps[";
  ctmp2->Print(psname.str().c_str());
  psname.str(""); psname << "plots/dedx_" << particle << "_pub_" << corname << ".ps";
  for( int i = 0 ; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      ctmp2->cd(j%9+1);
      hdedx_costh_pub[i][j]->Draw();
      if((j+1)%9==0)
        ctmp2->Print(psname.str().c_str());
    }
  }
  psname.str(""); psname << "plots/dedx_" << particle << "_pub_" << corname << ".ps]";
  ctmp2->Print(psname.str().c_str());
  delete ctmp2;


  // Print the histograms for quality control
  TCanvas* ctmp3 = new TCanvas("tmp3","tmp3",900,900);
  ctmp3->Divide(3,3);
  psname.str(""); psname << "plots/chi_" << particle << "_" << corname << ".ps[";
  ctmp3->Print(psname.str().c_str());
  psname.str(""); psname << "plots/chi_" << particle << "_" << corname << ".ps";
  for( int i = 0 ; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      ctmp3->cd(j%9+1);
      hchi_costh[i][j]->Draw();
      if((j+1)%9==0)
        ctmp3->Print(psname.str().c_str());
    }
  }
  psname.str(""); psname << "plots/chi_" << particle << "_" << corname << ".ps]";
  ctmp3->Print(psname.str().c_str());
  delete ctmp3;


  // Print the histograms for quality control
  TCanvas* ctmp4 = new TCanvas("tmp4","tmp4",900,900);
  ctmp4->Divide(3,3);
  psname.str(""); psname << "plots/chi_" << particle << "_pub_" << corname << ".ps[";
  ctmp4->Print(psname.str().c_str());
  psname.str(""); psname << "plots/chi_" << particle << "_pub_" << corname << ".ps";
  for( int i = 0 ; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      ctmp4->cd(j%9+1);
      hchi_costh_pub[i][j]->Draw();
      if((j+1)%9==0)
        ctmp4->Print(psname.str().c_str());
    }
  }
  psname.str(""); psname << "plots/chi_" << particle << "_pub_" << corname << ".ps]";
  ctmp4->Print(psname.str().c_str());
  delete ctmp4;

  double cbcenters[m_cosbins], cberrors[m_cosbins];
  for( int j = 0; j < m_cosbins; ++j ){
    cbcenters[j] = m_lowercos+0.5*cosstep+j*cosstep;
    cberrors[j] = 0;
  }

  // Plot the dE/dx means vs. cos(theta) for validation
  TGraphErrors gdedx_costh[m_bgbins];
  for( int i = 0; i < m_bgbins; ++i ){
    char graphname[100];
    satbg = m_lowerbg+0.5*bgstep+i*bgstep;
    sprintf(graphname, "dedx_costh_p_%d",i);
    gdedx_costh[i] = TGraphErrors(m_cosbins,cbcenters,means[i],cberrors,errors[i]);
  }
  TH1F* base = new TH1F("base","dE/dx vs. cos(#theta);cos(#theta);dE/dx",100,0,1.0);
  base->GetXaxis()->SetTitleOffset(1.5);
  base->SetMaximum(dedxmax*1.1);
  base->SetMinimum(dedxmin*0.9);
  base->SetStats(0);

  // Print the histograms for quality control
  TCanvas* ctmp6 = new TCanvas("tmp6","tmp6",900,900);
  base->DrawCopy();
  for( int i = 0 ; i < m_bgbins; ++i ){
    gdedx_costh[i].SetMarkerSize(0.8);
    gdedx_costh[i].SetMarkerColor(50+i);
    gdedx_costh[i].SetLineColor(50+i);
    gdedx_costh[i].DrawClone("same");
  }
  psname.str(""); psname << "plots/costh_" << particle << "_" << corname << ".eps";
  ctmp6->SaveAs(psname.str().c_str());
  delete ctmp6;

  std::cout << "HadronPrep: saving output to " << outfile->GetName() << std::endl;

  // write out the data to file
  outfile->cd();
  satTree->Write();

  delete base;
  for( int i = 0; i < m_bgbins; ++i ){
    for( int j = 0; j < m_cosbins; ++j ){
      delete hdedx_costh[i][j];
      delete hdedx_costh_pub[i][j];
      delete hchi_costh[i][j];
      delete hchi_costh_pub[i][j];
    }
  }
  infile->Close();
}


void
HadronPrep::bgFits( TString particle, TFile* outfile, bool correct, std::string paramfile ){

  // supress TCanvas::Print messages
  gErrorIgnoreLevel = kWarning;

  m_gpar.setParameters(paramfile);

  double mass = 0.0;
  if( particle == "pion" ) mass = Widget::mpion;
  else if( particle == "kaon" ) mass = Widget::mkaon;
  else if( particle == "proton" ) mass = Widget::mproton;
  else if( particle == "muon" ) mass = Widget::mmuon;
  else if( particle == "electron" ) mass = Widget::melectron;
  if( mass == 0.0 ) exit(1);

  TFile* infile = new TFile(m_filename);
  TTree* hadron = (TTree*)infile->Get("track");
  
  std::cout << "HadronPrep: reading in file " << m_filename << std::endl;
  std::cout << "\t beta-gamma range: " << m_lowerbg << " to " << m_upperbg << std::endl;

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

  if( m_type == 0 )
    hadron->SetBranchAddress("numGoodHits", &b3nhit);
  else if( m_type == 1 )
    hadron->SetBranchAddress("lNHitsUsed", &b2nhit);
    //    hadron->SetBranchAddress("numGoodLayerHits", &b2nhit);

  double bgstep = (m_upperbg-m_lowerbg)/m_bgbins;

  // Create the histograms to be fit (before correction)
  TH1F* hdedx_bg[m_bgbins];
  TH1F* hdedx_bg_pub[m_bgbins];
  TH1F* hchi_bg[m_bgbins];
  TH1F* hchi_bg_pub[m_bgbins];
  TH1F* hsigma_bg[m_bgbins];

  // initialize the histograms
  for( int i = 0; i < m_bgbins; ++i ){
    char histname[100], histname2[100], histname3[100];
    char histname4[100], histname5[100], histname6[100];
    sprintf(histname, "dedx_bg_%d",i);
    sprintf(histname2, "chi_bg_%d",i);
    sprintf(histname3, "chipub_bg_%d",i);
    sprintf(histname4, "dedxpub_bg_%d",i);
    sprintf(histname5, "dedxnew_bg_%d",i);
    sprintf(histname6, "sigma_bg_%d",i);

    if( particle == "electron" ){
      hdedx_bg[i] = new TH1F(histname,histname,200,0,2);
      hchi_bg[i] = new TH1F(histname2,histname2,200,-2,2);
      hchi_bg_pub[i] = new TH1F(histname3,histname3,200,-5,5);
      hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,2);
    }
    else if( particle == "proton" ){
      if( i < 5 ){
	hdedx_bg[i] = new TH1F(histname,histname,200,0,6);
	hchi_bg[i] = new TH1F(histname2,histname2,200,-30,30);
	hchi_bg_pub[i] = new TH1F(histname3,histname3,200,0,40);
	hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,6);
      }
      else if( i < 9 ){
	hdedx_bg[i] = new TH1F(histname,histname,200,0,4);
	hchi_bg[i] = new TH1F(histname2,histname2,200,-5,5);
	hchi_bg_pub[i] = new TH1F(histname3,histname3,200,-20,20);
	hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,4);
      }
      else{
	hdedx_bg[i] = new TH1F(histname,histname,200,0,2);
	hchi_bg[i] = new TH1F(histname2,histname2,200,-5,5);
	hchi_bg_pub[i] = new TH1F(histname3,histname3,200,-5,5);
	hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,2);
      }
    }
    else if( i < 2 ){
      hdedx_bg[i] = new TH1F(histname,histname,200,0,4);
      hchi_bg[i] = new TH1F(histname2,histname2,200,-5,5);
      hchi_bg_pub[i] = new TH1F(histname3,histname3,200,-20,20);
      hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,4);
    }
    else{
      hdedx_bg[i] = new TH1F(histname,histname,200,0,2);
      hchi_bg[i] = new TH1F(histname2,histname2,200,-5,5);
      hchi_bg_pub[i] = new TH1F(histname3,histname3,100,-5,5);
      hdedx_bg_pub[i] = new TH1F(histname4,histname4,200,0,2);
    }
    hsigma_bg[i] = new TH1F(histname6,histname6,100,-3,3);
  }

  // create some containers to calculate averages
  double sumcos[m_bgbins];
  double sumsin[m_bgbins];
  double sumbg[m_bgbins];
  double sumressq[m_bgbins];
  int sumsize[m_bgbins];
  for( int i = 0; i < m_bgbins; ++i ){
    sumcos[i] = 0;
    sumsin[i] = 0;
    sumbg[i] = 0;
    sumressq[i] = 0;
    sumsize[i] = 0;
  }

  // create some containers for means and errors
  double means[m_bgbins];
  double errors[m_bgbins];
  double widths[m_bgbins];
  for( int i = 0; i < m_bgbins; ++i ){
    means[i] = 0;
    errors[i] = 0;
    widths[i] = 0;
  }

  // create some containers for checking the residual saturation vs. cos(theta)
  const int cosbins = 20;
  const double cosstep = 2.0 / cosbins;
  double cosArray[cosbins], cosArrayErr[cosbins];
  for( int i = 0; i < cosbins; ++i ){
    cosArray[i] = -1 + (i*cosstep + cosstep/2.0);
    cosArrayErr[i] = 0.0;
  }

  TH1F* hchi_costh_all[2][cosbins];
  TH1F* hchi_costh_0[2][cosbins];
  TH1F* hchi_costh_1[2][cosbins];
  TH1F* hchi_costh_2[2][cosbins];

  for( int i = 0; i < cosbins; ++i ){
    char histname[100], histname0[100], histname1[100], histname2[100];
    sprintf(histname, "chi_costh_pos_%d",i);
    sprintf(histname0, "chi_costh0_pos_%d",i);
    sprintf(histname1, "chi_costh1_pos_%d",i);
    sprintf(histname2, "chi_costh2_pos_%d",i);

    hchi_costh_all[0][i] = new TH1F(histname,histname,100,-5,5);	
    hchi_costh_0[0][i] = new TH1F(histname0,histname0,100,-5,5);
    hchi_costh_1[0][i] = new TH1F(histname1,histname1,100,-5,5);
    hchi_costh_2[0][i] = new TH1F(histname2,histname2,100,-5,5);

    sprintf(histname, "chi_costh_neg_%d",i);
    sprintf(histname0, "chi_costh0_neg_%d",i);
    sprintf(histname1, "chi_costh1_neg_%d",i);
    sprintf(histname2, "chi_costh2_neg_%d",i);
    
    hchi_costh_all[1][i] = new TH1F(histname,histname,100,-5,5);	
    hchi_costh_0[1][i] = new TH1F(histname0,histname0,100,-5,5);
    hchi_costh_1[1][i] = new TH1F(histname1,histname1,100,-5,5);
    hchi_costh_2[1][i] = new TH1F(histname2,histname2,100,-5,5);
  }

  // get the hadron saturation parameters
  // if the parameters do not exist, use the values in the default constructor
  m_hadsat.setParameters("sat-pars.fit.txt");

  // --------------------------------------------------
  // LOOP OVER EVENTS AND FILL CONTAINERS
  // --------------------------------------------------

  // Fill the histograms to be fitted
  for( unsigned int index = 0; index < hadron->GetEntries(); ++index ){
    hadron->GetEvent(index);

    int charge = (p < 0)? 1 : 0;
    bg = fabs(p)/mass;

    int nhit;
    if( m_type == 0 )
      nhit = std::floor(b3nhit);
    else if( m_type == 1 )
      nhit = b2nhit;

    // clean up bad events and restrict the momentum range
    if( nhit < 0 || nhit > 100 || dedx <= 0 || costh != costh ||
	bg <= m_lowerbg || bg >= m_upperbg )
      continue;

    // apply an E/p cut for pions
    if( particle == "pion" && eop > 0.75 ) continue;

    // use loose dE/dx restrictions to remove contaminations
    if( particle == "proton" && (dedxpub < 0 || dedxpub > 100 || dedxpub < 0.85/fabs(p)) )
      continue;
    if( particle == "kaon" && dedxpub < 0.4/fabs(p) )
      continue;

    int bgBin = (int)((bg-m_lowerbg)/(m_upperbg-m_lowerbg) * m_bgbins);

    double dedx_pre = dedx;
    double dedx_new = dedx_pre;
    if( !m_mcFlag ) dedx_new = m_hadsat.D2I(costh,m_hadsat.I2D(costh,1.0)*dedx);
    double dedx_cur = m_gpar.dedxPrediction(bg);
    double dedx_res = m_gpar.sigmaPrediction(dedx_cur,nhit,sqrt(1-costh*costh));

    double chi_new  = (dedx_new-dedx_cur)/dedx_res;
    if( particle == "electron" ) chi_new = (dedx_new-1)/dedx_res;
    double res_cor = m_gpar.resPrediction(nhit,sqrt(1-costh*costh));
    int i = (int) ( (fabs(bg)-m_lowerbg)/bgstep );
    hsigma_bg[i]->Fill((dedx_new-dedx_cur)/res_cor);

    if( correct ) hdedx_bg[bgBin]->Fill(dedx_new);
    else hdedx_bg[bgBin]->Fill(dedx_pre);
    hdedx_bg_pub[bgBin]->Fill(dedxpub);
    if( correct ) hchi_bg[bgBin]->Fill(chi_new);
    hchi_bg_pub[bgBin]->Fill(chiPi);

    sumcos[bgBin] += costh;
    sumsin[bgBin] += sqrt(1-costh*costh);
    sumbg[bgBin] += bg;
    sumressq[bgBin] += pow(dedx_res,2);
    sumsize[bgBin] += 1;

    // make histograms of dE/dx vs. cos(theta) for validation
    int icos = (int)((costh+1)/cosstep);
    hchi_costh_all[charge][icos]->Fill(chi_new);
    if( bgBin < int(m_bgbins/3) )
      hchi_costh_0[charge][icos]->Fill(chi_new);
    else if( bgBin < 2*int(m_bgbins/3) )
      hchi_costh_1[charge][icos]->Fill(chi_new);
    else
      hchi_costh_2[charge][icos]->Fill(chi_new);

  }// end of event loop

  // --------------------------------------------------
  // FIT IN BINS OF BETA-GAMMA
  // --------------------------------------------------

  // fit the histograms with Gaussian functions
  // and extract the means and errors
  TTree* satTree = new TTree(particle,"dE/dx means and errors");

  double satbg;          // beta-gamma value for this bin
  double satcosth;       // cos(theta) value for this bin
  double satdedx;        // mean dE/dx value for this bin
  double saterror;       // error on ^
  double satbg_avg;      // average beta-gamma value for this sample
  double satcosth_avg;   // average cos(theta) value for this sample
  double satsinth_avg;   // average sin(theta) value for this sample
  double satdedxres_avg; // average dE/dx error squared for this sample
  double satpubdedx;     // mean "public" dE/dx value for this bin
  double satpuberror;    // error on ^
  double satpubwidth;    // width of ^ distribution
  double satchi;         // mean chi value for this bin
  double satchierr;      // error on ^
  double satchiwidth;    // width of ^ distribution
  double satchipub;      // mean "public" chi value for this bin
  double satchipuberr;   // error on ^
  double satchipubwidth; // width of ^ distribution
  double ratio;          // ratio of the predicted mean to that of the average

  satTree->Branch("bg",&satbg,"bg/D");
  satTree->Branch("costh",&satcosth,"costh/D");
  satTree->Branch("dedx",&satdedx,"dedx/D");
  satTree->Branch("error",&saterror,"error/D");
  satTree->Branch("bg_avg",&satbg_avg,"bg_avg/D");
  satTree->Branch("costh_avg",&satcosth_avg,"costh_avg/D");
  satTree->Branch("sinth_avg",&satsinth_avg,"sinth_avg/D");
  satTree->Branch("dedxres_avg",&satdedxres_avg,"dedxres_avg/D");
  satTree->Branch("dedx_pub",&satpubdedx,"dedx_pub/D");
  satTree->Branch("dedxerr_pub",&satpuberror,"dedxerr_pub/D");
  satTree->Branch("chi",&satchi,"chi/D");
  satTree->Branch("chi_pub",&satchipub,"chi_pub/D");
  satTree->Branch("sigma",&satchiwidth,"sigma/D");
  satTree->Branch("sigma_pub",&satchipubwidth,"sigma_pub/D");
  satTree->Branch("ratio",&ratio,"ratio/D");

  // Fit the histograms
  for( int i = 0; i < m_bgbins; ++i ){

    // fill some details for this bin
    satbg = m_lowerbg+0.5*bgstep+i*bgstep;

    satbg_avg = sumbg[i]/sumsize[i];
    satcosth_avg = sumcos[i]/sumsize[i];
    satsinth_avg = sumsin[i]/sumsize[i];
    //    satdedxres_avg = sumressq[i]/sumsize[i];

    hsigma_bg[i]->Fit("gaus","ql");
    satdedxres_avg = hsigma_bg[i]->GetFunction("gaus")->GetParameter(2);

    int fs;

    // fit the dE/dx distribution in bins of beta-gamma
    fs = hdedx_bg[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;
    double mean = hdedx_bg[i]->GetFunction("gaus")->GetParameter(1);
    double width = hdedx_bg[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_bg[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN FIT STATUS " << fs << std::endl;

    satdedx = hdedx_bg[i]->GetFunction("gaus")->GetParameter(1);
    saterror = hdedx_bg[i]->GetFunction("gaus")->GetParError(1);
    means[i] = satdedx;
    errors[i] = saterror;
    widths[i] = hdedx_bg[i]->GetFunction("gaus")->GetParameter(2);


    // fit the dE/dx distribution without correction in bins of beta-gamma
    fs = hdedx_bg_pub[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN PUB FIT STATUS " << fs << std::endl;
    mean = hdedx_bg_pub[i]->GetFunction("gaus")->GetParameter(1);
    width = hdedx_bg_pub[i]->GetFunction("gaus")->GetParameter(2);
    fs = hdedx_bg_pub[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"MEAN PUB FIT STATUS " << fs << std::endl;

    satpubdedx = hdedx_bg_pub[i]->GetFunction("gaus")->GetParameter(1);
    satpuberror = hdedx_bg_pub[i]->GetFunction("gaus")->GetParError(1);
    satpubwidth = hdedx_bg_pub[i]->GetFunction("gaus")->GetParameter(2);


    // fit the chi distribution
    fs = hchi_bg[i]->Fit("gaus","ql");
    if( fs != 0 ) std::cout << "\t\t" <<"CHI FIT STATUS " << fs << std::endl;
    mean = hchi_bg[i]->GetFunction("gaus")->GetParameter(1);
    width = hchi_bg[i]->GetFunction("gaus")->GetParameter(2);
    fs = hchi_bg[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"CHI FIT STATUS " << fs << std::endl;

    satchi = hchi_bg[i]->GetFunction("gaus")->GetParameter(1);
    satchierr = hchi_bg[i]->GetFunction("gaus")->GetParError(1);
    satchiwidth = hchi_bg[i]->GetFunction("gaus")->GetParameter(2);


    // fit the chi distribution without correction
    /*
    fs = hchi_bg_pub[i]->Fit("gaus","ql");	    
    if( fs != 0 ) std::cout << "\t\t" <<"CHI PUB FIT STATUS " << fs << std::endl;
    mean = hchi_bg_pub[i]->GetFunction("gaus")->GetParameter(1);
    width = hchi_bg_pub[i]->GetFunction("gaus")->GetParameter(2);
    fs = hchi_bg_pub[i]->Fit("gaus","ql","",mean-2.5*width,mean+2.5*width);
    if( fs != 0 ) std::cout << "\t\t" <<"CHI PUB FIT STATUS " << fs << std::endl;

    satchipub = hchi_bg_pub[i]->GetFunction("gaus")->GetParameter(1);
    satchipuberr = hchi_bg_pub[i]->GetFunction("gaus")->GetParError(1);
    satchipubwidth = hchi_bg_pub[i]->GetFunction("gaus")->GetParameter(2);
    */

    // determine the ratio of the predicted mean at a given bg to that of the average
    ratio = m_gpar.dedxPrediction(satbg_avg)/m_gpar.dedxPrediction(satbg);

    // fill the tree for this bin
    satTree->Fill();
  }

  std::string corname("uncorrected");
  if( correct == true ) corname = "corrected";

  // Print the histograms for quality control
  TCanvas* ctmp = new TCanvas("tmp","tmp",900,900);
  ctmp->Divide(3,3);
  std::stringstream psname; psname << "plots/dedx_" << particle << "_" << corname << ".ps[";
  ctmp->Print(psname.str().c_str());
  psname.str(""); psname << "plots/dedx_" << particle << "_" << corname << ".ps";
  for( int i = 0 ; i < m_bgbins; ++i ){
    ctmp->cd(i%9+1);
    hdedx_bg[i]->Draw();
    if((i+1)%9==0)
      ctmp->Print(psname.str().c_str());
  }
  psname.str(""); psname << "plots/dedx_" << particle << "_" << corname << ".ps]";
  ctmp->Print(psname.str().c_str());
  delete ctmp;

  std::cout << "HadronPrep: saving output to " << outfile->GetName() << std::endl;

  // write out the data to file
  outfile->cd();
  satTree->Write();

  std::cout << "making validation plots..." << std::endl;

  // --------------------------------------------------
  // FIT IN BINS OF COS(THETA) FOR VALIDATION
  // --------------------------------------------------

  double chicos[2][cosbins], sigmacos[2][cosbins];
  double chicoserr[2][cosbins], sigmacoserr[2][cosbins];

  double chicos0[2][cosbins], sigmacos0[2][cosbins];
  double chicos0err[2][cosbins], sigmacos0err[2][cosbins];

  double chicos1[2][cosbins], sigmacos1[2][cosbins];
  double chicos1err[2][cosbins], sigmacos1err[2][cosbins];

  double chicos2[2][cosbins], sigmacos2[2][cosbins];
  double chicos2err[2][cosbins], sigmacos2err[2][cosbins];

  for( int c = 0; c < 2; ++c ){
    for( int i = 0; i < cosbins; ++i ){
      if( hchi_costh_all[c][i]->Integral(1,100) == 0 ) continue;
      hchi_costh_all[c][i]->Fit("gaus","ql");
      chicos[c][i] = hchi_costh_all[c][i]->GetFunction("gaus")->GetParameter(1);
      chicoserr[c][i] = hchi_costh_all[c][i]->GetFunction("gaus")->GetParError(1);
      sigmacos[c][i] = hchi_costh_all[c][i]->GetFunction("gaus")->GetParameter(2);
      sigmacoserr[c][i] = hchi_costh_all[c][i]->GetFunction("gaus")->GetParError(2);
    }
    for( int i = 0; i < cosbins; ++i ){
      if( hchi_costh_0[c][i]->Integral(1,100) == 0 ) continue;
      hchi_costh_0[c][i]->Fit("gaus","ql");
      chicos0[c][i] = hchi_costh_0[c][i]->GetFunction("gaus")->GetParameter(1);
      chicos0err[c][i] = hchi_costh_0[c][i]->GetFunction("gaus")->GetParError(1);
      sigmacos0[c][i] = hchi_costh_0[c][i]->GetFunction("gaus")->GetParameter(2);
      sigmacos0err[c][i] = hchi_costh_0[c][i]->GetFunction("gaus")->GetParError(2);
    }
    for( int i = 0; i < cosbins; ++i ){
      if( hchi_costh_1[c][i]->Integral(1,100) == 0 ) continue;
      hchi_costh_1[c][i]->Fit("gaus","ql");
      chicos1[c][i] = hchi_costh_1[c][i]->GetFunction("gaus")->GetParameter(1);
      chicos1err[c][i] = hchi_costh_1[c][i]->GetFunction("gaus")->GetParError(1);
      sigmacos1[c][i] = hchi_costh_1[c][i]->GetFunction("gaus")->GetParameter(2);
      sigmacos1err[c][i] = hchi_costh_1[c][i]->GetFunction("gaus")->GetParError(2);
    }
    for( int i = 0; i < cosbins; ++i ){
      if( hchi_costh_2[c][i]->Integral(1,100) == 0 ) continue;
      hchi_costh_2[c][i]->Fit("gaus","ql");
      chicos2[c][i] = hchi_costh_2[c][i]->GetFunction("gaus")->GetParameter(1);
      chicos2err[c][i] = hchi_costh_2[c][i]->GetFunction("gaus")->GetParError(1);
      sigmacos2[c][i] = hchi_costh_2[c][i]->GetFunction("gaus")->GetParameter(2);
      sigmacos2err[c][i] = hchi_costh_2[c][i]->GetFunction("gaus")->GetParError(2);
    }
  }

  // Print the histograms for quality control
  TCanvas* ctmp2 = new TCanvas("tmp2","tmp2",900,900);
  ctmp2->Divide(3,3);
  psname.str(""); psname << "plots/chivscos_fits_" << particle << ".ps[";
  ctmp2->Print(psname.str().c_str());
  psname.str(""); psname << "plots/chivscos_fits_" << particle << ".ps";
  for( int i = 0 ; i < cosbins; ++i ){
    ctmp2->cd(i%9+1);
    hchi_costh_all[0][i]->Draw();
    hchi_costh_all[1][i]->SetMarkerColor(kRed);
    hchi_costh_all[1][i]->Draw("same");
    if((i+1)%9==0)
      ctmp2->Print(psname.str().c_str());
  }
  psname.str(""); psname << "plots/chivscos_fits_" << particle << ".ps]";
  ctmp2->Print(psname.str().c_str());
  delete ctmp2;

  TGraphErrors* grchicos = new TGraphErrors(cosbins,cosArray,chicos[0],cosArrayErr,chicoserr[0]);
  TGraphErrors* grchicos0 = new TGraphErrors(cosbins,cosArray,chicos0[0],cosArrayErr,chicos0err[0]);
  TGraphErrors* grchicos1 = new TGraphErrors(cosbins,cosArray,chicos1[0],cosArrayErr,chicos1err[0]);
  TGraphErrors* grchicos2 = new TGraphErrors(cosbins,cosArray,chicos2[0],cosArrayErr,chicos2err[0]);

  TGraphErrors* grchicosn = new TGraphErrors(cosbins,cosArray,chicos[1],cosArrayErr,chicoserr[1]);
  TGraphErrors* grchicos0n = new TGraphErrors(cosbins,cosArray,chicos0[1],cosArrayErr,chicos0err[1]);
  TGraphErrors* grchicos1n = new TGraphErrors(cosbins,cosArray,chicos1[1],cosArrayErr,chicos1err[1]);
  TGraphErrors* grchicos2n = new TGraphErrors(cosbins,cosArray,chicos2[1],cosArrayErr,chicos2err[1]);

  TGraphErrors* grsigmacos = new TGraphErrors(cosbins,cosArray,sigmacos[0],cosArrayErr,sigmacoserr[0]);
  TGraphErrors* grsigmacos0 = new TGraphErrors(cosbins,cosArray,sigmacos0[0],cosArrayErr,sigmacos0err[0]);
  TGraphErrors* grsigmacos1 = new TGraphErrors(cosbins,cosArray,sigmacos1[0],cosArrayErr,sigmacos1err[0]);
  TGraphErrors* grsigmacos2 = new TGraphErrors(cosbins,cosArray,sigmacos2[0],cosArrayErr,sigmacos2err[0]);

  TGraphErrors* grsigmacosn = new TGraphErrors(cosbins,cosArray,sigmacos[1],cosArrayErr,sigmacoserr[1]);
  TGraphErrors* grsigmacos0n = new TGraphErrors(cosbins,cosArray,sigmacos0[1],cosArrayErr,sigmacos0err[1]);
  TGraphErrors* grsigmacos1n = new TGraphErrors(cosbins,cosArray,sigmacos1[1],cosArrayErr,sigmacos1err[1]);
  TGraphErrors* grsigmacos2n = new TGraphErrors(cosbins,cosArray,sigmacos2[1],cosArrayErr,sigmacos2err[1]);

  TLine* line0 = new TLine(-1,0,1,0);
  line0->SetLineStyle(kDashed);
  line0->SetLineColor(kRed);
  TLine* line1 = new TLine(-1,1,1,1);
  line1->SetLineStyle(kDashed);
  line1->SetLineColor(kRed);

  TString ptype("");
  if( particle == "pion" ) ptype += "#pi";
  else if( particle == "kaon" ) ptype += "K";
  else if( particle == "proton" ) ptype += "p";
  else if( particle == "muon" ) ptype += "#mu";
  else if( particle == "electron" ) ptype += "e";

  TLegend* lchi = new TLegend(0.6,0.75,0.8,0.85);
  lchi->SetBorderSize(0);
  lchi->SetFillColor(0);
  lchi->AddEntry(grchicos,ptype+"^{+}","p");   
  lchi->AddEntry(grchicosn,ptype+"^{-}","p");
  TCanvas* cchi = new TCanvas("cchi","cchi",700,600);

  FormatGraph(grchicos,0);
  FormatGraph(grchicosn,1);
  grchicos->Draw("AP");
  grchicosn->Draw("P,same");
  line0->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/chiall"+particle+".eps");

  FormatGraph(grchicos0,0);
  FormatGraph(grchicos0n,1);
  grchicos0->Draw("AP");
  grchicos0n->Draw("P,same");
  line0->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/chi0"+particle+".eps");

  FormatGraph(grchicos1,0);
  FormatGraph(grchicos1n,1);
  grchicos1->Draw("AP");
  grchicos1n->Draw("P,same");
  line0->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/chi1"+particle+".eps");

  FormatGraph(grchicos2,0);
  FormatGraph(grchicos2n,1);
  grchicos2->Draw("AP");
  grchicos2n->Draw("P,same");
  line0->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/chi2"+particle+".eps");

  FormatGraph(grsigmacos,2);
  FormatGraph(grsigmacosn,1);
  grsigmacos->Draw("AP");
  grsigmacosn->Draw("P,same");
  line1->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/sigmaall"+particle+".eps");

  FormatGraph(grsigmacos0,2);
  FormatGraph(grsigmacos0n,1);
  grsigmacos0->Draw("AP");
  grsigmacos0n->Draw("P,same");
  line1->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/sigma0"+particle+".eps");

  FormatGraph(grsigmacos1,2);
  FormatGraph(grsigmacos1n,1);
  grsigmacos1->Draw("AP");
  grsigmacos1n->Draw("P,same");
  line1->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/sigma1"+particle+".eps");

  FormatGraph(grsigmacos2,2);
  FormatGraph(grsigmacos2n,1);
  grsigmacos2->Draw("AP");
  grsigmacos2n->Draw("P,same");
  line1->Draw("same");
  lchi->Draw("same");
  cchi->SaveAs("plots/sigma2"+particle+".eps");

  for( int i = 0; i < 2; ++i ){
    for( int j = 0; j < cosbins; ++j ){
      delete hchi_costh_all[i][j];
      delete hchi_costh_0[i][j];
      delete hchi_costh_1[i][j];
      delete hchi_costh_2[i][j];
    }
  }

  delete grchicos;
  delete grchicos0;
  delete grchicos1;
  delete grchicos2;

  delete grchicosn;
  delete grchicos0n;
  delete grchicos1n;
  delete grchicos2n;

  delete grsigmacos;
  delete grsigmacos0;
  delete grsigmacos1;
  delete grsigmacos2;

  delete grsigmacosn;
  delete grsigmacos0n;
  delete grsigmacos1n;
  delete grsigmacos2n;

  delete line0;
  delete line1;

  delete lchi;
  delete cchi;
}
