#include "WidgetQuality.h"

WidgetQuality::WidgetQuality(){
  m_type = 0;
  m_mass[0] = Widget::mpion;
  m_mass[1] = Widget::mkaon;
  m_mass[2] = Widget::mproton;
  m_mass[3] = Widget::mmuon;
  m_mass[4] = Widget::melectron;
  for( int i = 0; i < 5; ++i ){
    m_upperp[i] = 3.0;
    m_lowerp[i] = 0.0;
  }
}

WidgetQuality::WidgetQuality( int type, double* bgmax, double* bgmin ){
  m_type = type;
  m_mass[0] = Widget::mpion;
  m_mass[1] = Widget::mkaon;
  m_mass[2] = Widget::mproton;
  m_mass[3] = Widget::mmuon;
  m_mass[4] = Widget::melectron;
  for( int i = 0; i < 5; ++i ){
    m_upperp[i] = bgmax[i]*m_mass[i];
    m_lowerp[i] = bgmin[i]*m_mass[i];
  }
}

void
WidgetQuality::makeHistograms( std::vector<TString> filenames, TString outfilename ){

  std::cout << "Making histograms..." << std::endl;
  
  // check if the output file exists, do not overwrite it by default
  std::string overwrite;
  struct stat buffer;
  if( stat (outfilename.Data(), &buffer) == 0 ){
    TFile* testfile = new TFile(outfilename);
    if( testfile && testfile->GetListOfKeys()->IsEmpty() == false ){
      std::cout << outfilename << " exists. Are you sure you want to overwrite it (y/n)? ";
      std::cin >> overwrite;
      if( overwrite == "n" ) return;
      testfile->Close();
    }
  }

  // create the output file
  TFile* outfile = new TFile(outfilename,"RECREATE");

  std::cout << "Making histograms..." << std::endl;
    
  // for each file, create and save the various quality histograms
  for( int i = 0; i < filenames.size(); ++i ){
    TFile* infile = new TFile(filenames[i]);
    TTree* intree = (TTree*)infile->Get("track");
    infile->cd();

    TString cuts = TString::Format("dedxsat>0&&costh==costh");
    
    if( m_type == 0 )
      cuts += TString::Format("&&numGoodHits>=0&&numGoodHits<=100");
    if( m_type == 1 )
      cuts += TString::Format("&&lNHitsUsed>=0&&lNHitsUsed<=100");
      //      cuts += TString::Format("&&numGoodLayerHits>=0&&numGoodLayerHits<=100");

    TString name;
    if( i == 0 ){
      name = "pion";
      cuts += "&&eopst<0.75";
    }
    else if( i == 1 ){
      name = "kaon";
      cuts += "&&dedx>=0.4/TMath::Abs(pF)";
    }
    else if( i == 2 ){
      name = "proton";
      cuts += "&&dedx>0&&dedx<100&&dedx>0.85/TMath::Abs(pF)";
    }
    else if( i == 3 ) name = "muon";
    else if( i == 4 ) name = "electron";
  
    // make a histogram of momentum
    intree->Project(TString::Format("hp(100,0,%f)",m_upperp[i]),"abs(pF)",cuts);
    TH1F* hp = (TH1F*)infile->Get("hp");
    hp->SetName(name+"_hp");
    hp->SetTitle(""); hp->SetStats(0);
    hp->SetXTitle("p  [GeV/c]"); hp->SetYTitle("Events");

    // make a histogram of cos(theta)
    intree->Project("hcosth(100,-1,1)","costh",cuts);
    TH1F* hcosth = (TH1F*)infile->Get("hcosth");
    hcosth->SetName(name+"_hcosth");
    hcosth->SetTitle(""); hcosth->SetStats(0);
    hcosth->SetXTitle("cos(#theta)"); hcosth->SetYTitle("Events");

    // make a histogram of dE/dx vs. momentum
    if( i == 1 || i == 2 ) intree->Project(TString::Format("dedx_p(100,0,%f,100,0,10)",m_upperp[i]+(m_upperp[i]-m_lowerp[i])*0.1),"dedx:abs(pF)",cuts);
    else intree->Project(TString::Format("dedx_p(100,0,%f,100,0,3)",m_upperp[i]+(m_upperp[i]-m_lowerp[i])*0.1),"dedx:abs(pF)",cuts);
    TH2F* dedx_p = (TH2F*)infile->Get("dedx_p");
    dedx_p->SetName(name+"_dedx_p");
    dedx_p->SetTitle(""); dedx_p->SetStats(0);
    dedx_p->SetXTitle("p  [GeV/c]"); dedx_p->SetYTitle("dE/dx");

    // make a histogram of momentum vs. cos(theta)
    intree->Project(TString::Format("p_costh(100,-1,1,100,0,%f)",m_upperp[i]+(m_upperp[i]-m_lowerp[i])*0.1),"abs(pF):costh",cuts);
    TH2F* p_costh = (TH2F*)infile->Get("p_costh");
    p_costh->SetName(name+"_p_costh");
    p_costh->SetTitle(""); p_costh->SetStats(0);
    p_costh->SetYTitle("p  [GeV/c]"); p_costh->SetXTitle("cos(#theta)");

    // make a histogram of nhit vs. pt
    if( m_type == 1 )
      intree->Project("nhit_pt(100,0,2.0,60,0,60)","lNHitsUsed:abs(pF)*sqrt(1-costh*costh)",cuts);
      //      intree->Project("nhit_pt(100,0,2.0,40,0,40)","numGoodLayerHits:abs(pF)*sqrt(1-costh*costh)",cuts);
    else
      intree->Project("nhit_pt(100,0,2.0,60,0,60)","numGoodHits:abs(pF)*sqrt(1-costh*costh)",cuts);

    TH2F* nhit_pt = (TH2F*)infile->Get("nhit_pt");
    nhit_pt->SetName(name+"_nhit_pt");
    nhit_pt->SetTitle(""); nhit_pt->SetStats(0);
    nhit_pt->SetXTitle("p_t  [GeV/c]"); nhit_pt->SetYTitle("Number of Hits");
    
    // save it to the output file
    outfile->cd();
    hp->Write();
    hcosth->Write();
    dedx_p->Write();
    p_costh->Write();
    nhit_pt->Write();

    infile->Close();
  }

  outfile->Close();
}

void
WidgetQuality::plotHistograms( TString filename ){

  std::cout << "Plotting histograms..." << std::endl;
    
  TFile* f = new TFile(filename);

  // five different histograms are plotted for each particle type
  TString varname[5];
  varname[0] = "_hp";
  varname[1] = "_hcosth";
  varname[2] = "_dedx_p";
  varname[3] = "_p_costh";
  varname[4] = "_nhit_pt";

  TCanvas* can = new TCanvas("can","Quality plots",800,800);
  for( int i = 0; i < 5; ++i ){

    TString name;
    if( i == 0 ) name = "pion";
    else if( i == 1 ) name = "kaon";
    else if( i == 2 ) name = "proton";
    else if( i == 3 ) name = "muon";
    else if( i == 4 ) name = "electron";

    // plot 5 histograms for each particle type
    for( int j = 0; j < 5; ++j ){
      TString histname = name+varname[j];
      TH1F* h1 = (TH1F*)f->Get(histname);

      h1->DrawCopy("colz");
      if( j == 2 ){
	// draw lines that define the momentum range
	double max = 3.0; if( i == 1 || i == 2 ) max = 10.0;
	TLine* upper = new TLine(m_upperp[i],0,m_upperp[i],max);
	TLine* lower = new TLine(m_lowerp[i],0,m_lowerp[i],max);
	upper->SetLineColor(kRed); upper->SetLineStyle(kDashed);
	lower->SetLineColor(kRed); lower->SetLineStyle(kDashed);
	upper->DrawClone("same");
	lower->DrawClone("same");
	delete upper;
	delete lower;
      }
      can->SaveAs("plots/"+histname+".eps");
    }
  }

  f->Close();
  delete can;
}
