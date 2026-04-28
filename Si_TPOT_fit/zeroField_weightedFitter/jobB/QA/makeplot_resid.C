void fit_gauss_Rslice(TH2* h, double deltax, TString name, bool verbose=0);

namespace
{

  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  template <class T>
  inline constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    else if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }
    else
    {
      return phi;
    }
  }

}

float getPhi(float phi)
{
  float m_phiMin = 0;
  float m_phiMax = 2. * M_PI;
  while (phi < m_phiMin)
  {
    phi += 2. * M_PI;
  }
  while (phi >= m_phiMax)
  {
    phi -= 2. * M_PI;
  }
  return phi;
}

void makeplot_resid()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  double twopi = 2*TMath::Pi();
  std::pair<double, double> phirange_central = {-1.73246, -1.43608};
  std::pair<double, double> phirange_west = {-1.21241, -0.909953};
  std::pair<double, double> phirange_east = {-2.26272, -1.96089};

  const int nrun = 1;
  int mbdrates[nrun] = {400};
  int runs[nrun] = {52077};

  for (int k=0; k<nrun; k++)
  {
    TChain* intree = new TChain("h_QAG4SimulationDistortions_residTree");
    for (int i=0; i<200; i++) {intree->Add(Form("../../root/Reconstructed/%d/clusters_seeds_%d-%d-0.root_qa.root",runs[k],runs[k],i));}
    float drphi, dz, clusZErr, stateZErr, clusZ, clusR, clusRPhiErr, stateRPhiErr, clusPhi;
    float tanAlpha, tanBeta;
    float trackPt;
    int layer, charge;
    int event;
    std::vector<int> *layer_mvtx = 0;
    std::vector<float> *sclusgr_mvtx = 0;
    std::vector<float> *sclusgz_mvtx = 0;
    std::vector<float> *sclusphi_mvtx = 0;
    std::vector<float> *stategr_mvtx = 0;
    std::vector<float> *stategz_mvtx = 0;
    std::vector<float> *statephi_mvtx = 0;
    std::vector<int> *layer_intt = 0;
    std::vector<float> *sclusgr_intt = 0;
    std::vector<float> *sclusgz_intt = 0;
    std::vector<float> *sclusphi_intt = 0;
    std::vector<float> *stategr_intt = 0;
    std::vector<float> *stategz_intt = 0;
    std::vector<float> *statephi_intt = 0;
    std::vector<int> *tileid_tpot = 0;
    std::vector<int> *layer_tpot = 0;
    std::vector<float> *sclusgr_tpot = 0;
    std::vector<float> *sclusgz_tpot = 0;
    std::vector<float> *sclusphi_tpot = 0;
    std::vector<float> *stategr_tpot = 0;
    std::vector<float> *stategz_tpot = 0;
    std::vector<float> *statephi_tpot = 0;
 
    intree->SetBranchAddress("dz",&dz);
    intree->SetBranchAddress("drphi",&drphi);
    intree->SetBranchAddress("clusZErr",&clusZErr);
    intree->SetBranchAddress("stateZErr",&stateZErr);
    intree->SetBranchAddress("clusZ",&clusZ);
    intree->SetBranchAddress("clusR",&clusR);
    intree->SetBranchAddress("clusPhi",&clusPhi);
    intree->SetBranchAddress("clusRPhiErr",&clusRPhiErr);
    intree->SetBranchAddress("stateRPhiErr",&stateRPhiErr);
    intree->SetBranchAddress("tanAlpha",&tanAlpha);
    intree->SetBranchAddress("trackPt",&trackPt);
    intree->SetBranchAddress("tanBeta",&tanBeta);
    intree->SetBranchAddress("layer",&layer);
    intree->SetBranchAddress("event",&event);
    intree->SetBranchAddress("charge",&charge);
    intree->SetBranchAddress("layer_mvtx", &layer_mvtx);
    intree->SetBranchAddress("sclusgr_mvtx", &sclusgr_mvtx);
    intree->SetBranchAddress("sclusgz_mvtx", &sclusgz_mvtx);
    intree->SetBranchAddress("sclusphi_mvtx", &sclusphi_mvtx);
    intree->SetBranchAddress("stategr_mvtx", &stategr_mvtx);
    intree->SetBranchAddress("stategz_mvtx", &stategz_mvtx);
    intree->SetBranchAddress("statephi_mvtx", &statephi_mvtx);
    intree->SetBranchAddress("layer_intt", &layer_intt);
    intree->SetBranchAddress("sclusgr_intt", &sclusgr_intt);
    intree->SetBranchAddress("sclusgz_intt", &sclusgz_intt);
    intree->SetBranchAddress("sclusphi_intt", &sclusphi_intt);
    intree->SetBranchAddress("stategr_intt", &stategr_intt);
    intree->SetBranchAddress("stategz_intt", &stategz_intt);
    intree->SetBranchAddress("statephi_intt", &statephi_intt);
    intree->SetBranchAddress("tileid_tpot", &tileid_tpot);
    intree->SetBranchAddress("layer_tpot", &layer_tpot);
    intree->SetBranchAddress("sclusgr_tpot", &sclusgr_tpot);
    intree->SetBranchAddress("sclusgz_tpot", &sclusgz_tpot);
    intree->SetBranchAddress("sclusphi_tpot", &sclusphi_tpot);
    intree->SetBranchAddress("stategr_tpot", &stategr_tpot);
    intree->SetBranchAddress("stategz_tpot", &stategz_tpot);
    intree->SetBranchAddress("statephi_tpot", &statephi_tpot);
 
    TH1* h1_drphi = new TH1F("h1_drphi","; #delta r#phi (clus-state) (cm); Entries",100,-5,5);
    TH1* h1_dz = new TH1F("h1_dz","; #delta z (clus-state) (cm); Entries",100,-8,8);
    TH2* h2_drphi_layer = new TH2F("h2_drphi_layer","TPC RPhi residual vs layer; Layer; #delta r#phi (clus-state) (cm)",60,0,60,100,-1,1);
    TH2* h2_dz_layer = new TH2F("h2_dz_layer","TPC Z residual vs layer; Layer; #delta z (clus-state) (cm)",60,0,60,100,-2,2);
    TH1* h1_drphi_mvtx = new TH1F("h1_drphi_mvtx",";#delta r#phi (clus-state) (cm); Entries",100,-0.05,0.05);
    TH1* h1_dz_mvtx = new TH1F("h1_dz_mvtx",";#delta z (clus-state) (cm); Entries",100,-0.02,0.02);
    TH1* h1_drphi_intt = new TH1F("h1_drphi_intt",";#delta r#phi (clus-state) (cm); Entries",100,-0.2,0.2);
    TH1* h1_dz_intt = new TH1F("h1_dz_intt",";#delta z (clus-state) (cm); Entries",100,-1,1);
    TH1* h1_drphi_tpot = new TH1F("h1_drphi_tpot",";Layer 55 #delta r#phi (clus-state) (cm); Entries",100,-1,1);
    TH1* h1_dz_tpot = new TH1F("h1_dz_tpot",";Layer 56 #delta z (clus-state) (cm); Entries",100,-5,5);
    h1_drphi->SetLineColor(kBlack);
    h1_dz->SetLineColor(kBlack);
    h1_drphi_mvtx->SetLineColor(kBlack);
    h1_dz_mvtx->SetLineColor(kBlack);
    h1_drphi_intt->SetLineColor(kBlack);
    h1_dz_intt->SetLineColor(kBlack);
    h1_drphi_tpot->SetLineColor(kBlack);
    h1_dz_tpot->SetLineColor(kBlack);

    int nevent = intree->GetEntries();
    int event_last = -1;
    for (int i=0; i<nevent; i++)
    {
      if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
      intree->GetEntry(i);
      h1_drphi->Fill(drphi);
      h1_dz->Fill(dz);
      h2_drphi_layer->Fill(layer,drphi);
      h2_dz_layer->Fill(layer,dz);

      if (event!=event_last)
      {
	      for (int imvtx=0; imvtx<(sclusgr_mvtx->size()); imvtx++)
	      {
          h1_drphi_mvtx->Fill(sclusgr_mvtx->at(imvtx) * (sclusphi_mvtx->at(imvtx) - statephi_mvtx->at(imvtx)));
      	  h1_dz_mvtx->Fill(sclusgz_mvtx->at(imvtx) - stategz_mvtx->at(imvtx));
      	}
	      for (int iintt=0; iintt<(sclusgr_intt->size()); iintt++)
	      {
          h1_drphi_intt->Fill(sclusgr_intt->at(iintt) * (sclusphi_intt->at(iintt) - statephi_intt->at(iintt)));
      	  h1_dz_intt->Fill(sclusgz_intt->at(iintt) - stategz_intt->at(iintt));
      	}
	      for (int itpot=0; itpot<(sclusgr_tpot->size()); itpot++)
	      {
          if (layer_tpot->at(itpot)==55) {h1_drphi_tpot->Fill(sclusgr_tpot->at(itpot) * (sclusphi_tpot->at(itpot) - statephi_tpot->at(itpot)));}
      	  if (layer_tpot->at(itpot)==56) {h1_dz_tpot->Fill(sclusgz_tpot->at(itpot) - stategz_tpot->at(itpot));}
      	}
        event_last = event;
      }
    }

    TCanvas *can_drphi = new TCanvas("can_drphi","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_drphi->Draw("hist");
    TArrow* arrow_drphi_left = new TArrow(-2, 0.5*(h1_drphi->GetMaximum()), -2, 0, 0.03, "->");
    arrow_drphi_left->SetLineColor(kRed);
    arrow_drphi_left->Draw();
    TArrow* arrow_drphi_right = new TArrow(2, 0.5*(h1_drphi->GetMaximum()), 2, 0, 0.03, "->");
    arrow_drphi_right->SetLineColor(kRed);
    arrow_drphi_right->Draw();
    can_drphi->Update();
    can_drphi->SaveAs(Form("figure/%d_drphi.pdf",runs[k]));

    TCanvas *can_dz = new TCanvas("can_dz","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_dz->Draw("hist");
    TArrow* arrow_dz_left = new TArrow(-5, 0.5*(h1_dz->GetMaximum()), -5, 0, 0.03, "->");
    arrow_dz_left->SetLineColor(kRed);
    arrow_dz_left->Draw();
    TArrow* arrow_dz_right = new TArrow(5, 0.5*(h1_dz->GetMaximum()), 5, 0, 0.03, "->");
    arrow_dz_right->SetLineColor(kRed);
    arrow_dz_right->Draw();
    can_dz->Update();
    can_dz->SaveAs(Form("figure/%d_dz.pdf",runs[k]));

    /*
    TCanvas *can_drphi_layer = new TCanvas("can_drphi_layer","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h2_drphi_layer->Draw("colz");
    can_drphi_layer->Update();
    can_drphi_layer->SaveAs(Form("figure/%d_drphi_layer.pdf",runs[k]));

    TCanvas *can_dz_layer = new TCanvas("can_dz_layer","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h2_dz_layer->Draw("colz");
    can_dz_layer->Update();
    can_dz_layer->SaveAs(Form("figure/%d_dz_layer.pdf",runs[k]));
    */

    fit_gauss_Rslice(h2_drphi_layer, 0.2, Form("figure/%d_drphi_layer.pdf",runs[k]),1); 
    fit_gauss_Rslice(h2_dz_layer, 0.5, Form("figure/%d_dz_layer.pdf",runs[k]),1); 

    TCanvas *can_drphi_mvtx = new TCanvas("can_drphi_mvtx","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_drphi_mvtx->Draw("hist");
    can_drphi_mvtx->Update();
    can_drphi_mvtx->SaveAs(Form("figure/%d_drphi_mvtx.pdf",runs[k]));

    TCanvas *can_dz_mvtx = new TCanvas("can_dz_mvtx","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_dz_mvtx->Draw("hist");
    can_dz_mvtx->Update();
    can_dz_mvtx->SaveAs(Form("figure/%d_dz_mvtx.pdf",runs[k]));

    TCanvas *can_drphi_intt = new TCanvas("can_drphi_intt","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_drphi_intt->Draw("hist");
    can_drphi_intt->Update();
    can_drphi_intt->SaveAs(Form("figure/%d_drphi_intt.pdf",runs[k]));

    TCanvas *can_dz_intt = new TCanvas("can_dz_intt","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_dz_intt->Draw("hist");
    can_dz_intt->Update();
    can_dz_intt->SaveAs(Form("figure/%d_dz_intt.pdf",runs[k]));

    TCanvas *can_drphi_tpot = new TCanvas("can_drphi_tpot","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_drphi_tpot->Draw("hist");
    TArrow* arrow_drphi_tpot_left = new TArrow(-0.5, 0.5*(h1_drphi_tpot->GetMaximum()), -0.5, 0, 0.03, "->");
    arrow_drphi_tpot_left->SetLineColor(kRed);
    arrow_drphi_tpot_left->Draw();
    TArrow* arrow_drphi_tpot_right = new TArrow(0.5, 0.5*(h1_drphi_tpot->GetMaximum()), 0.5, 0, 0.03, "->");
    arrow_drphi_tpot_right->SetLineColor(kRed);
    arrow_drphi_tpot_right->Draw();
    can_drphi_tpot->Update();
    can_drphi_tpot->SaveAs(Form("figure/%d_drphi_tpot.pdf",runs[k]));

    TCanvas *can_dz_tpot = new TCanvas("can_dz_tpot","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_dz_tpot->Draw("hist");
    TArrow* arrow_dz_tpot_left = new TArrow(-2, 0.5*(h1_dz_tpot->GetMaximum()), -2, 0, 0.03, "->");
    arrow_dz_tpot_left->SetLineColor(kRed);
    arrow_dz_tpot_left->Draw();
    TArrow* arrow_dz_tpot_right = new TArrow(2, 0.5*(h1_dz_tpot->GetMaximum()), 2, 0, 0.03, "->");
    arrow_dz_tpot_right->SetLineColor(kRed);
    arrow_dz_tpot_right->Draw();
    can_dz_tpot->Update();
    can_dz_tpot->SaveAs(Form("figure/%d_dz_tpot.pdf",runs[k]));

  }
}

void fit_gauss_Rslice(TH2* h, double deltax, TString name, bool verbose=0)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 0;
  for (int i = 1; i <= h->GetNbinsX(); i+=(n+1))
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d Layer %d",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-deltax;
    double fitmax=maxBinCenter+deltax;
    TF1 *gausFit_sub = new TF1("gausFit_sub", "gaus", fitmin, fitmax);
    projection->Fit(gausFit_sub, "Q", "", fitmin, fitmax);

    if (verbose>0)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      can->SetTopMargin(0.12);
      projection->Draw("hist");
      gausFit_full->SetLineColor(kRed);
      gausFit_full->Draw("same");
      gausFit_sub->SetLineColor(kBlue);
      gausFit_sub->Draw("same");
      myText(0.2,0.95,kBlack,Form("%s , Bin %d Layer %d",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));
      can->Update();
      can->SaveAs(name + Form(".gausFit_%d.pdf",i));
      delete can;
    }

    Double_t mean_full = gausFit_full->GetParameter(1);
    Double_t error_full = gausFit_full->GetParError(1);

    Double_t mean_sub = gausFit_sub->GetParameter(1);
    Double_t error_sub = gausFit_sub->GetParError(1);

    if (verbose>0)
    {
      cout<<"Bin "<<i<<": full fit mean = "<<mean_full<<" +/- "<<error_full<<" , sub fit mean = "<<mean_sub<<" +/- "<<error_sub<<endl;
    }

    double center = 0;
    double width = 0;
    for (int j=0; j<(n+1); j++)
    {
      center += h->GetXaxis()->GetBinCenter(i+j) / (n+1);
      width += h->GetXaxis()->GetBinWidth(i+j) / 2;
    }
    graph_full->SetPoint(i-1, center, mean_full);
    graph_full->SetPointError(i-1, width, error_full);
    graph_sub->SetPoint(i-1, center, mean_sub);
    graph_sub->SetPointError(i-1, width, error_sub);
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  c1->SetTopMargin(0.12);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  gPad->SetLogz(1);
  h->Draw("COLZ");

  //graph_full->SetMarkerStyle(20);
  //graph_full->SetMarkerColor(kRed);
  //graph_full->Draw("P,same");
  graph_sub->SetMarkerStyle(20);
  graph_sub->SetMarkerColor(kBlack);
  graph_sub->Draw("P,same");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  myText(0.4,0.95,kBlack,h->GetTitle());

  /*
  TFile* ofile = new TFile("hist_rdphi_vs_r.root","recreate");
  ofile->cd();
  graph_full->Write("graph_full");
  graph_sub->Write("graph_sub");
  h->Write();
  ofile->Close();
  */

  c1->Update();
  c1->SaveAs(name);
  delete c1;

}
