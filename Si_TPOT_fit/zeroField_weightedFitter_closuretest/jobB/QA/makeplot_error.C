void fit_gauss_Rslice(TH2* h, TString name, bool verbose=0);
double fit_gauss_Zslice(TH2* h, TString name, int verbose=0);

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

void makeplot_error()
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
    std::vector<int> *tileid_tpot = 0;
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
    intree->SetBranchAddress("charge",&charge);
    intree->SetBranchAddress("tileid_tpot", &tileid_tpot);
  
    TH1* h1_cluster_zerr = new TH1F("h1_cluster_zerr","; Zerr (cm); Entries",100,0,0.08);
    TH1* h1_state_zerr = new TH1F("h1_state_zerr","; Zerr (cm); Entries",100,0,0.08);
    TH1* h1_cluster_rphierr = new TH1F("h1_cluster_rphierr","; RPhierr (cm); Entries",100,0,0.05);
    TH1* h1_state_rphierr = new TH1F("h1_state_rphierr","; RPhierr (cm); Entries",100,0,0.05);
    h1_cluster_zerr->SetLineColor(kBlack);
    h1_cluster_rphierr->SetLineColor(kBlack);
    h1_cluster_zerr->SetLineColor(kRed);
    h1_cluster_rphierr->SetLineColor(kRed);

    int nevent = intree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
      intree->GetEntry(i);
      h1_cluster_zerr->Fill(clusZErr);
      h1_cluster_rphierr->Fill(clusRPhiErr);
      h1_state_zerr->Fill(stateZErr);
      h1_state_rphierr->Fill(stateRPhiErr);
    }
  
    TCanvas *can_zerr = new TCanvas("can_zerr","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_cluster_zerr->Scale(1/(h1_cluster_zerr->Integral()));
    h1_state_zerr->Scale(1/(h1_state_zerr->Integral()));
    h1_state_zerr->Draw("hist");
    h1_cluster_zerr->Draw("hist,same");
    TLegend *legend_zerr = new TLegend(0.50, 0.65, 0.80, 0.9);
    legend_zerr->AddEntry(h1_cluster_zerr, "Cluster", "F");
    legend_zerr->AddEntry(h1_state_zerr, "State", "F");
    //legend_zerr->SetTextSize(0.02);
    legend_zerr->Draw("same");
    can_zerr->Update();
    can_zerr->SaveAs(Form("figure/%d_zerr.pdf",runs[k]));

    TCanvas *can_rphierr = new TCanvas("can_rphierr","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_cluster_rphierr->Scale(1/(h1_cluster_rphierr->Integral()));
    h1_state_rphierr->Scale(1/(h1_state_rphierr->Integral()));
    h1_state_rphierr->Draw("hist");
    h1_cluster_rphierr->Draw("hist,same");
    TLegend *legend_rphierr = new TLegend(0.50, 0.65, 0.80, 0.9);
    legend_rphierr->AddEntry(h1_cluster_rphierr, "Cluster", "F");
    legend_rphierr->AddEntry(h1_state_rphierr, "State", "F");
    //legend_rphierr->SetTextSize(0.02);
    legend_rphierr->Draw("same");
    can_rphierr->Update();
    can_rphierr->SaveAs(Form("figure/%d_rphierr.pdf",runs[k]));

  }
}
