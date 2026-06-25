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

void makeplot_angle()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  double twopi = 2*TMath::Pi();
  std::pair<double, double> phirange_central = {-1.73246, -1.43608};
  std::pair<double, double> phirange_west = {-1.21241, -0.909953};
  std::pair<double, double> phirange_east = {-2.26272, -1.96089};

  const int nrun = 1;
  int mbdrates[nrun] = {400};
  int runs[nrun] = {79516};

  for (int k=0; k<nrun; k++)
  {
    TChain* intree = new TChain("h_QAG4SimulationDistortions_residTree");
    for (int i=0; i<500; i++) {intree->Add(Form("../../root/Reconstructed/%d/clusters_seeds_%d-%d-0.root_qa.root",runs[k],runs[k],i));}
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
  
    TH1* h1_tanAlpha = new TH1F("h1_tanAlpha","tanAlpha; tanAlpha; Entries",100,-0.5,0.5);
    TH1* h1_tanBeta = new TH1F("h1_tanBeta","tanBeta; tanBeta; Entries",100,-1.5,1.5);
    TH1* h1_tanBeta_tpottile[8];
    for (int i=0; i<8; i++)
    {
      h1_tanBeta_tpottile[i] = new TH1F(Form("h1_tanBeta_tpottile%d",i),"tanBeta; tanBeta; Entries",100,-1.5,1.5);
    }

    int nevent = intree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
      intree->GetEntry(i);
      double erp = square(clusRPhiErr) + square(stateRPhiErr);
      double ez = square(clusZErr) + square(stateZErr);
      if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
      if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
      if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;

      h1_tanAlpha->Fill(tanAlpha);
      h1_tanBeta->Fill(tanBeta);
      for (int j=0; j<8; j++)
      {
        if (tileid_tpot->at(0)==j)
	{
	  h1_tanBeta_tpottile[j]->Fill(tanBeta);
	}
      }
    }
  
    TCanvas *can_tanAlpha = new TCanvas("can_tanAlpha","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_tanAlpha->Draw();
    can_tanAlpha->Update();
    can_tanAlpha->SaveAs(Form("figure/%d_tanAlpha.pdf",runs[k]));
 
    TCanvas *can_tanBeta = new TCanvas("can_tanBeta","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_tanBeta->Draw();
    can_tanBeta->Update();
    can_tanBeta->SaveAs(Form("figure/%d_tanBeta.pdf",runs[k]));

    TCanvas *can_tanBeta_tpottile = new TCanvas("can_tanBeta_tpottile","",1600*4,1200*2);
    can_tanBeta_tpottile->Divide(4,2);
    for (int i=0; i<8; i++)
    {
      can_tanBeta_tpottile->cd(i+1);
      gPad->SetTopMargin(0.20);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.15);
      h1_tanBeta_tpottile[i]->Draw();

      TPaveText *pt = new TPaveText(.05, .90, .95, 1.00, "NDC");
      pt->SetFillColor(0);
      //pt->SetFillStyle(0);//transparent
      pt->SetLineColor(0);
      pt->SetBorderSize(0);
      pt->SetTextColor(kBlack);
      pt->AddText(Form("TPOT tile ID %d",i));
      pt->Draw("same");
    }
    can_tanBeta_tpottile->Update();
    can_tanBeta_tpottile->SaveAs(Form("figure/%d_tanBeta_tpottile.pdf",runs[k]));

  }
}
