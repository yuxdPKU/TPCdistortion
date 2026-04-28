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
    TH1* h1_drphi_tpot = new TH1F("h1_drphi_tpot",";Layer 55 #delta r#phi (clus-state) (cm); Entries",100,-0.5,0.5);
    TH1* h1_dz_tpot = new TH1F("h1_dz_tpot",";Layer 56 #delta z (clus-state) (cm); Entries",100,-5,5);
    h1_drphi->SetLineColor(kBlack);
    h1_dz->SetLineColor(kBlack);
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

      if (event!=event_last)
      {
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

    TCanvas *can_drphi_tpot = new TCanvas("can_drphi_tpot","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_drphi_tpot->Draw("hist");
    TArrow* arrow_drphi_tpot_left = new TArrow(-0.1, 0.5*(h1_drphi_tpot->GetMaximum()), -0.1, 0, 0.03, "->");
    arrow_drphi_tpot_left->SetLineColor(kRed);
    arrow_drphi_tpot_left->Draw();
    TArrow* arrow_drphi_tpot_right = new TArrow(0.1, 0.5*(h1_drphi_tpot->GetMaximum()), 0.1, 0, 0.03, "->");
    arrow_drphi_tpot_right->SetLineColor(kRed);
    arrow_drphi_tpot_right->Draw();
    can_drphi_tpot->Update();
    can_drphi_tpot->SaveAs(Form("figure/%d_drphi_tpot.pdf",runs[k]));

    TCanvas *can_dz_tpot = new TCanvas("can_dz_tpot","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h1_dz_tpot->Draw("hist");
    TArrow* arrow_dz_tpot_left = new TArrow(-1, 0.5*(h1_dz_tpot->GetMaximum()), -1, 0, 0.03, "->");
    arrow_dz_tpot_left->SetLineColor(kRed);
    arrow_dz_tpot_left->Draw();
    TArrow* arrow_dz_tpot_right = new TArrow(1, 0.5*(h1_dz_tpot->GetMaximum()), 1, 0, 0.03, "->");
    arrow_dz_tpot_right->SetLineColor(kRed);
    arrow_dz_tpot_right->Draw();
    can_dz_tpot->Update();
    can_dz_tpot->SaveAs(Form("figure/%d_dz_tpot.pdf",runs[k]));

  }
}
