void fit_gauss(TH2* h, TString name, bool verbose=0);

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

void ntpcstates()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 1;
  //int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
  //int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
  //int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
  //int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
  int mbdrates[nrun] = {400};
  int runs[nrun] = {53877};

  for (int k=0; k<nrun; k++)
  {
    TFile* infile = new TFile(Form("allqa_%d.root",runs[k]));
    TTree* intree = (TTree*) infile->Get("h_QAG4SimulationDistortions_residTree");
    float drphi, clusR;
    int track_nmvtx, track_nintt, track_ntpc, track_ntpot;
    int track_nmvtxstate, track_ninttstate, track_ntpcstate, track_ntpotstate;
    int event;
    intree->SetBranchAddress("clusR",&clusR);
    intree->SetBranchAddress("drphi",&drphi);
    intree->SetBranchAddress("track_nmvtx",&track_nmvtx);
    intree->SetBranchAddress("track_nintt",&track_nintt);
    intree->SetBranchAddress("track_ntpc",&track_ntpc);
    intree->SetBranchAddress("track_ntpot",&track_ntpot);
    intree->SetBranchAddress("track_nmvtxstate",&track_nmvtxstate);
    intree->SetBranchAddress("track_ninttstate",&track_ninttstate);
    intree->SetBranchAddress("track_ntpcstate",&track_ntpcstate);
    intree->SetBranchAddress("track_ntpotstate",&track_ntpotstate);
    intree->SetBranchAddress("event",&event);
  
    TH2* h2_tpc_clus_state = new TH2F("h2_tpc_clus_state","TPC nclus vs nstate;nstates;nclus",50,0,50,50,0,50);
  
    int nevent = intree->GetEntries();
    int eventbefore = 0;
    for (int i=0; i<nevent; i++)
    {
      intree->GetEntry(i);
      if (event!=eventbefore)
      {
	if (!(clusR>50 && drphi<=0.5)) continue;
        //if (track_ntpcstate>30) continue;
        //if ((track_ntpc-track_ntpcstate)>5) continue;
	//if (!(track_ntpc>45 && track_ntpcstate>28 && track_ntpcstate<39)) continue;
        h2_tpc_clus_state->Fill(track_ntpcstate,track_ntpc);
	eventbefore=event;
      }
    }
  
    TCanvas *can = new TCanvas("can","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h2_tpc_clus_state->Draw("colz");
    can->Update();
    can->SaveAs(Form("figure/%d_tpc_clus_state.pdf",runs[k]));
   
  }
}
