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

void debug()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TFile* infile = new TFile("allqa_53534.root");
  TTree* intree = (TTree*) infile->Get("h_QAG4SimulationDistortions_residTree");
  float drphi, dz, clusZErr, stateZErr, clusZ, clusR, clusRPhiErr, stateRPhiErr, clusPhi;
  float tanAlpha, tanBeta;
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
  intree->SetBranchAddress("tanBeta",&tanBeta);
  int nevent = intree->GetEntries();
  TH1* h1_dz = new TH1F("h1_dz","", 40, 0, 105.5);
  h1_dz->GetXaxis()->SetTitle("clusZ (cm)");
  h1_dz->GetYaxis()->SetTitle("#Sum_{cl,tr} dz/#sigma_{z}^{2}");
  TH1* h1_drphi = new TH1F("h1_drphi","", 40, 0, 105.5);
  h1_drphi->GetXaxis()->SetTitle("clusZ (cm)");
  h1_drphi->GetYaxis()->SetTitle("#Sum_{cl,tr} rd#phi/#sigma_{r#phi}^{2}");
  TH1* h1_dr = new TH1F("h1_dr","", 40, 0, 105.5);
  h1_dr->GetXaxis()->SetTitle("clusZ (cm)");
  h1_dr->GetYaxis()->SetTitle("#Sum_{cl,tr} rd#phi#tan#alpha/#sigma_{r#phi}^{2}+dz#tan#beta/#sigma_{z}^{2}");
  for (int i=0; i<40; i++)
  {
    h1_dz->SetBinContent(i+1,0);
    h1_drphi->SetBinContent(i+1,0);
    h1_dr->SetBinContent(i+1,0);
  }
  for (int i=0; i<nevent; i++)
  {
    intree->GetEntry(i);
    if (getPhi(clusPhi)<4.887 || getPhi(clusPhi)>5.061) continue;
    if (clusR<66.3 || clusR>71.6) continue;
    int ibin = h1_dz->GetXaxis()->FindBin(clusZ);
    float Dz = (dz)/(clusZErr*clusZErr+stateZErr*stateZErr);
    float Drphi = (drphi)/(clusRPhiErr*clusRPhiErr+stateRPhiErr*stateRPhiErr);
    float Dr = (drphi*tanAlpha)/(clusRPhiErr*clusRPhiErr+stateRPhiErr*stateRPhiErr)+(dz*tanBeta)/(clusZErr*clusZErr+stateZErr*stateZErr);
    if (isnan(Dz) || isnan(Drphi) || isnan(Dr))
    {
	    //cout<<"ientry = "<<i<<" , dz = "<<dz<<" , clusZErr = "<<clusZErr<<" , stateZErr = "<<stateZErr<<endl;
	    //cout<<"ientry = "<<i<<" , drphi = "<<drphi<<" , clusRPhiErr = "<<clusRPhiErr<<" , stateRPhiErr = "<<stateRPhiErr<<endl;
	    continue;
    }
    h1_dz->SetBinContent(ibin, h1_dz->GetBinContent(ibin)+Dz);
    h1_drphi->SetBinContent(ibin, h1_drphi->GetBinContent(ibin)+Drphi);
    h1_dr->SetBinContent(ibin, h1_dr->GetBinContent(ibin)+Dr);
  }
  for (int i=0; i<40; i++)
  {
    //cout<<"bin "<<i<<": bincontent "<<h1_dz->GetBinContent(i+1)<<endl;
  }
  TCanvas *can = new TCanvas("can","",1800,600);
  can->Divide(3,1);
  can->cd(1); gPad->SetLeftMargin(0.20); h1_dz->Draw("hist");
  can->cd(2); gPad->SetLeftMargin(0.20); h1_drphi->Draw("hist");
  can->cd(3); gPad->SetLeftMargin(0.20); h1_dr->Draw("hist");
  can->Update();
  can->SaveAs("figure/53534_residual_clusz.pdf");
}
