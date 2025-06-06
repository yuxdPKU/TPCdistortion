void fit_gauss(TTree* t, TH2* h, TString name, bool verbose=0);

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

void plot()
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
    float drphi, dz, clusZErr, stateZErr, clusZ, stateZ, clusR, stateR, clusRPhiErr, stateRPhiErr, clusPhi, statePhi;
    float tanAlpha, tanBeta;
    float trackPt;
    int layer, charge;
    int track_nmvtx, track_nintt, track_ntpc, track_ntpot;
    int track_nmvtxstate, track_ninttstate, track_ntpcstate, track_ntpotstate;
    intree->SetBranchAddress("dz",&dz);
    intree->SetBranchAddress("drphi",&drphi);
    intree->SetBranchAddress("clusZErr",&clusZErr);
    intree->SetBranchAddress("stateZErr",&stateZErr);
    intree->SetBranchAddress("clusZ",&clusZ);
    intree->SetBranchAddress("stateZ",&stateZ);
    intree->SetBranchAddress("clusR",&clusR);
    intree->SetBranchAddress("stateR",&stateR);
    intree->SetBranchAddress("clusPhi",&clusPhi);
    intree->SetBranchAddress("statePhi",&statePhi);
    intree->SetBranchAddress("clusRPhiErr",&clusRPhiErr);
    intree->SetBranchAddress("stateRPhiErr",&stateRPhiErr);
    intree->SetBranchAddress("tanAlpha",&tanAlpha);
    intree->SetBranchAddress("trackPt",&trackPt);
    intree->SetBranchAddress("tanBeta",&tanBeta);
    intree->SetBranchAddress("layer",&layer);
    intree->SetBranchAddress("charge",&charge);
    intree->SetBranchAddress("track_nmvtx",&track_nmvtx);
    intree->SetBranchAddress("track_nintt",&track_nintt);
    intree->SetBranchAddress("track_ntpc",&track_ntpc);
    intree->SetBranchAddress("track_ntpot",&track_ntpot);
    intree->SetBranchAddress("track_nmvtxstate",&track_nmvtxstate);
    intree->SetBranchAddress("track_ninttstate",&track_ninttstate);
    intree->SetBranchAddress("track_ntpcstate",&track_ntpcstate);
    intree->SetBranchAddress("track_ntpotstate",&track_ntpotstate);
  
    TH2* h2_drphi_r = new TH2F("h2_drphi_r","Rd#phi vs. R @ |Z|<20;R (cm);Rd#phi (cm)",16,20,78,100,-2,2);
  
    int layer_last=0;
    int layer_now=0;
    std::vector<double> vec_clusZ, vec_clusR, vec_drphi;
    vec_clusZ.clear();
    vec_clusR.clear();
    vec_drphi.clear();
    int nevent = intree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      intree->GetEntry(i);
      double erp = square(clusRPhiErr) + square(stateRPhiErr);
      double ez = square(clusZErr) + square(stateZErr);
      layer_now = layer;
      /*
      if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
      if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
      if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;
      if (track_nmvtxstate<3) continue;
      if (track_ninttstate<2) continue;
      if (track_ntpc<35) continue;
      if (track_ntpotstate<2) continue;
  
      if (std::fabs(clusZ)>=20) continue;
      */
      {
        //if (std::fabs(clusZ)>=25) continue;
        if (std::fabs(stateZ)>=10) continue;
        if (trackPt<0.5) continue;
	if (track_nmvtx<3) continue;
	if (track_nintt<2) continue;
	if (track_ntpot<2) continue;
	if (std::fabs(stateR * deltaPhi(clusPhi - statePhi))>2) continue;
      }
      if (layer_now>layer_last)
      {
	// still the same track
        vec_clusZ.push_back(clusZ);
        //vec_clusR.push_back(clusR);
        vec_clusR.push_back(stateR);
        //vec_drphi.push_back(drphi);
        vec_drphi.push_back(stateR * deltaPhi(clusPhi - statePhi));
      }
      else
      {
        // another track
	int nclus = vec_clusR.size();
	//if (vec_clusZ[0]*vec_clusZ[nclus-1]>0)
	{
	  for (int ii=0; ii<nclus; ii++)
	  {
	    h2_drphi_r->Fill(vec_clusR[ii],vec_drphi[ii]);
	  }
	}
	vec_clusZ.clear();
	vec_clusR.clear();
	vec_drphi.clear();
        vec_clusZ.push_back(clusZ);
        //vec_clusR.push_back(clusR);
        vec_clusR.push_back(stateR);
        //vec_drphi.push_back(drphi);
        vec_drphi.push_back(stateR * deltaPhi(clusPhi - statePhi));
      }
      layer_last = layer_now;
    }
  
    /*
    TCanvas *can = new TCanvas("can","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h2_drphi_r->Draw("colz");
    can->Update();
    can->SaveAs(Form("figure/%d_drphi_r.pdf",runs[k]));
    */
  
    fit_gauss(intree, h2_drphi_r, Form("figure/%d_drphi_r.pdf",runs[k]),1);
   
  }
}

void fit_gauss(TTree* t, TH2* h, TString name, bool verbose=0)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 0;
  for (int i = 1; i <= h->GetNbinsX(); i+=(n+1))
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d R=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-0.5;
    double fitmax=maxBinCenter+0.5;
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
      myText(0.2,0.95,kBlack,Form("%s , Bin %d R=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));
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
  gPad->SetLogz(1);
  c1->SetTopMargin(0.12);
  h->Draw("COLZ");
  {
  TProfile *profX = h->ProfileX("profX");
  profX->SetMarkerStyle(20);
  profX->SetMarkerColor(kViolet);
  profX->Draw("P,same");
  auto p = new TProfile( "p", "", 16, 20, 80 );
  //t->Draw("stateR * deltaPhi(clusPhi - statePhi):stateR>>p","stateZ<10 && stateZ>-10 && trackPt>0.5 && track_nmvtx>=3 && track_nintt>=2 && track_ntpot>=2","goff");
  //t->Draw("stateR*(clusPhi - statePhi):stateR>>p","stateZ<10 && stateZ>-10 && trackPt>0.5 && track_nmvtx>=3 && track_nintt>=2 && track_ntpot>=2");
  t->Draw("drphi:clusR>>p","stateZ<10 && stateZ>-10 && trackPt>0.5 && track_nmvtx>=3 && track_nintt>=2 && track_ntpot>=2");
  p->SetMarkerStyle(20);
  p->SetMarkerColor(kGreen);
  p->Draw("P,same");
  }


  graph_full->SetMarkerStyle(20);
  graph_full->SetMarkerColor(kRed);
  graph_full->Draw("P,same");
  graph_sub->SetMarkerStyle(20);
  graph_sub->SetMarkerColor(kBlack);
  graph_sub->Draw("P,same");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  myText(0.4,0.95,kBlack,Form("Rd#phi vs. R @ |Z|<20"));

  TFile* ofile = new TFile("hist_rdphi_vs_r.root","recreate");
  ofile->cd();
  graph_full->Write("graph_full");
  graph_sub->Write("graph_sub");
  h->Write();
  ofile->Close();

  c1->Update();
  c1->SaveAs(name);
  delete c1;

}
