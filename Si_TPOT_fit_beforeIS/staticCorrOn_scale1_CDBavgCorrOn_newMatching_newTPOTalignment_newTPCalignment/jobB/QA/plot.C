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

void plot()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  double twopi = 2*TMath::Pi();
  std::pair<double, double> phirange_central = {-1.73246, -1.43608};
  std::pair<double, double> phirange_west = {-1.21241, -0.909953};
  std::pair<double, double> phirange_east = {-2.26272, -1.96089};

  const int nrun = 1;
  int mbdrates[nrun] = {400};
  int runs[nrun] = {53877};

  for (int k=0; k<nrun; k++)
  {
    //TFile* infile = new TFile(Form("allqa_%d.root",runs[k]));
    //TTree* intree = (TTree*) infile->Get("h_QAG4SimulationDistortions_residTree");
    TChain* intree = new TChain("h_QAG4SimulationDistortions_residTree");
    for (int i=0; i<200; i++) {intree->Add(Form("../../Reconstructed/%d/clusters_seeds_%d-%d.root_qa.root",runs[k],runs[k],i));}
    //intree->Add(Form("../../Reconstructed/%d/*qa.root",runs[k]));
    float drphi, dz, clusZErr, stateZErr, clusZ, clusR, clusRPhiErr, stateRPhiErr, clusPhi;
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
    intree->SetBranchAddress("clusR",&clusR);
    intree->SetBranchAddress("clusPhi",&clusPhi);
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
  
    TH2* h2_drphi_r = new TH2F("h2_drphi_r","Rd#phi vs. R @ |Z|<20;R (cm);Rd#phi (cm)",16,20,78,50,-2,2);
    const int nbin_z = 6;
    double arr_z[nbin_z] = {5,10,15,30,60,80};
    TH2* arr_h2_drphi_r_posz[nbin_z];
    TH2* arr_h2_drphi_r_negz[nbin_z];
    TH2* arr_h2_dz_r_posz[nbin_z];
    TH2* arr_h2_dz_r_negz[nbin_z];
    for (int i=0; i<nbin_z; i++)
    {
      arr_h2_drphi_r_posz[i] = new TH2F(Form("h2_drphi_r_atposZ%d",(int)arr_z[i]),Form("Rd#phi vs. R @ Z=%d#pm2cm;R (cm);Rd#phi (cm)",(int)arr_z[i]),16,20,78,50,-2,2);
      arr_h2_drphi_r_negz[i] = new TH2F(Form("h2_drphi_r_atnegZ%d",(int)arr_z[i]),Form("Rd#phi vs. R @ Z=%d#pm2cm;R (cm);Rd#phi (cm)",-1*(int)arr_z[i]),16,20,78,50,-2,2);
      arr_h2_dz_r_posz[i] = new TH2F(Form("h2_dz_r_atposZ%d",(int)arr_z[i]),Form("dz vs. R @ Z=%d#pm2cm;R (cm);dz (cm)",(int)arr_z[i]),16,20,78,50,-1,1);
      arr_h2_dz_r_negz[i] = new TH2F(Form("h2_dz_r_atnegZ%d",(int)arr_z[i]),Form("dz vs. R @ Z=%d#pm2cm;R (cm);dz (cm)",-1*(int)arr_z[i]),16,20,78,50,-1,1);
    }
    const int nbin_r = 10;
    double arr_r[nbin_r] = {75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    TH2* arr_h2_dz_z_central[nbin_r];
    TH2* arr_h2_dz_z_west[nbin_r];
    TH2* arr_h2_dz_z_east[nbin_r];
    for (int i=0; i<nbin_r; i++)
    {
      arr_h2_dz_z_central[i] = new TH2F(Form("h2_dz_z_atR%d_central",(int)arr_r[i]),Form("Central tiles dz vs. Z @ R=%dcm;Z (cm);dz (cm)",(int)arr_r[i]),80,-105.5,105.5,50,-3,3);
      arr_h2_dz_z_west[i] = new TH2F(Form("h2_dz_z_atR%d_west",(int)arr_r[i]),Form("West tiles dz vs. Z @ R=%dcm;Z (cm);dz (cm)",(int)arr_r[i]),80,-105.5,105.5,50,-3,3);
      arr_h2_dz_z_east[i] = new TH2F(Form("h2_dz_z_atR%d_east",(int)arr_r[i]),Form("East tiles dz vs. Z @ R=%dcm;Z (cm);dz (cm)",(int)arr_r[i]),80,-105.5,105.5,50,-3,3);
    }
    int layer_last=0;
    int layer_now=0;
    std::vector<double> vec_clusZ, vec_clusR, vec_drphi;
    vec_clusZ.clear();
    vec_clusR.clear();
    vec_drphi.clear();
    int nevent = intree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
      intree->GetEntry(i);
      double erp = square(clusRPhiErr) + square(stateRPhiErr);
      double ez = square(clusZErr) + square(stateZErr);
      layer_now = layer;
      if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
      if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
      if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;
      if (track_nmvtxstate<3) continue;
      if (track_ninttstate<2) continue;
      if (track_ntpc<35) continue;
      if (track_ntpotstate<2) continue;

      for (int j=0; j<nbin_z; j++)
      {
        if (fabs(clusZ)>arr_z[j]-2 && fabs(clusZ)<arr_z[j]+2)
	{
	  if (clusZ>0)
	  {
	    arr_h2_drphi_r_posz[j]->Fill(clusR,drphi);
	    arr_h2_dz_r_posz[j]->Fill(clusR,dz);
	  }
	  else if (clusZ<0)
	  {
	    arr_h2_drphi_r_negz[j]->Fill(clusR,drphi);
	    arr_h2_dz_r_negz[j]->Fill(clusR,dz);
	  }
	}
      }

      double rbinwidth = (78 - 20) / 16.;
      for (int j=0; j<nbin_r; j++)
      {
        if (fabs(clusR)>arr_r[j]-rbinwidth/2 && fabs(clusR)<arr_r[j]+rbinwidth/2)
	{
          if (clusPhi>phirange_central.first && clusPhi<phirange_central.second)
	  {
	    arr_h2_dz_z_central[j]->Fill(clusZ,dz);
	  }
	  else if (clusPhi>phirange_west.first && clusPhi<phirange_west.second)
	  {
	    arr_h2_dz_z_west[j]->Fill(clusZ,dz);
	  }
	  else if (clusPhi>phirange_east.first && clusPhi<phirange_east.second)
	  {
	    arr_h2_dz_z_east[j]->Fill(clusZ,dz);
	  }
	}
      }

      if (std::fabs(clusZ)>=20) continue;
      if (layer_now>layer_last)
      {
	// still the same track
        vec_clusZ.push_back(clusZ);
        vec_clusR.push_back(clusR);
        vec_drphi.push_back(drphi);
      }
      else
      {
        // another track
	int nclus = vec_clusR.size();
	if (vec_clusZ[0]*vec_clusZ[nclus-1]>0)
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
        vec_clusR.push_back(clusR);
        vec_drphi.push_back(drphi);
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
  
    fit_gauss_Rslice(h2_drphi_r, Form("figure/%d_drphi_r.pdf",runs[k]),1);
    for (int i=0; i<nbin_z; i++)
    {
      fit_gauss_Rslice(arr_h2_drphi_r_posz[i], Form("figure/%d_drphi_r_posz%d.pdf",runs[k],(int)arr_z[i]),0);
      fit_gauss_Rslice(arr_h2_drphi_r_negz[i], Form("figure/%d_drphi_r_negz%d.pdf",runs[k],(int)arr_z[i]),0);
      fit_gauss_Rslice(arr_h2_dz_r_posz[i], Form("figure/%d_dz_r_posz%d.pdf",runs[k],(int)arr_z[i]),0);
      fit_gauss_Rslice(arr_h2_dz_r_negz[i], Form("figure/%d_dz_r_negz%d.pdf",runs[k],(int)arr_z[i]),0);
    }

    double slope_central[nbin_r], slope_west[nbin_r], slope_east[nbin_r];
    for (int i=0; i<nbin_r; i++)
    {
      slope_central[i] = fit_gauss_Zslice(arr_h2_dz_z_central[i], Form("figure/%d_dz_z_R%d_central.pdf",runs[k],(int)arr_r[i]),1);
      slope_west[i] = fit_gauss_Zslice(arr_h2_dz_z_west[i], Form("figure/%d_dz_z_R%d_west.pdf",runs[k],(int)arr_r[i]),1);
      slope_east[i] = fit_gauss_Zslice(arr_h2_dz_z_east[i], Form("figure/%d_dz_z_R%d_east.pdf",runs[k],(int)arr_r[i]),1);
    }
    TGraph* graph_slope_central = new TGraph(nbin_r,arr_r,slope_central);
    TGraph* graph_slope_west = new TGraph(nbin_r,arr_r,slope_west);
    TGraph* graph_slope_east = new TGraph(nbin_r,arr_r,slope_east);

    double ymin = min(min(graph_slope_central->GetHistogram()->GetMinimum(),graph_slope_west->GetHistogram()->GetMinimum()),graph_slope_east->GetHistogram()->GetMinimum());
    double ymax = max(max(graph_slope_central->GetHistogram()->GetMaximum(),graph_slope_west->GetHistogram()->GetMaximum()),graph_slope_east->GetHistogram()->GetMaximum());
    graph_slope_central->SetMinimum(ymin);
    graph_slope_west->SetMinimum(ymin);
    graph_slope_east->SetMinimum(ymin);
    graph_slope_central->SetMaximum(ymax);
    graph_slope_west->SetMaximum(ymax);
    graph_slope_east->SetMaximum(ymax);

    graph_slope_central->GetXaxis()->SetTitle("Raidus (cm)");
    graph_slope_central->GetYaxis()->SetTitle("Slope");

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetTopMargin(0.05);
    graph_slope_central->SetMarkerStyle(20);
    graph_slope_central->SetMarkerColor(kBlack);
    graph_slope_central->Draw("APL");
    graph_slope_west->SetMarkerStyle(20);
    graph_slope_west->SetMarkerColor(kBlue);
    graph_slope_west->Draw("PL,same");
    graph_slope_east->SetMarkerStyle(20);
    graph_slope_east->SetMarkerColor(kRed);
    graph_slope_east->Draw("PL,same");
    TLegend *legend = new TLegend(0.35, 0.65, 0.50, 0.9);
    legend->AddEntry(graph_slope_central, "Central tiles", "p");
    legend->AddEntry(graph_slope_west, "West tiles", "p");
    legend->AddEntry(graph_slope_east, "East tiles", "p");
    legend->SetTextSize(0.02);
    legend->Draw("same");
    c1->Update();
    c1->SaveAs("figure/slope_radius.pdf");

  }
}

void fit_gauss_Rslice(TH2* h, TString name, bool verbose=0)
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
  c1->SetTopMargin(0.12);
  gPad->SetLogz(1);
  h->Draw("COLZ");

  graph_full->SetMarkerStyle(20);
  graph_full->SetMarkerColor(kRed);
  graph_full->Draw("P,same");
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

double fit_gauss_Zslice(TH2* h, TString name, int verbose=0)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 0;
  double xmin_south = 100, xmax_south = -100;
  double xmin_north = 100, xmax_north = -100;
  for (int i = 1; i <= h->GetNbinsX(); i+=(n+1))
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d Z=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    if (projection->Integral()>10)
    {
      double xnow = h->GetXaxis()->GetBinCenter(i);
      if (xnow > 0)
      {
        //north side
        if (xmin_north > xnow)
	{
	  xmin_north = xnow;
	}
	if (xmax_north < xnow)
	{
	  xmax_north = xnow;
	}
      }
      else
      {
        //south side
        if (xmin_south > xnow)
	{
	  xmin_south = xnow;
	}
	if (xmax_south < xnow)
	{
	  xmax_south = xnow;
	}
      } 
    }

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-0.5;
    double fitmax=maxBinCenter+0.5;
    TF1 *gausFit_sub = new TF1("gausFit_sub", "gaus", fitmin, fitmax);
    projection->Fit(gausFit_sub, "Q", "", fitmin, fitmax);

    if (verbose>1)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      can->SetTopMargin(0.12);
      projection->Draw("hist");
      gausFit_full->SetLineColor(kRed);
      gausFit_full->Draw("same");
      gausFit_sub->SetLineColor(kBlue);
      gausFit_sub->Draw("same");
      myText(0.2,0.95,kBlack,Form("%s , Bin %d Z=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));
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

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 800);
  c1->SetTopMargin(0.20);
  gPad->SetLogz(1);
  h->Draw("COLZ");

  graph_full->SetMarkerStyle(20);
  graph_full->SetMarkerColor(kRed);
  graph_full->Draw("P,same");
  graph_sub->SetMarkerStyle(20);
  graph_sub->SetMarkerColor(kBlack);
  graph_sub->Draw("P,same");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  myText(0.2,0.95,kBlack,h->GetTitle());

  /*
  TFile* ofile = new TFile("hist_rdphi_vs_r.root","recreate");
  ofile->cd();
  graph_full->Write("graph_full");
  graph_sub->Write("graph_sub");
  h->Write();
  ofile->Close();
  */

  /*
  TF1 *linearFit_south = new TF1("linearFit_south", "[0]+[1]*x", xmin_south, xmax_south);
  linearFit_south->SetParameter(0,-0.01);
  linearFit_south->SetParameter(1,0);
  linearFit_south->SetParLimits(0,-1,1);
  linearFit_south->SetParLimits(1,-0.05,0.05);
  linearFit_south->SetLineColor(kViolet);
  graph_sub->Fit("linearFit_south", "QR");

  TF1 *linearFit_north = new TF1("linearFit_north", "[0]+[1]*x", xmin_north, xmax_north);
  linearFit_north->SetParameter(0,0.01);
  linearFit_north->SetParameter(1,5e-3);
  linearFit_north->SetParLimits(0,-1,1);
  linearFit_north->SetParLimits(1,-0.05,0.05);
  linearFit_north->SetLineColor(kViolet);
  graph_sub->Fit("linearFit_north", "QR");
  */

  TF1 *linearFit_south = new TF1("linearFit_south", "[0]+[1]*x", xmin_south, xmax_south);
  TF1 *linearFit_north = new TF1("linearFit_north", "[0]+[2]*x", xmin_north, xmax_north);

  linearFit_south->SetLineColor(kViolet);
  linearFit_north->SetLineColor(kViolet);

  linearFit_south->SetParameter(0, -0.01);
  linearFit_south->SetParameter(1, 0.01);
  linearFit_north->SetParameter(0, 0.01);
  linearFit_north->SetParameter(2, 0.01);
  
  linearFit_south->SetParLimits(0, -1, 1);
  linearFit_south->SetParLimits(1, -0.05, 0.05);
  linearFit_north->SetParLimits(0, -1, 1);
  linearFit_north->SetParLimits(2, -0.05, 0.05);

  linearFit_north->SetParameter(2, linearFit_south->GetParameter(1));

  graph_sub->Fit(linearFit_south, "QR");
  graph_sub->Fit(linearFit_north, "QR+");

  graph_sub->Fit(linearFit_south, "QR");
  graph_sub->Fit(linearFit_north, "QR F+");

  linearFit_south->Draw("same");
  linearFit_north->Draw("same");

  if (verbose>0)
  {
    cout<<"xmin_south = "<<xmin_south<<" , xmax_south = "<<xmax_south<<endl;
    cout << "south fit intercept = " << linearFit_south->GetParameter(0) << " +/- " << linearFit_south->GetParError(0) << endl;
    cout << "south fit slope = " << linearFit_south->GetParameter(1) << " +/- " << linearFit_south->GetParError(1) << endl;
    cout<<"xmin_north = "<<xmin_north<<" , xmax_north = "<<xmax_north<<endl;
    cout << "north fit intercept = " << linearFit_north->GetParameter(0) << " +/- " << linearFit_north->GetParError(0) << endl;
    cout << "north fit slope = " << linearFit_north->GetParameter(2) << " +/- " << linearFit_north->GetParError(2) << endl;
  }

  double intercept_south = linearFit_south->GetParameter(0);
  double intercept_south_err = linearFit_south->GetParError(0);
  double slope_south = linearFit_south->GetParameter(1);
  double slope_south_err = linearFit_south->GetParError(1);

  double intercept_north = linearFit_north->GetParameter(0);
  double intercept_north_err = linearFit_north->GetParError(0);
  double slope_north = linearFit_north->GetParameter(2);
  double slope_north_err = linearFit_north->GetParError(2);

  double chi2_south = linearFit_south->GetChisquare();
  int ndf_south = linearFit_south->GetNDF();
  double chi2_ndf_south = chi2_south / ndf_south;

  double chi2_north = linearFit_north->GetChisquare();
  int ndf_north = linearFit_north->GetNDF();
  double chi2_ndf_north = chi2_north / ndf_north;

  cout<<"chi2_ndf_south = "<<chi2_ndf_south<<" , chi2_ndf_north = "<<chi2_ndf_north<<endl;

  TPaveText *pt = new TPaveText(.05, .82, .95, .92, "NDC");
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText(Form("Slope: (%.2f#pm%.2f)#times10^{-2} (south), (%.2f#pm%.2f)#times10^{-2} (north)",100*slope_south,100*slope_south_err,100*slope_north,100*slope_north_err));
  pt->AddText(Form("Intercept: (%.2f#pm%.2f)#times10^{-2} (south), (%.2f#pm%.2f)#times10^{-2} (north)",100*intercept_south,100*intercept_south_err,100*intercept_north,100*intercept_north_err));
  pt->Draw("same");

  c1->Update();
  c1->SaveAs(name);
  delete c1;

  if (chi2_ndf_south < 50 && chi2_ndf_north < 50)
  {
    return (slope_south + slope_north) / 2;
  }
  else if (chi2_ndf_south < 50 && chi2_ndf_north >= 50)
  {
    return slope_south;
  }
  else if (chi2_ndf_south >= 50 && chi2_ndf_north < 50)
  {
    return slope_north;
  }
  else
  {
    return -100;
  }

}
