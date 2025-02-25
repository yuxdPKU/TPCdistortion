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

void plot()
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

  TH2* h2_drphi_r = new TH2F("h2_drphi_r","Rd#phi vs. R when |Z|<20;R (cm);Rd#phi (cm)",16,20,78,50,-2,2);

  int nevent = intree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    intree->GetEntry(i);
    double erp = square(clusRPhiErr) + square(stateRPhiErr);
    double ez = square(clusZErr) + square(stateZErr);
    //if (erp < 0.005) continue; if (ez < 0.01) continue; // bug
    if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
    if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
    if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;

    if (std::fabs(clusZ)<20) h2_drphi_r->Fill(clusR,drphi);
  }

  /*
  TCanvas *can = new TCanvas("can","",800,600);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  h2_drphi_r->Draw("colz");
  can->Update();
  can->SaveAs("figure/53534_drphi_r.pdf");
  */

  fit_gauss(h2_drphi_r, Form("figure/53534_drphi_r.pdf"),1);
}

void fit_gauss(TH2* h, TString name, bool verbose=0)
{
  TGraphErrors *graph = new TGraphErrors();

  int n = 1;
  for (int i = 1; i <= h->GetNbinsX(); i+=n) {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d R=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    //TF1 *gausFit = new TF1("gausFit", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-0.5;
    double fitmax=maxBinCenter+0.5;
    TF1 *gausFit = new TF1("gausFit", "gaus", fitmin, fitmax);

    projection->Fit(gausFit, "Q", "", fitmin, fitmax);

    if (verbose>0)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      projection->Draw("hist");
      gausFit->SetLineColor(kRed);
      gausFit->Draw("same");
      can->Update();
      can->SaveAs(name + Form(".gausFit_%d.pdf",i));
      delete can;
    }

    Double_t mean = gausFit->GetParameter(1);
    Double_t error = gausFit->GetParError(1);

    double center = 0;
    double width = 0;
    for (int j=0; j<n; j++)
    {
      center += h->GetXaxis()->GetBinCenter(i+j) / n;
      width += h->GetXaxis()->GetBinWidth(i+j) / 2;
    }
    graph->SetPoint(i-1, center, mean);
    graph->SetPointError(i-1, width, error);
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  h->Draw("COLZ");

  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kRed);
  graph->Draw("P");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  c1->Update();
  c1->SaveAs(name);
  delete c1;

}
