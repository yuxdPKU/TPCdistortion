#include <TGaxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPad.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

struct threeDHist {
    TH3* h_N_pos;
    TH3* h_P_pos;
    TH3* h_R_pos;
    TH3* h_Z_pos;
    TH3* h_N_neg;
    TH3* h_P_neg;
    TH3* h_R_neg;
    TH3* h_Z_neg;
};

struct twoDHist {
    TH2* h_N_pos;
    TH2* h_P_pos;
    TH2* h_R_pos;
    TH2* h_Z_pos;
    TH2* h_N_neg;
    TH2* h_P_neg;
    TH2* h_R_neg;
    TH2* h_Z_neg;
};

TH2* MakeRPhiHist(TH3* h3, TString name, TString title)
{
  return new TH2F(name, title,
                  h3->GetXaxis()->GetNbins(), h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax(),
                  h3->GetZaxis()->GetNbins(), h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax());
}

// 3D x: phi, y: r, z: z
// 2D x: phi, y: z
void plot2D_Ybin(TH3* h3, TH2* h2, float r)
{
  int rbin = h3->GetYaxis()->FindBin(r);
  for (int i = 1; i <= h3->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h3->GetNbinsZ(); j++)
    {
      double value = h3->GetBinContent(i, rbin, j);
      double error = h3->GetBinError(i, rbin, j);
      if (std::isnan(value))
      {
        value = 0.0;
        error = 0.0;
      }
      h2->SetBinContent(i, j, value);
      h2->SetBinError(i, j, error);
    }
  }
}

int Mask2DWithN(TH2* h_N, TH2* h)
{
  int nMasked = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h->GetNbinsY(); j++)
    {
      if (h_N->GetBinContent(i, j) <= 0)
      {
        h->SetBinContent(i, j, std::numeric_limits<double>::quiet_NaN());
        h->SetBinError(i, j, 0.0);
        nMasked++;
      }
    }
  }
  return nMasked;
}

void SetCommonZRange(const std::vector<TH2*>& histograms, bool positiveOnly=false)
{
  double zmin = 1e30;
  double zmax = -1e30;
  for (TH2* h : histograms)
  {
    if (!h) continue;
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
      for (int j = 1; j <= h->GetNbinsY(); j++)
      {
        double value = h->GetBinContent(i, j);
        if (!std::isfinite(value)) continue;
        zmin = std::min(zmin, value);
        zmax = std::max(zmax, value);
      }
    }
  }

  if (zmin == 1e30) zmin = 0;
  if (zmax == -1e30) zmax = 1;

  if (positiveOnly)
  {
    zmin = 0;
    if (zmax <= 0) zmax = 1;
  }
  else
  {
    if (zmin < 0) zmin *= 1.1; else zmin = -0.1;
    if (zmax > 0) zmax *= 1.1; else zmax = 0.1;
  }

  for (TH2* h : histograms)
  {
    if (!h) continue;
    h->SetMinimum(zmin);
    h->SetMaximum(zmax);
  }
}

void DrawMaskedBinsWhite(TH2* h_N)
{
  for (int i = 1; i <= h_N->GetNbinsX(); i++)
  {
    double xlow = h_N->GetXaxis()->GetBinLowEdge(i);
    double xup = h_N->GetXaxis()->GetBinUpEdge(i);
    for (int j = 1; j <= h_N->GetNbinsY(); j++)
    {
      if (h_N->GetBinContent(i, j) > 0) continue;
      double ylow = h_N->GetYaxis()->GetBinLowEdge(j);
      double yup = h_N->GetYaxis()->GetBinUpEdge(j);
      TBox* box = new TBox(xlow, ylow, xup, yup);
      box->SetFillColor(kWhite);
      box->SetLineColor(kWhite);
      box->Draw("same");
    }
  }
}

void Draw2D(TH2* h, TH2* h_N=nullptr)
{
  gPad->SetRightMargin(0.15);
  h->Draw(h_N ? "colz0" : "colz");
  if (h_N) DrawMaskedBinsWhite(h_N);
  gPad->RedrawAxis();
}

bool Load3DHists(TFile* file_3D_map, threeDHist& hists3D)
{
  hists3D.h_N_pos = (TH3*) file_3D_map->Get("hentries_posz");
  hists3D.h_P_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  hists3D.h_R_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  hists3D.h_Z_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  hists3D.h_N_neg = (TH3*) file_3D_map->Get("hentries_negz");
  hists3D.h_P_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  hists3D.h_R_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  hists3D.h_Z_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");

  return hists3D.h_N_pos && hists3D.h_P_pos && hists3D.h_R_pos && hists3D.h_Z_pos &&
         hists3D.h_N_neg && hists3D.h_P_neg && hists3D.h_R_neg && hists3D.h_Z_neg;
}

void draw2D_phiz_from3D_sim_closure(int run=29, TString tag="MI_sim_reco_truth_notpot", int selectR=0)
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);
gSystem->mkdir("figure", true);
gSystem->mkdir("Rootfiles", true);

TString filename = Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, tag.Data());
TFile* file_3D_map = new TFile(filename, "");
if (!file_3D_map || file_3D_map->IsZombie())
{
  std::cerr << "Error opening file: " << filename << std::endl;
  return;
}

threeDHist hists3D;
if (!Load3DHists(file_3D_map, hists3D))
{
  std::cerr << "Error: missing N/P/R/Z 3D histograms in " << filename << std::endl;
  delete file_3D_map;
  return;
}

int rbin = hists3D.h_N_pos->GetYaxis()->FindBin(selectR);
double rBinCenter = hists3D.h_N_pos->GetYaxis()->GetBinCenter(rbin);
std::cout << ", selectR = " << selectR
          << ", rbin = " << rbin
          << ", bin center = " << rBinCenter
          << std::endl;

twoDHist hists2D;
hists2D.h_N_pos = MakeRPhiHist(hists3D.h_N_pos, Form("hentries_posz_%s",tag.Data()), Form("N posz @ r=%d;#phi (rad);Z (cm);N",selectR));
hists2D.h_P_pos = MakeRPhiHist(hists3D.h_P_pos, Form("hIntDistortionP_posz_%s",tag.Data()), Form("d#phi posz @ r=%d;#phi (rad);Z (cm);d#phi (rad)",selectR));
hists2D.h_R_pos = MakeRPhiHist(hists3D.h_R_pos, Form("hIntDistortionR_posz_%s",tag.Data()), Form("dR posz @ r=%d;#phi (rad);Z (cm);dR (cm)",selectR));
hists2D.h_Z_pos = MakeRPhiHist(hists3D.h_Z_pos, Form("hIntDistortionZ_posz_%s",tag.Data()), Form("dz posz @ r=%d;#phi (rad);Z (cm);dz (cm)",selectR));
hists2D.h_N_neg = MakeRPhiHist(hists3D.h_N_neg, Form("hentries_negz_%s",tag.Data()), Form("N negz @ r=%d;#phi (rad);Z (cm);N",selectR));
hists2D.h_P_neg = MakeRPhiHist(hists3D.h_P_neg, Form("hIntDistortionP_negz_%s",tag.Data()), Form("d#phi negz @ r=%d;#phi (rad);Z (cm);d#phi (rad)",selectR));
hists2D.h_R_neg = MakeRPhiHist(hists3D.h_R_neg, Form("hIntDistortionR_negz_%s",tag.Data()), Form("dR negz @ r=%d;#phi (rad);Z (cm);dR (cm)",selectR));
hists2D.h_Z_neg = MakeRPhiHist(hists3D.h_Z_neg, Form("hIntDistortionZ_negz_%s",tag.Data()), Form("dz negz @ r=%d;#phi (rad);Z (cm);dz (cm)",selectR));

hists2D.h_N_pos->SetDirectory(0);
hists2D.h_P_pos->SetDirectory(0);
hists2D.h_R_pos->SetDirectory(0);
hists2D.h_Z_pos->SetDirectory(0);
hists2D.h_N_neg->SetDirectory(0);
hists2D.h_P_neg->SetDirectory(0);
hists2D.h_R_neg->SetDirectory(0);
hists2D.h_Z_neg->SetDirectory(0);

plot2D_Ybin(hists3D.h_N_pos, hists2D.h_N_pos, selectR);
plot2D_Ybin(hists3D.h_P_pos, hists2D.h_P_pos, selectR);
plot2D_Ybin(hists3D.h_R_pos, hists2D.h_R_pos, selectR);
plot2D_Ybin(hists3D.h_Z_pos, hists2D.h_Z_pos, selectR);
plot2D_Ybin(hists3D.h_N_neg, hists2D.h_N_neg, selectR);
plot2D_Ybin(hists3D.h_P_neg, hists2D.h_P_neg, selectR);
plot2D_Ybin(hists3D.h_R_neg, hists2D.h_R_neg, selectR);
plot2D_Ybin(hists3D.h_Z_neg, hists2D.h_Z_neg, selectR);

int nMaskedPos = Mask2DWithN(hists2D.h_N_pos, hists2D.h_P_pos);
Mask2DWithN(hists2D.h_N_pos, hists2D.h_R_pos);
Mask2DWithN(hists2D.h_N_pos, hists2D.h_Z_pos);
int nMaskedNeg = Mask2DWithN(hists2D.h_N_neg, hists2D.h_P_neg);
Mask2DWithN(hists2D.h_N_neg, hists2D.h_R_neg);
Mask2DWithN(hists2D.h_N_neg, hists2D.h_Z_neg);
std::cout << "Mask N=0 bins in P/R/Z: posz = " << nMaskedPos
          << ", negz = " << nMaskedNeg
          << std::endl;

SetCommonZRange({hists2D.h_N_pos, hists2D.h_N_neg}, true);
SetCommonZRange({hists2D.h_P_pos, hists2D.h_P_neg});
SetCommonZRange({hists2D.h_R_pos, hists2D.h_R_neg});
SetCommonZRange({hists2D.h_Z_pos, hists2D.h_Z_neg});

TCanvas* can = new TCanvas("can","",6400,2600);
can->Divide(4,2);
can->cd(1); Draw2D(hists2D.h_N_pos);
can->cd(2); Draw2D(hists2D.h_P_pos, hists2D.h_N_pos);
can->cd(3); Draw2D(hists2D.h_R_pos, hists2D.h_N_pos);
can->cd(4); Draw2D(hists2D.h_Z_pos, hists2D.h_N_pos);
can->cd(5); Draw2D(hists2D.h_N_neg);
can->cd(6); Draw2D(hists2D.h_P_neg, hists2D.h_N_neg);
can->cd(7); Draw2D(hists2D.h_R_neg, hists2D.h_N_neg);
can->cd(8); Draw2D(hists2D.h_Z_neg, hists2D.h_N_neg);
gPad->RedrawAxis();
can->Update();
can->SaveAs(Form("figure/phiz_NPRZ_from3D_%d_%s_R%d.pdf", run, tag.Data(), (int)selectR));
delete can;

TFile* outfile = new TFile(Form("Rootfiles/hist_NPRZ_2Dphiz_%d_%s_R%d.root", run, tag.Data(), (int)selectR), "recreate");
outfile->cd();
hists2D.h_N_pos->Write();
hists2D.h_P_pos->Write();
hists2D.h_R_pos->Write();
hists2D.h_Z_pos->Write();
hists2D.h_N_neg->Write();
hists2D.h_P_neg->Write();
hists2D.h_R_neg->Write();
hists2D.h_Z_neg->Write();
outfile->Write();
outfile->Close();

delete file_3D_map;
}
