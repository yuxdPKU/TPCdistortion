#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void plot1D(TH2* h2, TH1* h1, double phi0 = 4.55073)
{
  int xbin = h2->GetXaxis()->FindBin(phi0);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      float bincontent = h2->GetBinContent(xbin, i);
      float binerror = h2->GetBinError(xbin, i);
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void draw1D_r()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

const int nrun = 1;
int runs[nrun] = {79516};

double phi0 = 4.7;

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  cdbfilename[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  cout<<"CDB file path: "<<cdbfilename[i].c_str()<<endl;
  cdbfile[i] = new TFile(cdbfilename[i].c_str(),"");
  hcdb_R_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_posz");
  hcdb_P_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_posz");
  hcdb_Z_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_posz");
  hcdb_N_pr_pos[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");
  hcdb_R_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_negz");
  hcdb_P_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_negz");
  hcdb_Z_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_negz");
  hcdb_N_pr_neg[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");

  hcdb_N_pos[i] = new TH1F(Form("hentries_posz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (rad)"),hcdb_P_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_N_neg[i] = new TH1F(Form("hentries_negz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (rad)"),hcdb_P_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmax());
  plot1D(hcdb_N_pr_pos[i],hcdb_N_pos[i],phi0);
  plot1D(hcdb_R_pr_pos[i],hcdb_R_pos[i],phi0);
  plot1D(hcdb_P_pr_pos[i],hcdb_P_pos[i],phi0);
  plot1D(hcdb_Z_pr_pos[i],hcdb_Z_pos[i],phi0);
  plot1D(hcdb_N_pr_neg[i],hcdb_N_neg[i],phi0);
  plot1D(hcdb_R_pr_neg[i],hcdb_R_neg[i],phi0);
  plot1D(hcdb_P_pr_neg[i],hcdb_P_neg[i],phi0);
  plot1D(hcdb_Z_pr_neg[i],hcdb_Z_neg[i],phi0);
  hcdb_N_pos[i]->SetLineColor(1); hcdb_N_pos[i]->SetLineWidth(1); hcdb_N_pos[i]->SetFillColor(0); hcdb_N_pos[i]->SetMarkerColor(1);
  hcdb_R_pos[i]->SetLineColor(1); hcdb_R_pos[i]->SetLineWidth(1); hcdb_R_pos[i]->SetFillColor(0); hcdb_R_pos[i]->SetMarkerColor(1);
  hcdb_P_pos[i]->SetLineColor(1); hcdb_P_pos[i]->SetLineWidth(1); hcdb_P_pos[i]->SetFillColor(0); hcdb_P_pos[i]->SetMarkerColor(1);
  hcdb_Z_pos[i]->SetLineColor(1); hcdb_Z_pos[i]->SetLineWidth(1); hcdb_Z_pos[i]->SetFillColor(0); hcdb_Z_pos[i]->SetMarkerColor(1);
  hcdb_N_neg[i]->SetLineColor(1); hcdb_N_neg[i]->SetLineWidth(1); hcdb_N_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_R_neg[i]->SetLineColor(1); hcdb_R_neg[i]->SetLineWidth(1); hcdb_R_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_P_neg[i]->SetLineColor(1); hcdb_P_neg[i]->SetLineWidth(1); hcdb_P_neg[i]->SetFillColor(0); hcdb_P_neg[i]->SetMarkerColor(1);
  hcdb_Z_neg[i]->SetLineColor(1); hcdb_Z_neg[i]->SetLineWidth(1); hcdb_Z_neg[i]->SetFillColor(0); hcdb_Z_neg[i]->SetMarkerColor(1);
}

//get map from my lamination study
std::string myfilename[nrun];
TFile* myfile[nrun];
TH2 *hmy_N_pr_pos[nrun], *hmy_R_pr_pos[nrun], *hmy_P_pr_pos[nrun], *hmy_Z_pr_pos[nrun];
TH2 *hmy_N_pr_neg[nrun], *hmy_R_pr_neg[nrun], *hmy_P_pr_neg[nrun], *hmy_Z_pr_neg[nrun];
TH1 *hmy_N_pos[nrun], *hmy_R_pos[nrun], *hmy_P_pos[nrun], *hmy_Z_pos[nrun];
TH1 *hmy_N_neg[nrun], *hmy_R_neg[nrun], *hmy_P_neg[nrun], *hmy_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  //myfile[i] = new TFile(Form("/sphenix/u/xyu3/hftg01/Lamination_DST/LaminationFitOut_20260131/rename/Laminations_run3pp_ana532_2025p009_v001-000%d.root",runs[i]),"");
  myfile[i] = new TFile(Form("/sphenix/u/xyu3/hftg01/Lamination_DST/LaminationFitOut_20260131_round2/rename/Laminations_run3pp_ana532_2025p009_v001-000%d.root",runs[i]),"");
  hmy_R_pr_pos[i] = (TH2*) myfile[i]->Get("hIntDistortionR_posz");
  hmy_P_pr_pos[i] = (TH2*) myfile[i]->Get("hIntDistortionP_posz");
  hmy_Z_pr_pos[i] = (TH2*) myfile[i]->Get("hIntDistortionZ_posz");
  hmy_N_pr_pos[i] = (TH2*) myfile[i]->Get("hEntries_negz");
  hmy_R_pr_neg[i] = (TH2*) myfile[i]->Get("hIntDistortionR_negz");
  hmy_P_pr_neg[i] = (TH2*) myfile[i]->Get("hIntDistortionP_negz");
  hmy_Z_pr_neg[i] = (TH2*) myfile[i]->Get("hIntDistortionZ_negz");
  hmy_N_pr_neg[i] = (TH2*) myfile[i]->Get("hEntries_negz");

  hmy_N_pos[i] = new TH1F(Form("hentries_posz_%d_my",runs[i]),Form("N;R (cm);N"),hmy_N_pr_pos[i]->GetYaxis()->GetNbins(),hmy_N_pr_pos[i]->GetYaxis()->GetXmin(),hmy_N_pr_pos[i]->GetYaxis()->GetXmax());
  hmy_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_my",runs[i]),Form("dR;R (cm);dR (cm)"),hmy_R_pr_pos[i]->GetYaxis()->GetNbins(),hmy_R_pr_pos[i]->GetYaxis()->GetXmin(),hmy_R_pr_pos[i]->GetYaxis()->GetXmax());
  hmy_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_my",runs[i]),Form("Rdphi;R (cm);Rdphi (rad)"),hmy_P_pr_pos[i]->GetYaxis()->GetNbins(),hmy_P_pr_pos[i]->GetYaxis()->GetXmin(),hmy_P_pr_pos[i]->GetYaxis()->GetXmax());
  hmy_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d_my",runs[i]),Form("dz;R (cm);dz (cm)"),hmy_Z_pr_pos[i]->GetYaxis()->GetNbins(),hmy_Z_pr_pos[i]->GetYaxis()->GetXmin(),hmy_Z_pr_pos[i]->GetYaxis()->GetXmax());
  hmy_N_neg[i] = new TH1F(Form("hentries_negz_%d_my",runs[i]),Form("N;R (cm);N"),hmy_N_pr_neg[i]->GetYaxis()->GetNbins(),hmy_N_pr_neg[i]->GetYaxis()->GetXmin(),hmy_N_pr_neg[i]->GetYaxis()->GetXmax());
  hmy_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_my",runs[i]),Form("dR;R (cm);dR (cm)"),hmy_R_pr_neg[i]->GetYaxis()->GetNbins(),hmy_R_pr_neg[i]->GetYaxis()->GetXmin(),hmy_R_pr_neg[i]->GetYaxis()->GetXmax());
  hmy_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_my",runs[i]),Form("Rdphi;R (cm);Rdphi (rad)"),hmy_P_pr_neg[i]->GetYaxis()->GetNbins(),hmy_P_pr_neg[i]->GetYaxis()->GetXmin(),hmy_P_pr_neg[i]->GetYaxis()->GetXmax());
  hmy_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d_my",runs[i]),Form("dz;R (cm);dz (cm)"),hmy_Z_pr_neg[i]->GetYaxis()->GetNbins(),hmy_Z_pr_neg[i]->GetYaxis()->GetXmin(),hmy_Z_pr_neg[i]->GetYaxis()->GetXmax());
  plot1D(hmy_N_pr_pos[i],hmy_N_pos[i],phi0);
  plot1D(hmy_R_pr_pos[i],hmy_R_pos[i],phi0);
  plot1D(hmy_P_pr_pos[i],hmy_P_pos[i],phi0);
  plot1D(hmy_Z_pr_pos[i],hmy_Z_pos[i],phi0);
  plot1D(hmy_N_pr_neg[i],hmy_N_neg[i],phi0);
  plot1D(hmy_R_pr_neg[i],hmy_R_neg[i],phi0);
  plot1D(hmy_P_pr_neg[i],hmy_P_neg[i],phi0);
  plot1D(hmy_Z_pr_neg[i],hmy_Z_neg[i],phi0);
  hmy_N_pos[i]->SetLineColor(2); hmy_N_pos[i]->SetLineWidth(1); hmy_N_pos[i]->SetFillColor(0); hmy_N_pos[i]->SetMarkerColor(2);
  hmy_R_pos[i]->SetLineColor(2); hmy_R_pos[i]->SetLineWidth(1); hmy_R_pos[i]->SetFillColor(0); hmy_R_pos[i]->SetMarkerColor(2);
  hmy_P_pos[i]->SetLineColor(2); hmy_P_pos[i]->SetLineWidth(1); hmy_P_pos[i]->SetFillColor(0); hmy_P_pos[i]->SetMarkerColor(2);
  hmy_Z_pos[i]->SetLineColor(2); hmy_Z_pos[i]->SetLineWidth(1); hmy_Z_pos[i]->SetFillColor(0); hmy_Z_pos[i]->SetMarkerColor(2);
  hmy_N_neg[i]->SetLineColor(2); hmy_N_neg[i]->SetLineWidth(1); hmy_N_neg[i]->SetFillColor(0); hmy_R_neg[i]->SetMarkerColor(2);
  hmy_R_neg[i]->SetLineColor(2); hmy_R_neg[i]->SetLineWidth(1); hmy_R_neg[i]->SetFillColor(0); hmy_R_neg[i]->SetMarkerColor(2);
  hmy_P_neg[i]->SetLineColor(2); hmy_P_neg[i]->SetLineWidth(1); hmy_P_neg[i]->SetFillColor(0); hmy_P_neg[i]->SetMarkerColor(2);
  hmy_Z_neg[i]->SetLineColor(2); hmy_Z_neg[i]->SetLineWidth(1); hmy_Z_neg[i]->SetFillColor(0); hmy_Z_neg[i]->SetMarkerColor(2);
}

for (int i=0; i<nrun; i++)
{
  TCanvas* can = new TCanvas("can","",1600,1200);
  can->Divide(2,2);
  can->cd(1);
  gPad->SetLogy(0);
  SetCommonYRange({hmy_P_pos[i],hcdb_P_pos[i]});
  hmy_P_pos[i]->Draw("hist");
  hcdb_P_pos[i]->Draw("hist,same");
  can->cd(2);
  gPad->SetLogy(0);
  SetCommonYRange({hmy_R_pos[i],hcdb_R_pos[i]});
  hmy_R_pos[i]->Draw("hist");
  hcdb_R_pos[i]->Draw("hist,same");
  can->cd(3);
  gPad->SetLogy(0);
  SetCommonYRange({hmy_P_neg[i],hcdb_P_neg[i]});
  hmy_P_neg[i]->Draw("hist");
  hcdb_P_neg[i]->Draw("hist,same");
  can->cd(4);
  gPad->SetLogy(0);
  hmy_R_neg[i]->Draw("hist");
  SetCommonYRange({hmy_R_neg[i],hcdb_R_neg[i]});
  hcdb_R_neg[i]->Draw("hist,same");

  gPad->RedrawAxis();

  can->Update();
  can->SaveAs(Form("figure/resid_vsR_%d.pdf",runs[i]));

  TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
  legend->SetHeader(Form("Run %d, phi = %.3f", runs[i], phi0));
  legend->AddEntry(hmy_P_pos[i], Form("My Lamination map"), "l");
  legend->AddEntry(hcdb_P_pos[i], Form("CDB Lamination map"), "l");
  legend->Draw();
  TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
  legend->Draw();
  can_leg->SaveAs(Form("figure/resid_vsR_%d_legend.pdf",runs[i]));

  delete can;
  delete can_leg;
}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin > 0) yMin*=0.9; else yMin*=1.1;
    if (yMax < 0) yMax*=0.9; else yMax*=1.1;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
