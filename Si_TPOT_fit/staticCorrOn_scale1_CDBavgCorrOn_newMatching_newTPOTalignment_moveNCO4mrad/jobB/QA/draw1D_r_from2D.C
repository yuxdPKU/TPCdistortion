void plot1D_Ybin(TH2* h2, TH1* h1, float y)
{
  int ybin = h2->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h2->GetNbinsX(); i++)
  {
      h1->SetBinContent(i, h2->GetBinContent(i, ybin));
      h1->SetBinError(i, h2->GetBinError(i, ybin));
  }
}

void draw1D_r_from2D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

// runNo - CAD-MBD-NS CAD-ZDC-NS
// 53877 - 400khz 4758.536231884057
// 53876 - 430khz 5082.326086956524
// 53756 - 380khz 4471.421428571427
// 53744 - 300khz 3581.862318840581
// 53630 - 550khz 6849.317241379308
// 53534 - 250khz 3013.5338345864657
// 53285 - 70khz 787.2765957446811

const int nrun = 1;
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

//get dr map from raw residual
TFile* rawdrfile[nrun];
TH2 *hraw_dr_rphi_posz[nrun], *hraw_rdphi_rphi_posz[nrun], *hraw_dr_zr_posz[nrun], *hraw_dz_zr_posz[nrun], *hraw_counts_posz[nrun];
TH2 *hraw_dr_rphi_negz[nrun], *hraw_rdphi_rphi_negz[nrun], *hraw_dr_zr_negz[nrun], *hraw_dz_zr_negz[nrun], *hraw_counts_negz[nrun];

for (int i = 0; i < nrun; i++)
{
  rawdrfile[i] = new TFile(Form("hist_residual_2Drz_%d.root",runs[i]),"");
  hraw_dr_rphi_posz[i] = (TH2*) rawdrfile[i]->Get("h2_dr_rphi_posz");
  hraw_rdphi_rphi_posz[i] = (TH2*) rawdrfile[i]->Get("h2_rdphi_rphi_posz");
  hraw_dr_zr_posz[i] = (TH2*) rawdrfile[i]->Get("h2_dr_zr_posz");
  hraw_dz_zr_posz[i] = (TH2*) rawdrfile[i]->Get("h2_dz_zr_posz");
  hraw_counts_posz[i] = (TH2*) rawdrfile[i]->Get("h2_counts_posz");
  hraw_dr_rphi_negz[i] = (TH2*) rawdrfile[i]->Get("h2_dr_rphi_negz");
  hraw_rdphi_rphi_negz[i] = (TH2*) rawdrfile[i]->Get("h2_rdphi_rphi_negz");
  hraw_dr_zr_negz[i] = (TH2*) rawdrfile[i]->Get("h2_dr_zr_negz");
  hraw_dz_zr_negz[i] = (TH2*) rawdrfile[i]->Get("h2_dz_zr_negz");
  hraw_counts_negz[i] = (TH2*) rawdrfile[i]->Get("h2_counts_negz");
}

TH1 *h1_N_posz[nrun], *h1_R_rphi_posz[nrun], *h1_R_zr_posz[nrun], *h1_RP_rphi_posz[nrun], *h1_Z_zr_posz[nrun];
TH1 *h1_N_negz[nrun], *h1_R_rphi_negz[nrun], *h1_R_zr_negz[nrun], *h1_RP_rphi_negz[nrun], *h1_Z_zr_negz[nrun];
for (int i = 0; i < nrun; i++)
{
  h1_N_posz[i] = new TH1F(Form("h1_N_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),hraw_counts_posz[i]->GetXaxis()->GetNbins(),hraw_counts_posz[i]->GetXaxis()->GetXmin(),hraw_counts_posz[i]->GetXaxis()->GetXmax());
  h1_R_rphi_posz[i] = new TH1F(Form("h1_R_rphi_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),hraw_dr_rphi_posz[i]->GetXaxis()->GetNbins(),hraw_dr_rphi_posz[i]->GetXaxis()->GetXmin(),hraw_dr_rphi_posz[i]->GetXaxis()->GetXmax());
  h1_R_zr_posz[i] = new TH1F(Form("h1_R_zr_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),hraw_dr_zr_posz[i]->GetXaxis()->GetNbins(),hraw_dr_zr_posz[i]->GetXaxis()->GetXmin(),hraw_dr_zr_posz[i]->GetXaxis()->GetXmax());
  h1_RP_rphi_posz[i] = new TH1F(Form("h1_RP_rphi_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),hraw_rdphi_rphi_posz[i]->GetXaxis()->GetNbins(),hraw_rdphi_rphi_posz[i]->GetXaxis()->GetXmin(),hraw_rdphi_rphi_posz[i]->GetXaxis()->GetXmax());
  h1_Z_zr_posz[i] = new TH1F(Form("h1_Z_zr_posz_%d",runs[i]),Form("dZ @ Z=%d cm;R (cm);dZ (cm)",(int)selectZ),hraw_dz_zr_posz[i]->GetXaxis()->GetNbins(),hraw_dz_zr_posz[i]->GetXaxis()->GetXmin(),hraw_dz_zr_posz[i]->GetXaxis()->GetXmax());
  plot1D_Ybin(hraw_counts_posz[i],h1_N_posz[i],selectZ);
  plot1D_Ybin(hraw_dr_rphi_posz[i],h1_R_rphi_posz[i],selectZ);
  plot1D_Ybin(hraw_dr_zr_posz[i],h1_R_zr_posz[i],selectZ);
  plot1D_Ybin(hraw_rdphi_rphi_posz[i],h1_RP_rphi_posz[i],selectZ);
  plot1D_Ybin(hraw_dz_zr_posz[i],h1_Z_zr_posz[i],selectZ);
  h1_N_posz[i]->SetLineColor(i+2); h1_N_posz[i]->SetLineWidth(1); h1_N_posz[i]->SetFillColor(0); h1_N_posz[i]->SetMarkerColor(i+2);
  h1_R_rphi_posz[i]->SetLineColor(kRed+i); h1_R_rphi_posz[i]->SetLineWidth(1); h1_R_rphi_posz[i]->SetFillColor(0); h1_R_rphi_posz[i]->SetMarkerColor(kRed+2);
  h1_R_zr_posz[i]->SetLineColor(kBlue+i); h1_R_zr_posz[i]->SetLineWidth(1); h1_R_zr_posz[i]->SetFillColor(0); h1_R_zr_posz[i]->SetMarkerColor(kBlue+i);
  h1_RP_rphi_posz[i]->SetLineColor(kRed+i); h1_RP_rphi_posz[i]->SetLineWidth(1); h1_RP_rphi_posz[i]->SetFillColor(0); h1_RP_rphi_posz[i]->SetMarkerColor(kRed+1);
  h1_Z_zr_posz[i]->SetLineColor(kBlue+i); h1_Z_zr_posz[i]->SetLineWidth(1); h1_Z_zr_posz[i]->SetFillColor(0); h1_Z_zr_posz[i]->SetMarkerColor(kBlue+i);
  h1_R_rphi_posz[i]->SetMinimum(-1); h1_R_rphi_posz[i]->SetMaximum(1);
  h1_RP_rphi_posz[i]->SetMinimum(-1); h1_RP_rphi_posz[i]->SetMaximum(1);
  h1_R_zr_posz[i]->SetMinimum(-1); h1_R_zr_posz[i]->SetMaximum(1);
  h1_Z_zr_posz[i]->SetMinimum(-1); h1_Z_zr_posz[i]->SetMaximum(1);

  h1_N_negz[i] = new TH1F(Form("h1_N_negz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",-1*(int)selectZ),hraw_counts_negz[i]->GetXaxis()->GetNbins(),hraw_counts_negz[i]->GetXaxis()->GetXmin(),hraw_counts_negz[i]->GetXaxis()->GetXmax());
  h1_R_rphi_negz[i] = new TH1F(Form("h1_R_rphi_negz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",-1*(int)selectZ),hraw_dr_rphi_negz[i]->GetXaxis()->GetNbins(),hraw_dr_rphi_negz[i]->GetXaxis()->GetXmin(),hraw_dr_rphi_negz[i]->GetXaxis()->GetXmax());
  h1_R_zr_negz[i] = new TH1F(Form("h1_R_zr_negz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",-1*(int)selectZ),hraw_dr_zr_negz[i]->GetXaxis()->GetNbins(),hraw_dr_zr_negz[i]->GetXaxis()->GetXmin(),hraw_dr_zr_negz[i]->GetXaxis()->GetXmax());
  h1_RP_rphi_negz[i] = new TH1F(Form("h1_RP_rphi_negz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",-1*(int)selectZ),hraw_rdphi_rphi_negz[i]->GetXaxis()->GetNbins(),hraw_rdphi_rphi_negz[i]->GetXaxis()->GetXmin(),hraw_rdphi_rphi_negz[i]->GetXaxis()->GetXmax());
  h1_Z_zr_negz[i] = new TH1F(Form("h1_Z_zr_negz_%d",runs[i]),Form("dZ @ Z=%d cm;R (cm);dZ (cm)",-1*(int)selectZ),hraw_dz_zr_negz[i]->GetXaxis()->GetNbins(),hraw_dz_zr_negz[i]->GetXaxis()->GetXmin(),hraw_dz_zr_negz[i]->GetXaxis()->GetXmax());
  plot1D_Ybin(hraw_counts_negz[i],h1_N_negz[i],-1*selectZ);
  plot1D_Ybin(hraw_dr_rphi_negz[i],h1_R_rphi_negz[i],-1*selectZ);
  plot1D_Ybin(hraw_dr_zr_negz[i],h1_R_zr_negz[i],-1*selectZ);
  plot1D_Ybin(hraw_rdphi_rphi_negz[i],h1_RP_rphi_negz[i],-1*selectZ);
  plot1D_Ybin(hraw_dz_zr_negz[i],h1_Z_zr_negz[i],-1*selectZ);
  h1_N_negz[i]->SetLineColor(i+2); h1_N_negz[i]->SetLineWidth(1); h1_N_negz[i]->SetFillColor(0); h1_N_negz[i]->SetMarkerColor(i+2);
  h1_R_rphi_negz[i]->SetLineColor(kRed+i); h1_R_rphi_negz[i]->SetLineWidth(1); h1_R_rphi_negz[i]->SetFillColor(0); h1_R_rphi_negz[i]->SetMarkerColor(kRed+i);
  h1_R_zr_negz[i]->SetLineColor(kBlue+i); h1_R_zr_negz[i]->SetLineWidth(1); h1_R_zr_negz[i]->SetFillColor(0); h1_R_zr_negz[i]->SetMarkerColor(kBlue+i);
  h1_RP_rphi_negz[i]->SetLineColor(kRed+i); h1_RP_rphi_negz[i]->SetLineWidth(1); h1_RP_rphi_negz[i]->SetFillColor(0); h1_RP_rphi_negz[i]->SetMarkerColor(kRed+i);
  h1_Z_zr_negz[i]->SetLineColor(kBlue+i); h1_Z_zr_negz[i]->SetLineWidth(1); h1_Z_zr_negz[i]->SetFillColor(0); h1_Z_zr_negz[i]->SetMarkerColor(kBlue+i);
  h1_R_rphi_negz[i]->SetMinimum(-1); h1_R_rphi_negz[i]->SetMaximum(1);
  h1_R_zr_negz[i]->SetMinimum(-1); h1_R_zr_negz[i]->SetMaximum(1);
  h1_RP_rphi_negz[i]->SetMinimum(-1); h1_RP_rphi_negz[i]->SetMaximum(1);
  h1_Z_zr_negz[i]->SetMinimum(-1); h1_Z_zr_negz[i]->SetMaximum(1);
}

TFile* ofile = new TFile(Form("hist_residual_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h1_R_rphi_posz[i]->Write();
  h1_RP_rphi_posz[i]->Write();
  h1_R_zr_posz[i]->Write();
  h1_Z_zr_posz[i]->Write();
  h1_N_posz[i]->Write();
  h1_R_rphi_negz[i]->Write();
  h1_RP_rphi_negz[i]->Write();
  h1_R_zr_negz[i]->Write();
  h1_Z_zr_negz[i]->Write();
  h1_N_negz[i]->Write();
}
ofile->Write();

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
for (int i=0; i<nrun; i++) h1_R_rphi_posz[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h1_R_zr_posz[i]->Draw("hist,e,same");
can->cd(2);
for (int i=0; i<nrun; i++) h1_RP_rphi_posz[i]->Draw("hist,e,same");
can->cd(3);
for (int i=0; i<nrun; i++) h1_Z_zr_posz[i]->Draw("hist,e,same");
can->cd(4);
for (int i=0; i<nrun; i++) h1_N_posz[i]->Draw("hist,e,same");
can->cd(5);
for (int i=0; i<nrun; i++) h1_R_rphi_negz[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h1_R_zr_negz[i]->Draw("hist,e,same");
can->cd(6);
for (int i=0; i<nrun; i++) h1_RP_rphi_negz[i]->Draw("hist,e,same");
can->cd(7);
for (int i=0; i<nrun; i++) h1_Z_zr_negz[i]->Draw("hist,e,same");
can->cd(8);
for (int i=0; i<nrun; i++) h1_N_negz[i]->Draw("hist,e,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/raw_residual_vsR_from2D_atZ%d.pdf",(int)selectZ));

delete can;
}

}
