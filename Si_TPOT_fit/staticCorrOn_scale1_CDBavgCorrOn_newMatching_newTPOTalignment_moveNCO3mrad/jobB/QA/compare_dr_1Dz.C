void compare_dr_1Dz()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
{

const int nrun = 1;
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

TH1 *h1_R_rphi_posz[nrun], *h1_R_zr_posz[nrun], *h1_R_matrix_posz[nrun];
TH1 *h1_R_rphi_negz[nrun], *h1_R_zr_negz[nrun], *h1_R_matrix_negz[nrun];

TFile* infile_raw = new TFile(Form("hist_residual_1Dz_R%d.root",(int)selectR));
infile_raw->cd();
for (int i=0; i<nrun; i++)
{
  h1_R_rphi_posz[i] = (TH1*) infile_raw->Get(Form("h1_R_rphi_posz_%d",runs[i]));
  h1_R_zr_posz[i] = (TH1*) infile_raw->Get(Form("h1_R_zr_posz_%d",runs[i]));
  h1_R_rphi_negz[i] = (TH1*) infile_raw->Get(Form("h1_R_rphi_negz_%d",runs[i]));
  h1_R_zr_negz[i] = (TH1*) infile_raw->Get(Form("h1_R_zr_negz_%d",runs[i]));

  h1_R_rphi_posz[i]->SetMarkerColor(kRed+i); h1_R_rphi_posz[i]->SetLineColor(kRed+i);
  h1_R_zr_posz[i]->SetMarkerColor(kBlue+i); h1_R_zr_posz[i]->SetLineColor(kBlue+i);
  h1_R_rphi_negz[i]->SetMarkerColor(kRed+i); h1_R_rphi_negz[i]->SetLineColor(kRed+i);
  h1_R_zr_negz[i]->SetMarkerColor(kBlue+i); h1_R_zr_negz[i]->SetLineColor(kBlue+i);
}

TFile* infile_matrix = new TFile(Form("../Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR));
infile_matrix->cd();
for (int i=0; i<nrun; i++)
{
  h1_R_matrix_posz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h1_R_matrix_negz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  h1_R_matrix_posz[i]->SetMarkerColor(kBlack+i); h1_R_matrix_posz[i]->SetLineColor(kBlack+i);
  h1_R_matrix_negz[i]->SetMarkerColor(kBlack+i); h1_R_matrix_negz[i]->SetLineColor(kBlack+i);
}

TLegend *legend_posz = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_posz->SetHeader(Form("PositiveZ R = %d cm",(int)selectR));
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_R_rphi_posz[i], Form("Slope of linear fit to R-#phi plane"), "l");
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_R_zr_posz[i], Form("Slope of linear fit to Z-R plane"), "l");
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_R_matrix_posz[i], Form("Matrix inversion"), "l");
legend_posz->Draw();

TLegend *legend_negz = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_negz->SetHeader(Form("NegativeZ R = %d cm",(int)selectR));
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_R_rphi_negz[i], Form("Slope of linear fit to R-#phi plane"), "l");
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_R_zr_negz[i], Form("Slope of linear fit to Z-R plane"), "l");
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_R_matrix_negz[i], Form("Matrix inversion"), "l");
legend_negz->Draw();

TCanvas* can = new TCanvas("can","",1600,1200);
can->Divide(2,2);
can->cd(1);
for (int i=0; i<nrun; i++)
{
  h1_R_rphi_posz[i]->Draw("e,hist,same");
  h1_R_zr_posz[i]->Draw("e,hist,same");
  h1_R_matrix_posz[i]->Draw("e,hist,same");
}
can->cd(2);
legend_posz->Draw();
can->cd(3);
for (int i=0; i<nrun; i++)
{
  h1_R_rphi_negz[i]->Draw("e,hist,same");
  h1_R_zr_negz[i]->Draw("e,hist,same");
  h1_R_matrix_negz[i]->Draw("e,hist,same");
}
can->cd(4);
legend_negz->Draw();
can->SaveAs(Form("figure/dr_vsZ_compare_R%d.pdf",(int)selectR));

}
}
