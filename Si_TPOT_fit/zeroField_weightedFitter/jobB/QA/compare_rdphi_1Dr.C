void compare_rdphi_1Dr()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

const int nrun = 1;
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

TH1 *h1_RP_rphi_posz[nrun], *h1_RP_matrix_posz[nrun], *h1_RP_cdb_posz[nrun];
TH1 *h1_RP_rphi_negz[nrun], *h1_RP_matrix_negz[nrun], *h1_RP_cdb_negz[nrun];

TFile* infile_raw = new TFile(Form("hist_residual_1Dr_Z%d.root",(int)selectZ));
infile_raw->cd();
for (int i=0; i<nrun; i++)
{
  h1_RP_rphi_posz[i] = (TH1*) infile_raw->Get(Form("h1_RP_rphi_posz_%d",runs[i]));
  h1_RP_rphi_negz[i] = (TH1*) infile_raw->Get(Form("h1_RP_rphi_negz_%d",runs[i]));

  h1_RP_rphi_posz[i]->SetMarkerColor(kRed+i); h1_RP_rphi_posz[i]->SetLineColor(kRed+i);
  h1_RP_rphi_negz[i]->SetMarkerColor(kRed+i); h1_RP_rphi_negz[i]->SetLineColor(kRed+i);
}

TFile* infile_matrix = new TFile(Form("../Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ));
infile_matrix->cd();
for (int i=0; i<nrun; i++)
{
  h1_RP_matrix_posz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionRP_posz_%d",runs[i]));
  h1_RP_cdb_posz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionP_posz_%d_cdb",runs[i]));
  h1_RP_matrix_negz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionRP_negz_%d",runs[i]));
  h1_RP_cdb_negz[i] = (TH1*) infile_matrix->Get(Form("hIntDistortionP_negz_%d_cdb",runs[i]));

  h1_RP_matrix_posz[i]->SetMarkerColor(kBlack+i); h1_RP_matrix_posz[i]->SetLineColor(kBlack+i);
  h1_RP_cdb_posz[i]->SetMarkerColor(kGreen+i); h1_RP_cdb_posz[i]->SetLineColor(kGreen+i);
  h1_RP_matrix_negz[i]->SetMarkerColor(kBlack+i); h1_RP_matrix_negz[i]->SetLineColor(kBlack+i);
  h1_RP_cdb_negz[i]->SetMarkerColor(kGreen+i); h1_RP_cdb_negz[i]->SetLineColor(kGreen+i);
}

TLegend *legend_posz = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_posz->SetHeader(Form("Z = %d cm",(int)selectZ));
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_RP_rphi_posz[i], Form("Intercept of linear fit to R-#phi plane"), "l");
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_RP_matrix_posz[i], Form("Matrix inversion"), "l");
for (int i=0; i<nrun; i++) legend_posz->AddEntry(h1_RP_cdb_posz[i], Form("Lamination Fitting in CDB"), "l");
legend_posz->Draw();

TLegend *legend_negz = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_negz->SetHeader(Form("Z = %d cm",-1*(int)selectZ));
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_RP_rphi_negz[i], Form("Intercept of linear fit to R-#phi plane"), "l");
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_RP_matrix_negz[i], Form("Matrix inversion"), "l");
for (int i=0; i<nrun; i++) legend_negz->AddEntry(h1_RP_cdb_negz[i], Form("Lamination Fitting in CDB"), "l");
legend_negz->Draw();

TCanvas* can = new TCanvas("can","",1600,1200);
can->Divide(2,2);
can->cd(1);
for (int i=0; i<nrun; i++)
{
  h1_RP_rphi_posz[i]->Draw("e,hist,same");
  h1_RP_matrix_posz[i]->Draw("e,hist,same");
  h1_RP_cdb_posz[i]->Draw("hist,same");
}
can->cd(2);
legend_posz->Draw();
can->cd(3);
for (int i=0; i<nrun; i++)
{
  h1_RP_rphi_negz[i]->Draw("e,hist,same");
  h1_RP_matrix_negz[i]->Draw("e,hist,same");
  h1_RP_cdb_negz[i]->Draw("hist,same");
}
can->cd(4);
legend_negz->Draw();
can->SaveAs(Form("figure/rdphi_vsR_compare_Z%d.pdf",(int)selectZ));

}
}
