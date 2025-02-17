void drawCombine()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

//53877 - 400khz
//53876 - 430khz
//53756 - 380khz
//53744 - 300khz
//53630 - 550khz
//53534 - 250khz
//53285 - 70khz

int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};

TFile* file_1D_map_sim = new TFile("default_map.root","");
TH1* h_R_pos_sim = (TH1*) file_1D_map_sim->Get("h_dr_r_pos");
TH1* h_P_pos_sim = (TH1*) file_1D_map_sim->Get("h_dphi_r_pos");
TH1* h_Z_pos_sim = (TH1*) file_1D_map_sim->Get("h_dz_r_pos");
TH1* h_RP_pos_sim = (TH1*) file_1D_map_sim->Get("h_rdphi_r_pos");
TH1* h_R_neg_sim = (TH1*) file_1D_map_sim->Get("h_dr_r_neg");
TH1* h_P_neg_sim = (TH1*) file_1D_map_sim->Get("h_dphi_r_neg");
TH1* h_Z_neg_sim = (TH1*) file_1D_map_sim->Get("h_dz_r_neg");
TH1* h_RP_neg_sim = (TH1*) file_1D_map_sim->Get("h_rdphi_r_neg");

h_R_pos_sim->SetLineColor(kBlack);
h_P_pos_sim->SetLineColor(kBlack);
h_Z_pos_sim->SetLineColor(kBlack);
h_RP_pos_sim->SetLineColor(kBlack);
h_R_neg_sim->SetLineColor(kBlack);
h_P_neg_sim->SetLineColor(kBlack);
h_Z_neg_sim->SetLineColor(kBlack);
h_RP_neg_sim->SetLineColor(kBlack);

h_R_pos_sim->SetMinimum(-0.6); h_R_pos_sim->SetMaximum(1.4);
h_P_pos_sim->SetMinimum(-0.005); h_P_pos_sim->SetMaximum(0.05);
h_Z_pos_sim->SetMinimum(-0.4); h_Z_pos_sim->SetMaximum(1.3);
h_RP_pos_sim->SetMinimum(-0.3); h_RP_pos_sim->SetMaximum(1.7);
h_R_pos_sim->SetTitle("dR vs R Positive Z");
h_P_pos_sim->SetTitle("dphi vs R Positive Z");
h_Z_pos_sim->SetTitle("dz vs R Positive Z");
h_RP_pos_sim->SetTitle("Rdphi vs R Positive Z");

h_R_neg_sim->SetMinimum(-0.6); h_R_neg_sim->SetMaximum(1.4);
h_P_neg_sim->SetMinimum(-0.005); h_P_neg_sim->SetMaximum(0.05);
h_Z_neg_sim->SetMinimum(-0.4); h_Z_neg_sim->SetMaximum(1.3);
h_RP_neg_sim->SetMinimum(-0.3); h_RP_neg_sim->SetMaximum(1.7);
h_R_neg_sim->SetTitle("dR vs R Negative Z");
h_P_neg_sim->SetTitle("dphi vs R Negative Z");
h_Z_neg_sim->SetTitle("dz vs  Negative Z");
h_RP_neg_sim->SetTitle("rdphi vs R Negative Z");

TFile* file_1D_map[7];
TH1* h_N_pos[7];
TH1* h_R_pos[7];
TH1* h_P_pos[7];
TH1* h_Z_pos[7];
TH1* h_RP_pos[7];
TH1* h_N_neg[7];
TH1* h_R_neg[7];
TH1* h_P_neg[7];
TH1* h_Z_neg[7];
TH1* h_RP_neg[7];
for (int i = 1; i < 7; i++)
{
  file_1D_map[i] = new TFile(Form("Rootfiles/Distortions_1D_mm_%d_radius.root",runs[i]),"");
  h_R_pos[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionR_posz");
  h_P_pos[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionP_posz");
  h_Z_pos[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionZ_posz");
  h_RP_pos[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionRP_posz");
  h_N_pos[i] = (TH1*) file_1D_map[i]->Get("hentries_posz");
  h_R_neg[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionR_negz");
  h_P_neg[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionP_negz");
  h_Z_neg[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionZ_negz");
  h_RP_neg[i] = (TH1*) file_1D_map[i]->Get("hIntDistortionRP_negz");
  h_N_neg[i] = (TH1*) file_1D_map[i]->Get("hentries_negz");
  h_R_pos[i]->SetLineColor(i+2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(i+2);
  h_P_pos[i]->SetLineColor(i+2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(i+2);
  h_Z_pos[i]->SetLineColor(i+2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(i+2);
  h_RP_pos[i]->SetLineColor(i+2); h_RP_pos[i]->SetLineWidth(1); h_RP_pos[i]->SetFillColor(0); h_RP_pos[i]->SetMarkerColor(i+2);
  h_R_neg[i]->SetLineColor(i+2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_P_neg[i]->SetLineColor(i+2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(i+2);
  h_Z_neg[i]->SetLineColor(i+2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(i+2);
  h_RP_neg[i]->SetLineColor(i+2); h_RP_neg[i]->SetLineWidth(1); h_RP_neg[i]->SetFillColor(0); h_RP_neg[i]->SetMarkerColor(i+2);
}

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
h_P_pos_sim->Draw("hist"); for (int i=1; i<7; i++) h_P_pos[i]->Draw("hist,e,same"); h_P_pos_sim->Draw("hist,same");
can->cd(2);
h_RP_pos_sim->Draw("hist"); for (int i=1; i<7; i++) h_RP_pos[i]->Draw("hist,e,same"); h_RP_pos_sim->Draw("hist,same");
can->cd(3);
h_R_pos_sim->Draw("hist"); for (int i=1; i<7; i++) h_R_pos[i]->Draw("hist,e,same"); h_R_pos_sim->Draw("hist,same");
can->cd(4);
h_Z_pos_sim->Draw("hist"); for (int i=1; i<7; i++) h_Z_pos[i]->Draw("hist,e,same"); h_Z_pos_sim->Draw("hist,same");
can->cd(5);
h_P_neg_sim->Draw("hist"); for (int i=1; i<7; i++) h_P_neg[i]->Draw("hist,e,same"); h_P_neg_sim->Draw("hist,same");
can->cd(6);
h_RP_neg_sim->Draw("hist"); for (int i=1; i<7; i++) h_RP_neg[i]->Draw("hist,e,same"); h_RP_neg_sim->Draw("hist,same");
can->cd(7);
h_R_neg_sim->Draw("hist"); for (int i=1; i<7; i++) h_R_neg[i]->Draw("hist,e,same"); h_R_neg_sim->Draw("hist,same");
can->cd(8);
h_Z_neg_sim->Draw("hist"); for (int i=1; i<7; i++) h_Z_neg[i]->Draw("hist,e,same"); h_Z_neg_sim->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs("figure/resid_radius.pdf");

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->AddEntry(h_P_pos_sim, "Default simulation distortion map", "l");
for (int i=1; i<7; i++) legend->AddEntry(h_P_pos[i], Form("Run %d, MBD rate %d kHz",runs[i],mbdrates[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs("figure/resid_radius_leg.pdf");
}
