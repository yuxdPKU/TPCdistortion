void plot1D_Ybin(TH2* h2, TH1* h1, float y)
{
  int ybin = h2->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h2->GetNbinsX(); i++)
  {
      h1->SetBinContent(i, h2->GetBinContent(i, ybin));
      h1->SetBinError(i, h2->GetBinError(i, ybin));
  }
}

void draw1D_from2D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectZs={5,10,15,30,60,80};
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

const int nrun = 6;
//int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
//int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};

TFile* file_1D_map_sim = new TFile("default_map.root","");
TH1* h_R_pos_sim = (TH1*) file_1D_map_sim->Get("h_dr_r_pos");
TH1* h_P_pos_sim = (TH1*) file_1D_map_sim->Get("h_dphi_r_pos");
TH1* h_Z_pos_sim = (TH1*) file_1D_map_sim->Get("h_dz_r_pos");
TH1* h_RP_pos_sim = (TH1*) file_1D_map_sim->Get("h_rdphi_r_pos");
TH1* h_R_neg_sim = (TH1*) file_1D_map_sim->Get("h_dr_r_neg");
TH1* h_P_neg_sim = (TH1*) file_1D_map_sim->Get("h_dphi_r_neg");
TH1* h_Z_neg_sim = (TH1*) file_1D_map_sim->Get("h_dz_r_neg");
TH1* h_RP_neg_sim = (TH1*) file_1D_map_sim->Get("h_rdphi_r_neg");
TH1* h_N_pos_sim = new TH1F("h_n_r_pos","",10,20,80);
TH1* h_N_neg_sim = new TH1F("h_n_r_neg","",10,20,80);

h_N_pos_sim->SetLineColor(kBlack);
h_R_pos_sim->SetLineColor(kBlack);
h_P_pos_sim->SetLineColor(kBlack);
h_Z_pos_sim->SetLineColor(kBlack);
h_RP_pos_sim->SetLineColor(kBlack);
h_N_neg_sim->SetLineColor(kBlack);
h_R_neg_sim->SetLineColor(kBlack);
h_P_neg_sim->SetLineColor(kBlack);
h_Z_neg_sim->SetLineColor(kBlack);
h_RP_neg_sim->SetLineColor(kBlack);
h_N_pos_sim->SetTitle(Form("N vs R @ Z=%.0f",selectZ));
h_R_pos_sim->SetTitle(Form("dR vs R @ Z=%.0f",selectZ));
h_P_pos_sim->SetTitle(Form("dphi vs R @ Z=%.0f",selectZ));
h_Z_pos_sim->SetTitle(Form("dz vs R Z=%.0f",selectZ));
h_RP_pos_sim->SetTitle(Form("Rdphi vs R Z=%.0f",selectZ));
h_N_neg_sim->SetTitle(Form("N vs R @ Z=%.0f",selectZ));
h_R_neg_sim->SetTitle(Form("dR vs R @ Z=-%.0f",selectZ));
h_P_neg_sim->SetTitle(Form("dphi vs R @ Z=-%.0f",selectZ));
h_Z_neg_sim->SetTitle(Form("dz vs R Z=-%.0f",selectZ));
h_RP_neg_sim->SetTitle(Form("Rdphi vs R Z=-%.0f",selectZ));

TFile* file_2D_map[nrun];
TH2 *h_N_rz_pos[nrun], *h_R_rz_pos[nrun], *h_P_rz_pos[nrun], *h_Z_rz_pos[nrun], *h_RP_rz_pos[nrun];
TH2 *h_N_rz_neg[nrun], *h_R_rz_neg[nrun], *h_P_rz_neg[nrun], *h_Z_rz_neg[nrun], *h_RP_rz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun], *h_RP_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun], *h_RP_neg[nrun];
double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
double ymax_RP=-100, ymin_RP=100;
for (int i = 0; i < nrun; i++)
{
  file_2D_map[i] = new TFile(Form("Rootfiles/Distortions_2D_mm_%d_rz.root",runs[i]),"");
  h_R_rz_pos[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionR_posz");
  h_P_rz_pos[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionP_posz");
  h_Z_rz_pos[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionZ_posz");
  h_RP_rz_pos[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionRP_posz");
  h_N_rz_pos[i] = (TH2*) file_2D_map[i]->Get("hentries_posz");
  h_R_rz_neg[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionR_negz");
  h_P_rz_neg[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionP_negz");
  h_Z_rz_neg[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionZ_negz");
  h_RP_rz_neg[i] = (TH2*) file_2D_map[i]->Get("hIntDistortionRP_negz");
  h_N_rz_neg[i] = (TH2*) file_2D_map[i]->Get("hentries_negz");

  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),"",h_N_rz_pos[i]->GetXaxis()->GetNbins(),h_N_rz_pos[i]->GetXaxis()->GetXmin(),h_N_rz_pos[i]->GetXaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),"",h_R_rz_pos[i]->GetXaxis()->GetNbins(),h_R_rz_pos[i]->GetXaxis()->GetXmin(),h_R_rz_pos[i]->GetXaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),"",h_P_rz_pos[i]->GetXaxis()->GetNbins(),h_P_rz_pos[i]->GetXaxis()->GetXmin(),h_P_rz_pos[i]->GetXaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),"",h_Z_rz_pos[i]->GetXaxis()->GetNbins(),h_Z_rz_pos[i]->GetXaxis()->GetXmin(),h_Z_rz_pos[i]->GetXaxis()->GetXmax());
  h_RP_pos[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),"",h_RP_rz_pos[i]->GetXaxis()->GetNbins(),h_RP_rz_pos[i]->GetXaxis()->GetXmin(),h_RP_rz_pos[i]->GetXaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),"",h_N_rz_neg[i]->GetXaxis()->GetNbins(),h_N_rz_neg[i]->GetXaxis()->GetXmin(),h_N_rz_neg[i]->GetXaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),"",h_R_rz_neg[i]->GetXaxis()->GetNbins(),h_R_rz_neg[i]->GetXaxis()->GetXmin(),h_R_rz_neg[i]->GetXaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),"",h_P_rz_neg[i]->GetXaxis()->GetNbins(),h_P_rz_neg[i]->GetXaxis()->GetXmin(),h_P_rz_neg[i]->GetXaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),"",h_Z_rz_neg[i]->GetXaxis()->GetNbins(),h_Z_rz_neg[i]->GetXaxis()->GetXmin(),h_Z_rz_neg[i]->GetXaxis()->GetXmax());
  h_RP_neg[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),"",h_RP_rz_neg[i]->GetXaxis()->GetNbins(),h_RP_rz_neg[i]->GetXaxis()->GetXmin(),h_RP_rz_neg[i]->GetXaxis()->GetXmax());
  plot1D_Ybin(h_N_rz_pos[i],h_N_pos[i],selectZ);
  plot1D_Ybin(h_R_rz_pos[i],h_R_pos[i],selectZ);
  plot1D_Ybin(h_P_rz_pos[i],h_P_pos[i],selectZ);
  plot1D_Ybin(h_Z_rz_pos[i],h_Z_pos[i],selectZ);
  plot1D_Ybin(h_RP_rz_pos[i],h_RP_pos[i],selectZ);
  plot1D_Ybin(h_N_rz_neg[i],h_N_neg[i],-1*selectZ);
  plot1D_Ybin(h_R_rz_neg[i],h_R_neg[i],-1*selectZ);
  plot1D_Ybin(h_P_rz_neg[i],h_P_neg[i],-1*selectZ);
  plot1D_Ybin(h_Z_rz_neg[i],h_Z_neg[i],-1*selectZ);
  plot1D_Ybin(h_RP_rz_neg[i],h_RP_neg[i],-1*selectZ);
  h_N_pos[i]->SetLineColor(i+2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(i+2);
  h_R_pos[i]->SetLineColor(i+2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(i+2);
  h_P_pos[i]->SetLineColor(i+2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(i+2);
  h_Z_pos[i]->SetLineColor(i+2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(i+2);
  h_RP_pos[i]->SetLineColor(i+2); h_RP_pos[i]->SetLineWidth(1); h_RP_pos[i]->SetFillColor(0); h_RP_pos[i]->SetMarkerColor(i+2);
  h_N_neg[i]->SetLineColor(i+2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_R_neg[i]->SetLineColor(i+2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_P_neg[i]->SetLineColor(i+2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(i+2);
  h_Z_neg[i]->SetLineColor(i+2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(i+2);
  h_RP_neg[i]->SetLineColor(i+2); h_RP_neg[i]->SetLineWidth(1); h_RP_neg[i]->SetFillColor(0); h_RP_neg[i]->SetMarkerColor(i+2);

  if (h_N_neg[i]->GetMaximum()>ymax_N) ymax_N = h_N_neg[i]->GetMaximum();
  if (h_N_pos[i]->GetMaximum()>ymax_N) ymax_N = h_N_pos[i]->GetMaximum();
  if (h_N_neg[i]->GetMinimum()<ymin_N) ymin_N = h_N_neg[i]->GetMinimum();
  if (h_N_pos[i]->GetMinimum()<ymin_N) ymin_N = h_N_pos[i]->GetMinimum();

  if (h_R_neg[i]->GetMaximum()>ymax_R) ymax_R = h_R_neg[i]->GetMaximum();
  if (h_R_pos[i]->GetMaximum()>ymax_R) ymax_R = h_R_pos[i]->GetMaximum();
  if (h_R_neg[i]->GetMinimum()<ymin_R) ymin_R = h_R_neg[i]->GetMinimum();
  if (h_R_pos[i]->GetMinimum()<ymin_R) ymin_R = h_R_pos[i]->GetMinimum();

  if (h_P_neg[i]->GetMaximum()>ymax_P) ymax_P = h_P_neg[i]->GetMaximum();
  if (h_P_pos[i]->GetMaximum()>ymax_P) ymax_P = h_P_pos[i]->GetMaximum();
  if (h_P_neg[i]->GetMinimum()<ymin_P) ymin_P = h_P_neg[i]->GetMinimum();
  if (h_P_pos[i]->GetMinimum()<ymin_P) ymin_P = h_P_pos[i]->GetMinimum();

  if (h_Z_neg[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_neg[i]->GetMaximum();
  if (h_Z_pos[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_pos[i]->GetMaximum();
  if (h_Z_neg[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_neg[i]->GetMinimum();
  if (h_Z_pos[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_pos[i]->GetMinimum();

  if (h_RP_neg[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_neg[i]->GetMaximum();
  if (h_RP_pos[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_pos[i]->GetMaximum();
  if (h_RP_neg[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_neg[i]->GetMinimum();
  if (h_RP_pos[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_pos[i]->GetMinimum();
}

if (ymax_N>0) ymax_N *= 1.1; else ymax_N = 1;
if (ymax_R>0) ymax_R *= 1.1; else ymax_R = 0.1;
if (ymax_P>0) ymax_P *= 1.1; else ymax_P = 0.1;
if (ymax_Z>0) ymax_Z *= 1.1; else ymax_Z = 0.1;
if (ymax_RP>0) ymax_RP *= 1.1; else ymax_RP = 0.1;

if (ymin_N<0) ymin_N *= 1.1; else ymin_N = 0.9;
if (ymin_R<0) ymin_R *= 1.1; else ymin_R = -0.1;
if (ymin_P<0) ymin_P *= 1.1; else ymin_P = -0.1;
if (ymin_Z<0) ymin_Z *= 1.1; else ymin_Z = -0.1;
if (ymin_RP<0) ymin_RP *= 1.1; else ymin_RP = -0.1;

h_N_pos_sim->SetMinimum(ymin_N); h_N_pos_sim->SetMaximum(ymax_N);
h_R_pos_sim->SetMinimum(ymin_R); h_R_pos_sim->SetMaximum(ymax_R);
h_P_pos_sim->SetMinimum(ymin_P); h_P_pos_sim->SetMaximum(ymax_P);
h_Z_pos_sim->SetMinimum(ymin_Z); h_Z_pos_sim->SetMaximum(ymax_Z);
h_RP_pos_sim->SetMinimum(ymin_RP); h_RP_pos_sim->SetMaximum(ymax_RP);

h_N_neg_sim->SetMinimum(ymin_N); h_N_neg_sim->SetMaximum(ymax_N);
h_R_neg_sim->SetMinimum(ymin_R); h_R_neg_sim->SetMaximum(ymax_R);
h_P_neg_sim->SetMinimum(ymin_P); h_P_neg_sim->SetMaximum(ymax_P);
h_Z_neg_sim->SetMinimum(ymin_Z); h_Z_neg_sim->SetMaximum(ymax_Z);
h_RP_neg_sim->SetMinimum(ymin_RP); h_RP_neg_sim->SetMaximum(ymax_RP);

TCanvas* can = new TCanvas("can","",4000,1200);
can->Divide(5,2);
can->cd(1);
h_P_pos_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_P_pos[i]->Draw("hist,e,same"); h_P_pos_sim->Draw("hist,same");
can->cd(2);
h_RP_pos_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_RP_pos[i]->Draw("hist,e,same"); h_RP_pos_sim->Draw("hist,same");
can->cd(3);
h_R_pos_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_R_pos[i]->Draw("hist,e,same"); h_R_pos_sim->Draw("hist,same");
can->cd(4);
h_Z_pos_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_Z_pos[i]->Draw("hist,e,same"); h_Z_pos_sim->Draw("hist,same");
can->cd(5);
gPad->SetLogy(1);
h_N_pos_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_N_pos[i]->Draw("hist,e,same"); h_N_pos_sim->Draw("hist,same");
can->cd(6);
h_P_neg_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_P_neg[i]->Draw("hist,e,same"); h_P_neg_sim->Draw("hist,same");
can->cd(7);
h_RP_neg_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_RP_neg[i]->Draw("hist,e,same"); h_RP_neg_sim->Draw("hist,same");
can->cd(8);
h_R_neg_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_R_neg[i]->Draw("hist,e,same"); h_R_neg_sim->Draw("hist,same");
can->cd(9);
h_Z_neg_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_Z_neg[i]->Draw("hist,e,same"); h_Z_neg_sim->Draw("hist,same");
can->cd(10);
gPad->SetLogy(1);
h_N_neg_sim->Draw("hist"); for (int i=0; i<nrun; i++) h_N_neg[i]->Draw("hist,e,same"); h_N_neg_sim->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_radius_from2D_atZ%d.pdf",(int)selectZ));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->AddEntry(h_P_pos_sim, "Default simulation distortion map", "l");
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos[i], Form("Run %d, MBD rate %d kHz",runs[i],mbdrates[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_radius_from2D_leg.pdf"));

for (int i=0; i<nrun; i++)
{
  TFile* outfile = new TFile(Form("Rootfiles/Distortions_vsR_atZ%d_%d.root",(int)selectZ,runs[i]),"recreate");
  outfile->cd();
  h_N_pos[i]->Write("hentries_posz");
  h_R_pos[i]->Write("hIntDistortionR_posz");
  h_P_pos[i]->Write("hIntDistortionP_posz");
  h_Z_pos[i]->Write("hIntDistortionZ_posz");
  h_RP_pos[i]->Write("hIntDistortionRP_posz");
  h_N_neg[i]->Write("hentries_negz");
  h_R_neg[i]->Write("hIntDistortionR_negz");
  h_P_neg[i]->Write("hIntDistortionP_negz");
  h_Z_neg[i]->Write("hIntDistortionZ_negz");
  h_RP_neg[i]->Write("hIntDistortionRP_negz");
  outfile->Write();
  outfile->Close();
}
delete can;
delete can_leg;
}

}
