void plot1D_Zbin(TH3* h3, TH1* h1, float z)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
	  //cout<<"pbin = "<<xbin<<", rbin = "<<i<<" , zbin = "<<zbin<<" , content = "<<h3->GetBinContent(xbin, i, zbin)<<endl;
      h1->SetBinContent(i, h3->GetBinContent(xbin, i, zbin));
      h1->SetBinError(i, h3->GetBinError(xbin, i, zbin));
  }
}

void draw1D_r_from3D()
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
//int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
//int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
//int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
//int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
int mbdrates[nrun] = {250};
int runs[nrun] = {53534};

TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun], *h_RP_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun], *h_RP_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun], *h_RP_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun], *h_RP_neg[nrun];
double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
double ymax_RP=-100, ymin_RP=100;
for (int i = 0; i < nrun; i++)
{
  file_3D_map[i] = new TFile(Form("Rootfiles/Distortions_2D_mm_%d_rz.root",runs[i]),"");
  h_R_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_posz");
  h_RP_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_posz");
  h_N_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hentries_posz");
  h_R_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_negz");
  h_RP_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_negz");
  h_N_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hentries_negz");

  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_RP_pos[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_RP_prz_pos[i]->GetYaxis()->GetNbins(),h_RP_prz_pos[i]->GetYaxis()->GetXmin(),h_RP_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());
  h_RP_neg[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_RP_prz_neg[i]->GetYaxis()->GetNbins(),h_RP_prz_neg[i]->GetYaxis()->GetXmin(),h_RP_prz_neg[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos[i],h_N_pos[i],selectZ);
  plot1D_Zbin(h_R_prz_pos[i],h_R_pos[i],selectZ);
  plot1D_Zbin(h_P_prz_pos[i],h_P_pos[i],selectZ);
  plot1D_Zbin(h_Z_prz_pos[i],h_Z_pos[i],selectZ);
  plot1D_Zbin(h_RP_prz_pos[i],h_RP_pos[i],selectZ);
  plot1D_Zbin(h_N_prz_neg[i],h_N_neg[i],-selectZ);
  plot1D_Zbin(h_R_prz_neg[i],h_R_neg[i],-selectZ);
  plot1D_Zbin(h_P_prz_neg[i],h_P_neg[i],-selectZ);
  plot1D_Zbin(h_Z_prz_neg[i],h_Z_neg[i],-selectZ);
  plot1D_Zbin(h_RP_prz_neg[i],h_RP_neg[i],-selectZ);
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
if (ymax_R>0) ymax_R *= 1.1; else ymax_R = 0;
if (ymax_P>0) ymax_P *= 1.1; else ymax_P = 0;
if (ymax_Z>0) ymax_Z *= 1.1; else ymax_Z = 0;
if (ymax_RP>0) ymax_RP *= 1.1; else ymax_RP = 0;

if (ymin_N<0) ymin_N *= 1.1; else ymin_N = 0.9;
if (ymin_R<0) ymin_R *= 1.1; else ymin_R = 0.0;
if (ymin_P<0) ymin_P *= 1.1; else ymin_P = 0.0;
if (ymin_Z<0) ymin_Z *= 1.1; else ymin_Z = 0.0;
if (ymin_RP<0) ymin_RP *= 1.1; else ymin_RP = 0.0;

for (int i = 0; i < nrun; i++)
{
  h_N_pos[i]->SetMinimum(ymin_N); h_N_pos[i]->SetMaximum(ymax_N);
  h_R_pos[i]->SetMinimum(ymin_R); h_R_pos[i]->SetMaximum(ymax_R);
  h_P_pos[i]->SetMinimum(ymin_P); h_P_pos[i]->SetMaximum(ymax_P);
  h_Z_pos[i]->SetMinimum(ymin_Z); h_Z_pos[i]->SetMaximum(ymax_Z);
  h_RP_pos[i]->SetMinimum(ymin_RP); h_RP_pos[i]->SetMaximum(ymax_RP);

  h_N_neg[i]->SetMinimum(ymin_N); h_N_neg[i]->SetMaximum(ymax_N);
  h_R_neg[i]->SetMinimum(ymin_R); h_R_neg[i]->SetMaximum(ymax_R);
  h_P_neg[i]->SetMinimum(ymin_P); h_P_neg[i]->SetMaximum(ymax_P);
  h_Z_neg[i]->SetMinimum(ymin_Z); h_Z_neg[i]->SetMaximum(ymax_Z);
  h_RP_neg[i]->SetMinimum(ymin_RP); h_RP_neg[i]->SetMaximum(ymax_RP);
}

TCanvas* can = new TCanvas("can","",4000,1200);
can->Divide(5,2);
can->cd(1);
for (int i=0; i<nrun; i++) h_P_pos[i]->Draw("hist,e,same");
can->cd(2);
for (int i=0; i<nrun; i++) h_RP_pos[i]->Draw("hist,e,same");
can->cd(3);
for (int i=0; i<nrun; i++) h_R_pos[i]->Draw("hist,e,same");
can->cd(4);
for (int i=0; i<nrun; i++) h_Z_pos[i]->Draw("hist,e,same");
can->cd(5);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_pos[i]->Draw("hist,e,same");
can->cd(6);
for (int i=0; i<nrun; i++) h_P_neg[i]->Draw("hist,e,same");
can->cd(7);
for (int i=0; i<nrun; i++) h_RP_neg[i]->Draw("hist,e,same");
can->cd(8);
for (int i=0; i<nrun; i++) h_R_neg[i]->Draw("hist,e,same");
can->cd(9);
for (int i=0; i<nrun; i++) h_Z_neg[i]->Draw("hist,e,same");
can->cd(10);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_neg[i]->Draw("hist,e,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d.pdf",(int)selectZ));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos[i], Form("Run %d, MBD rate %d kHz",runs[i],mbdrates[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsR_from3D_leg.pdf"));

/*
for (int i=0; i<nrun; i++)
{
  TFile* outfile = new TFile(Form("Rootfiles/Distortions_vsR_atR%d_%d.root",(int)selectR,runs[i]),"recreate");
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
*/

delete can;
delete can_leg;
}

}
