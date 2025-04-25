void plot2D_Zbin(TH3* h3, TH2* h2, float z)
{
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h3->GetNbinsY(); j++)
    {
      if (isnan(h3->GetBinContent(i, j, zbin)))
      {
        h2->SetBinContent(i, j, 0.0);
	h2->SetBinError(i, j, 0.0);
      }
      else
      {
        h2->SetBinContent(i, j, h3->GetBinContent(i, j, zbin));
        h2->SetBinError(i, j, h3->GetBinError(i, j, zbin));
      }
    }
  }
}

void SetTH2ZeroBin(TH2* h2)
{
  for (int i = 1; i <= h2->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h2->GetNbinsY(); j++)
    {
      if (isnan(h2->GetBinContent(i, j)))
      {
        h2->SetBinContent(i, j, -1000);
	h2->SetBinError(i, j, -1000);
      }
      else if (h2->GetBinContent(i, j)==0)
      {
        h2->SetBinContent(i, j, -1000);
        h2->SetBinError(i, j, -1000);
      }
    }
  }
}

void draw2D_phir_from3D()
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
TH2 *h_N_pr_pos[nrun], *h_R_pr_pos[nrun], *h_P_pr_pos[nrun], *h_Z_pr_pos[nrun], *h_RP_pr_pos[nrun];
TH2 *h_N_pr_neg[nrun], *h_R_pr_neg[nrun], *h_P_pr_neg[nrun], *h_Z_pr_neg[nrun], *h_RP_pr_neg[nrun];
for (int i = 0; i < nrun; i++)
{
  //file_3D_map[i] = new TFile(Form("Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
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

  h_N_pr_pos[i] = new TH2F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;#phi (rad);R (cm);N",(int)selectZ),
		  h_N_prz_pos[i]->GetXaxis()->GetNbins(),h_N_prz_pos[i]->GetXaxis()->GetXmin(),h_N_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pr_pos[i] = new TH2F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;#phi (rad);R (cm);dR (cm)",(int)selectZ),
		  h_R_prz_pos[i]->GetXaxis()->GetNbins(),h_R_prz_pos[i]->GetXaxis()->GetXmin(),h_R_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pr_pos[i] = new TH2F(Form("hIntDistortionP_posz_%d",runs[i]),Form("d#phi @ Z=%d cm;#phi (rad);R (cm);d#phi (rad)",(int)selectZ),
		  h_P_prz_pos[i]->GetXaxis()->GetNbins(),h_P_prz_pos[i]->GetXaxis()->GetXmin(),h_P_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pr_pos[i] = new TH2F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;#phi (rad);R (cm);dz (cm)",(int)selectZ),
		  h_Z_prz_pos[i]->GetXaxis()->GetNbins(),h_Z_prz_pos[i]->GetXaxis()->GetXmin(),h_Z_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_RP_pr_pos[i] = new TH2F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rd#phi @ Z=%d cm;#phi (rad);R (cm);Rd#phi (cm)",(int)selectZ),
		  h_RP_prz_pos[i]->GetXaxis()->GetNbins(),h_RP_prz_pos[i]->GetXaxis()->GetXmin(),h_RP_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_RP_prz_pos[i]->GetYaxis()->GetNbins(),h_RP_prz_pos[i]->GetYaxis()->GetXmin(),h_RP_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_pr_neg[i] = new TH2F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;#phi (rad);R (cm);N",(int)selectZ),
		  h_N_prz_neg[i]->GetXaxis()->GetNbins(),h_N_prz_neg[i]->GetXaxis()->GetXmin(),h_N_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_pr_neg[i] = new TH2F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;#phi (rad);R (cm);dR (cm)",(int)selectZ),
		  h_R_prz_neg[i]->GetXaxis()->GetNbins(),h_R_prz_neg[i]->GetXaxis()->GetXmin(),h_R_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_pr_neg[i] = new TH2F(Form("hIntDistortionP_negz_%d",runs[i]),Form("d#phi @ Z=-%d cm;#phi (rad);R (cm);d#phi (rad)",(int)selectZ),
		  h_P_prz_neg[i]->GetXaxis()->GetNbins(),h_P_prz_neg[i]->GetXaxis()->GetXmin(),h_P_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_pr_neg[i] = new TH2F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;#phi (rad);R (cm);dz (cm)",(int)selectZ),
		  h_Z_prz_neg[i]->GetXaxis()->GetNbins(),h_Z_prz_neg[i]->GetXaxis()->GetXmin(),h_Z_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());
  h_RP_pr_neg[i] = new TH2F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rd#phi @ Z=-%d cm;#phi (rad);R (cm);Rd#phi (cm)",(int)selectZ),
		  h_RP_prz_neg[i]->GetXaxis()->GetNbins(),h_RP_prz_neg[i]->GetXaxis()->GetXmin(),h_RP_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_RP_prz_neg[i]->GetYaxis()->GetNbins(),h_RP_prz_neg[i]->GetYaxis()->GetXmin(),h_RP_prz_neg[i]->GetYaxis()->GetXmax());
  plot2D_Zbin(h_N_prz_pos[i],h_N_pr_pos[i],selectZ);
  plot2D_Zbin(h_R_prz_pos[i],h_R_pr_pos[i],selectZ);
  plot2D_Zbin(h_P_prz_pos[i],h_P_pr_pos[i],selectZ);
  plot2D_Zbin(h_Z_prz_pos[i],h_Z_pr_pos[i],selectZ);
  plot2D_Zbin(h_RP_prz_pos[i],h_RP_pr_pos[i],selectZ);
  plot2D_Zbin(h_N_prz_neg[i],h_N_pr_neg[i],-selectZ);
  plot2D_Zbin(h_R_prz_neg[i],h_R_pr_neg[i],-selectZ);
  plot2D_Zbin(h_P_prz_neg[i],h_P_pr_neg[i],-selectZ);
  plot2D_Zbin(h_Z_prz_neg[i],h_Z_pr_neg[i],-selectZ);
  plot2D_Zbin(h_RP_prz_neg[i],h_RP_pr_neg[i],-selectZ);

  double zmax_N=-100, zmin_N=100;
  double zmax_R=-100, zmin_R=100;
  double zmax_P=-100, zmin_P=100;
  double zmax_Z=-100, zmin_Z=100;
  double zmax_RP=-100, zmin_RP=100;

  if (h_N_pr_neg[i]->GetMaximum()>zmax_N) zmax_N = h_N_pr_neg[i]->GetMaximum();
  if (h_N_pr_pos[i]->GetMaximum()>zmax_N) zmax_N = h_N_pr_pos[i]->GetMaximum();
  //if (h_N_pr_neg[i]->GetMinimum()<zmin_N) zmin_N = h_N_pr_neg[i]->GetMinimum();
  //if (h_N_pr_pos[i]->GetMinimum()<zmin_N) zmin_N = h_N_pr_pos[i]->GetMinimum();

  if (h_R_pr_neg[i]->GetMaximum()>zmax_R) zmax_R = h_R_pr_neg[i]->GetMaximum();
  if (h_R_pr_pos[i]->GetMaximum()>zmax_R) zmax_R = h_R_pr_pos[i]->GetMaximum();
  if (h_R_pr_neg[i]->GetMinimum()<zmin_R) zmin_R = h_R_pr_neg[i]->GetMinimum();
  if (h_R_pr_pos[i]->GetMinimum()<zmin_R) zmin_R = h_R_pr_pos[i]->GetMinimum();

  if (h_P_pr_neg[i]->GetMaximum()>zmax_P) zmax_P = h_P_pr_neg[i]->GetMaximum();
  if (h_P_pr_pos[i]->GetMaximum()>zmax_P) zmax_P = h_P_pr_pos[i]->GetMaximum();
  if (h_P_pr_neg[i]->GetMinimum()<zmin_P) zmin_P = h_P_pr_neg[i]->GetMinimum();
  if (h_P_pr_pos[i]->GetMinimum()<zmin_P) zmin_P = h_P_pr_pos[i]->GetMinimum();

  if (h_Z_pr_neg[i]->GetMaximum()>zmax_Z) zmax_Z = h_Z_pr_neg[i]->GetMaximum();
  if (h_Z_pr_pos[i]->GetMaximum()>zmax_Z) zmax_Z = h_Z_pr_pos[i]->GetMaximum();
  if (h_Z_pr_neg[i]->GetMinimum()<zmin_Z) zmin_Z = h_Z_pr_neg[i]->GetMinimum();
  if (h_Z_pr_pos[i]->GetMinimum()<zmin_Z) zmin_Z = h_Z_pr_pos[i]->GetMinimum();

  if (h_RP_pr_neg[i]->GetMaximum()>zmax_RP) zmax_RP = h_RP_pr_neg[i]->GetMaximum();
  if (h_RP_pr_pos[i]->GetMaximum()>zmax_RP) zmax_RP = h_RP_pr_pos[i]->GetMaximum();
  if (h_RP_pr_neg[i]->GetMinimum()<zmin_RP) zmin_RP = h_RP_pr_neg[i]->GetMinimum();
  if (h_RP_pr_pos[i]->GetMinimum()<zmin_RP) zmin_RP = h_RP_pr_pos[i]->GetMinimum();

  if (zmax_N>0) zmax_N *= 1.1; else zmax_N = 1;
  if (zmax_R>0) zmax_R *= 1.1; else zmax_R = 0.1;
  if (zmax_P>0) zmax_P *= 1.1; else zmax_P = 0.1;
  if (zmax_Z>0) zmax_Z *= 1.1; else zmax_Z = 0.1;
  if (zmax_RP>0) zmax_RP *= 1.1; else zmax_RP = 0.1;

  //if (zmin_N<0) zmin_N *= 1.1; else zmin_N = 0.9;
  if (zmin_R<0) zmin_R *= 1.1; else zmin_R = -0.1;
  if (zmin_P<0) zmin_P *= 1.1; else zmin_P = -0.1;
  if (zmin_Z<0) zmin_Z *= 1.1; else zmin_Z = -0.1;
  if (zmin_RP<0) zmin_RP *= 1.1; else zmin_RP = -0.1;

  h_N_pr_pos[i]->SetMinimum(0); h_N_pr_pos[i]->SetMaximum(zmax_N);
  h_R_pr_pos[i]->SetMinimum(zmin_R); h_R_pr_pos[i]->SetMaximum(zmax_R);
  h_P_pr_pos[i]->SetMinimum(zmin_P); h_P_pr_pos[i]->SetMaximum(zmax_P);
  h_Z_pr_pos[i]->SetMinimum(zmin_Z); h_Z_pr_pos[i]->SetMaximum(zmax_Z);
  h_RP_pr_pos[i]->SetMinimum(zmin_RP); h_RP_pr_pos[i]->SetMaximum(zmax_RP);

  h_N_pr_neg[i]->SetMinimum(0); h_N_pr_neg[i]->SetMaximum(zmax_N);
  h_R_pr_neg[i]->SetMinimum(zmin_R); h_R_pr_neg[i]->SetMaximum(zmax_R);
  h_P_pr_neg[i]->SetMinimum(zmin_P); h_P_pr_neg[i]->SetMaximum(zmax_P);
  h_Z_pr_neg[i]->SetMinimum(zmin_Z); h_Z_pr_neg[i]->SetMaximum(zmax_Z);
  h_RP_pr_neg[i]->SetMinimum(zmin_RP); h_RP_pr_neg[i]->SetMaximum(zmax_RP);

  SetTH2ZeroBin(h_N_pr_pos[i]);
  SetTH2ZeroBin(h_R_pr_pos[i]);
  SetTH2ZeroBin(h_P_pr_pos[i]);
  SetTH2ZeroBin(h_Z_pr_pos[i]);
  SetTH2ZeroBin(h_RP_pr_pos[i]);
  SetTH2ZeroBin(h_N_pr_neg[i]);
  SetTH2ZeroBin(h_R_pr_neg[i]);
  SetTH2ZeroBin(h_P_pr_neg[i]);
  SetTH2ZeroBin(h_Z_pr_neg[i]);
  SetTH2ZeroBin(h_RP_pr_neg[i]);

  TCanvas* can = new TCanvas("can","",8000,2400);
  can->Divide(5,2);
  can->cd(1); gPad->SetRightMargin(0.15);
  h_P_pr_pos[i]->Draw("colz");
  can->cd(2); gPad->SetRightMargin(0.15);
  h_RP_pr_pos[i]->Draw("colz");
  can->cd(3); gPad->SetRightMargin(0.15);
  h_R_pr_pos[i]->Draw("colz");
  can->cd(4); gPad->SetRightMargin(0.15);
  h_Z_pr_pos[i]->Draw("colz");
  can->cd(5); gPad->SetRightMargin(0.15);
  h_N_pr_pos[i]->Draw("colz");
  can->cd(6); gPad->SetRightMargin(0.15);
  h_P_pr_neg[i]->Draw("colz");
  can->cd(7); gPad->SetRightMargin(0.15);
  h_RP_pr_neg[i]->Draw("colz");
  can->cd(8); gPad->SetRightMargin(0.15);
  h_R_pr_neg[i]->Draw("colz");
  can->cd(9); gPad->SetRightMargin(0.15);
  h_Z_pr_neg[i]->Draw("colz");
  can->cd(10); gPad->SetRightMargin(0.15);
  h_N_pr_neg[i]->Draw("colz");
  
  gPad->RedrawAxis();
  
  can->Update();
  can->SaveAs(Form("figure/resid_phir_from3D_%d_atZ%d.pdf",runs[i],(int)selectZ));
  delete can;

}


for (int i=0; i<nrun; i++)
{
  TFile* outfile = new TFile(Form("Rootfiles/Distortions_vsPR_atR%d_%d.root",(int)selectZ,runs[i]),"recreate");
  outfile->cd();
  h_N_pr_pos[i]->Write("hentries_posz");
  h_R_pr_pos[i]->Write("hIntDistortionR_posz");
  h_P_pr_pos[i]->Write("hIntDistortionP_posz");
  h_Z_pr_pos[i]->Write("hIntDistortionZ_posz");
  h_RP_pr_pos[i]->Write("hIntDistortionRP_posz");
  h_N_pr_neg[i]->Write("hentries_negz");
  h_R_pr_neg[i]->Write("hIntDistortionR_negz");
  h_P_pr_neg[i]->Write("hIntDistortionP_negz");
  h_Z_pr_neg[i]->Write("hIntDistortionZ_negz");
  h_RP_pr_neg[i]->Write("hIntDistortionRP_negz");
  outfile->Write();
  outfile->Close();
}
}

}
