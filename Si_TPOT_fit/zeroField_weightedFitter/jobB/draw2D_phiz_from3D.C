void plot2D_Ybin(TH3* h3, TH2* h2, float y)
{
  int ybin = h3->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h3->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h3->GetNbinsZ(); j++)
    {
      if (isnan(h3->GetBinContent(i, ybin, j)))
      {
        h2->SetBinContent(i, j, 0.0);
	h2->SetBinError(i, j, 0.0);
      }
      else
      {
        h2->SetBinContent(i, j, h3->GetBinContent(i, ybin, j));
        h2->SetBinError(i, j, h3->GetBinError(i, ybin, j));
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

void draw2D_phiz_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
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
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun], *h_RP_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun], *h_RP_prz_neg[nrun];
TH2 *h_N_pz_pos[nrun], *h_R_pz_pos[nrun], *h_P_pz_pos[nrun], *h_Z_pz_pos[nrun], *h_RP_pz_pos[nrun];
TH2 *h_N_pz_neg[nrun], *h_R_pz_neg[nrun], *h_P_pz_neg[nrun], *h_Z_pz_neg[nrun], *h_RP_pz_neg[nrun];
for (int i = 0; i < nrun; i++)
{
  file_3D_map[i] = new TFile(Form("Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  //file_3D_map[i] = new TFile(Form("Rootfiles/Distortions_2D_mm_%d_rz.root",runs[i]),"");
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

  h_N_pz_pos[i] = new TH2F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;#phi (rad);Z (cm);N",(int)selectR),
		  h_N_prz_pos[i]->GetXaxis()->GetNbins(),h_N_prz_pos[i]->GetXaxis()->GetXmin(),h_N_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_N_prz_pos[i]->GetZaxis()->GetNbins(),h_N_prz_pos[i]->GetZaxis()->GetXmin(),h_N_prz_pos[i]->GetZaxis()->GetXmax());
  h_R_pz_pos[i] = new TH2F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;#phi (rad);Z (cm);dR (cm)",(int)selectR),
		  h_R_prz_pos[i]->GetXaxis()->GetNbins(),h_R_prz_pos[i]->GetXaxis()->GetXmin(),h_R_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_R_prz_pos[i]->GetZaxis()->GetNbins(),h_R_prz_pos[i]->GetZaxis()->GetXmin(),h_R_prz_pos[i]->GetZaxis()->GetXmax());
  h_P_pz_pos[i] = new TH2F(Form("hIntDistortionP_posz_%d",runs[i]),Form("d#phi @ R=%d cm;#phi (rad);Z (cm);d#phi (rad)",(int)selectR),
		  h_P_prz_pos[i]->GetXaxis()->GetNbins(),h_P_prz_pos[i]->GetXaxis()->GetXmin(),h_P_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_P_prz_pos[i]->GetZaxis()->GetNbins(),h_P_prz_pos[i]->GetZaxis()->GetXmin(),h_P_prz_pos[i]->GetZaxis()->GetXmax());
  h_Z_pz_pos[i] = new TH2F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;#phi (rad);Z (cm);dz (cm)",(int)selectR),
		  h_Z_prz_pos[i]->GetXaxis()->GetNbins(),h_Z_prz_pos[i]->GetXaxis()->GetXmin(),h_Z_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_Z_prz_pos[i]->GetZaxis()->GetNbins(),h_Z_prz_pos[i]->GetZaxis()->GetXmin(),h_Z_prz_pos[i]->GetZaxis()->GetXmax());
  h_RP_pz_pos[i] = new TH2F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rd#phi @ R=%d cm;#phi (rad);Z (cm);Rd#phi (cm)",(int)selectR),
		  h_RP_prz_pos[i]->GetXaxis()->GetNbins(),h_RP_prz_pos[i]->GetXaxis()->GetXmin(),h_RP_prz_pos[i]->GetXaxis()->GetXmax(),
		  h_RP_prz_pos[i]->GetZaxis()->GetNbins(),h_RP_prz_pos[i]->GetZaxis()->GetXmin(),h_RP_prz_pos[i]->GetZaxis()->GetXmax());
  h_N_pz_neg[i] = new TH2F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;#phi (rad);Z (cm);N",(int)selectR),
		  h_N_prz_neg[i]->GetXaxis()->GetNbins(),h_N_prz_neg[i]->GetXaxis()->GetXmin(),h_N_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_N_prz_neg[i]->GetZaxis()->GetNbins(),h_N_prz_neg[i]->GetZaxis()->GetXmin(),h_N_prz_neg[i]->GetZaxis()->GetXmax());
  h_R_pz_neg[i] = new TH2F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;#phi (rad);Z (cm);dR (cm)",(int)selectR),
		  h_R_prz_neg[i]->GetXaxis()->GetNbins(),h_R_prz_neg[i]->GetXaxis()->GetXmin(),h_R_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_R_prz_neg[i]->GetZaxis()->GetNbins(),h_R_prz_neg[i]->GetZaxis()->GetXmin(),h_R_prz_neg[i]->GetZaxis()->GetXmax());
  h_P_pz_neg[i] = new TH2F(Form("hIntDistortionP_negz_%d",runs[i]),Form("d#phi @ R=%d cm;#phi (rad);Z (cm);d#phi (rad)",(int)selectR),
		  h_P_prz_neg[i]->GetXaxis()->GetNbins(),h_P_prz_neg[i]->GetXaxis()->GetXmin(),h_P_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_P_prz_neg[i]->GetZaxis()->GetNbins(),h_P_prz_neg[i]->GetZaxis()->GetXmin(),h_P_prz_neg[i]->GetZaxis()->GetXmax());
  h_Z_pz_neg[i] = new TH2F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;#phi (rad);Z (cm);dz (cm)",(int)selectR),
		  h_Z_prz_neg[i]->GetXaxis()->GetNbins(),h_Z_prz_neg[i]->GetXaxis()->GetXmin(),h_Z_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_Z_prz_neg[i]->GetZaxis()->GetNbins(),h_Z_prz_neg[i]->GetZaxis()->GetXmin(),h_Z_prz_neg[i]->GetZaxis()->GetXmax());
  h_RP_pz_neg[i] = new TH2F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rd#phi @ R=%d cm;#phi (rad);Z (cm);Rd#phi (cm)",(int)selectR),
		  h_RP_prz_neg[i]->GetXaxis()->GetNbins(),h_RP_prz_neg[i]->GetXaxis()->GetXmin(),h_RP_prz_neg[i]->GetXaxis()->GetXmax(),
		  h_RP_prz_neg[i]->GetZaxis()->GetNbins(),h_RP_prz_neg[i]->GetZaxis()->GetXmin(),h_RP_prz_neg[i]->GetZaxis()->GetXmax());
  plot2D_Ybin(h_N_prz_pos[i],h_N_pz_pos[i],selectR);
  plot2D_Ybin(h_R_prz_pos[i],h_R_pz_pos[i],selectR);
  plot2D_Ybin(h_P_prz_pos[i],h_P_pz_pos[i],selectR);
  plot2D_Ybin(h_Z_prz_pos[i],h_Z_pz_pos[i],selectR);
  plot2D_Ybin(h_RP_prz_pos[i],h_RP_pz_pos[i],selectR);
  plot2D_Ybin(h_N_prz_neg[i],h_N_pz_neg[i],selectR);
  std::cout<<"print h_R_prz_neg "<<runs[i]<<std::endl;
  plot2D_Ybin(h_R_prz_neg[i],h_R_pz_neg[i],selectR);
  plot2D_Ybin(h_P_prz_neg[i],h_P_pz_neg[i],selectR);
  plot2D_Ybin(h_Z_prz_neg[i],h_Z_pz_neg[i],selectR);
  plot2D_Ybin(h_RP_prz_neg[i],h_RP_pz_neg[i],selectR);

  double zmax_N=-100, zmin_N=100;
  double zmax_R=-100, zmin_R=100;
  double zmax_P=-100, zmin_P=100;
  double zmax_Z=-100, zmin_Z=100;
  double zmax_RP=-100, zmin_RP=100;

  if (h_N_pz_neg[i]->GetMaximum()>zmax_N) zmax_N = h_N_pz_neg[i]->GetMaximum();
  if (h_N_pz_pos[i]->GetMaximum()>zmax_N) zmax_N = h_N_pz_pos[i]->GetMaximum();
  //if (h_N_pz_neg[i]->GetMinimum()<zmin_N) zmin_N = h_N_pz_neg[i]->GetMinimum();
  //if (h_N_pz_pos[i]->GetMinimum()<zmin_N) zmin_N = h_N_pz_pos[i]->GetMinimum();

  if (h_R_pz_neg[i]->GetMaximum()>zmax_R) zmax_R = h_R_pz_neg[i]->GetMaximum();
  if (h_R_pz_pos[i]->GetMaximum()>zmax_R) zmax_R = h_R_pz_pos[i]->GetMaximum();
  if (h_R_pz_neg[i]->GetMinimum()<zmin_R) zmin_R = h_R_pz_neg[i]->GetMinimum();
  if (h_R_pz_pos[i]->GetMinimum()<zmin_R) zmin_R = h_R_pz_pos[i]->GetMinimum();

  if (h_P_pz_neg[i]->GetMaximum()>zmax_P) zmax_P = h_P_pz_neg[i]->GetMaximum();
  if (h_P_pz_pos[i]->GetMaximum()>zmax_P) zmax_P = h_P_pz_pos[i]->GetMaximum();
  if (h_P_pz_neg[i]->GetMinimum()<zmin_P) zmin_P = h_P_pz_neg[i]->GetMinimum();
  if (h_P_pz_pos[i]->GetMinimum()<zmin_P) zmin_P = h_P_pz_pos[i]->GetMinimum();

  if (h_Z_pz_neg[i]->GetMaximum()>zmax_Z) zmax_Z = h_Z_pz_neg[i]->GetMaximum();
  if (h_Z_pz_pos[i]->GetMaximum()>zmax_Z) zmax_Z = h_Z_pz_pos[i]->GetMaximum();
  if (h_Z_pz_neg[i]->GetMinimum()<zmin_Z) zmin_Z = h_Z_pz_neg[i]->GetMinimum();
  if (h_Z_pz_pos[i]->GetMinimum()<zmin_Z) zmin_Z = h_Z_pz_pos[i]->GetMinimum();

  if (h_RP_pz_neg[i]->GetMaximum()>zmax_RP) zmax_RP = h_RP_pz_neg[i]->GetMaximum();
  if (h_RP_pz_pos[i]->GetMaximum()>zmax_RP) zmax_RP = h_RP_pz_pos[i]->GetMaximum();
  if (h_RP_pz_neg[i]->GetMinimum()<zmin_RP) zmin_RP = h_RP_pz_neg[i]->GetMinimum();
  if (h_RP_pz_pos[i]->GetMinimum()<zmin_RP) zmin_RP = h_RP_pz_pos[i]->GetMinimum();

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

  h_N_pz_pos[i]->SetMinimum(0); h_N_pz_pos[i]->SetMaximum(zmax_N);
  h_R_pz_pos[i]->SetMinimum(zmin_R); h_R_pz_pos[i]->SetMaximum(zmax_R);
  h_P_pz_pos[i]->SetMinimum(zmin_P); h_P_pz_pos[i]->SetMaximum(zmax_P);
  h_Z_pz_pos[i]->SetMinimum(zmin_Z); h_Z_pz_pos[i]->SetMaximum(zmax_Z);
  h_RP_pz_pos[i]->SetMinimum(zmin_RP); h_RP_pz_pos[i]->SetMaximum(zmax_RP);

  h_N_pz_neg[i]->SetMinimum(0); h_N_pz_neg[i]->SetMaximum(zmax_N);
  h_R_pz_neg[i]->SetMinimum(zmin_R); h_R_pz_neg[i]->SetMaximum(zmax_R);
  h_P_pz_neg[i]->SetMinimum(zmin_P); h_P_pz_neg[i]->SetMaximum(zmax_P);
  h_Z_pz_neg[i]->SetMinimum(zmin_Z); h_Z_pz_neg[i]->SetMaximum(zmax_Z);
  h_RP_pz_neg[i]->SetMinimum(zmin_RP); h_RP_pz_neg[i]->SetMaximum(zmax_RP);

  SetTH2ZeroBin(h_N_pz_pos[i]);
  SetTH2ZeroBin(h_R_pz_pos[i]);
  SetTH2ZeroBin(h_P_pz_pos[i]);
  SetTH2ZeroBin(h_Z_pz_pos[i]);
  SetTH2ZeroBin(h_RP_pz_pos[i]);
  SetTH2ZeroBin(h_N_pz_neg[i]);
  SetTH2ZeroBin(h_R_pz_neg[i]);
  SetTH2ZeroBin(h_P_pz_neg[i]);
  SetTH2ZeroBin(h_Z_pz_neg[i]);
  SetTH2ZeroBin(h_RP_pz_neg[i]);

  TCanvas* can = new TCanvas("can","",8000,2400);
  can->Divide(5,2);
  can->cd(1); gPad->SetRightMargin(0.15);
  h_P_pz_pos[i]->Draw("colz");
  can->cd(2); gPad->SetRightMargin(0.15);
  h_RP_pz_pos[i]->Draw("colz");
  can->cd(3); gPad->SetRightMargin(0.15);
  h_R_pz_pos[i]->Draw("colz");
  can->cd(4); gPad->SetRightMargin(0.15);
  h_Z_pz_pos[i]->Draw("colz");
  can->cd(5); gPad->SetRightMargin(0.15);
  h_N_pz_pos[i]->Draw("colz");
  can->cd(6); gPad->SetRightMargin(0.15);
  h_P_pz_neg[i]->Draw("colz");
  can->cd(7); gPad->SetRightMargin(0.15);
  h_RP_pz_neg[i]->Draw("colz");
  can->cd(8); gPad->SetRightMargin(0.15);
  h_R_pz_neg[i]->Draw("colz");
  can->cd(9); gPad->SetRightMargin(0.15);
  h_Z_pz_neg[i]->Draw("colz");
  can->cd(10); gPad->SetRightMargin(0.15);
  h_N_pz_neg[i]->Draw("colz");
  
  gPad->RedrawAxis();
  
  can->Update();
  can->SaveAs(Form("figure/resid_phiz_from3D_%d_atR%d.pdf",runs[i],(int)selectR));
  delete can;

}


for (int i=0; i<nrun; i++)
{
  TFile* outfile = new TFile(Form("Rootfiles/Distortions_vsPZ_atR%d_%d.root",(int)selectR,runs[i]),"recreate");
  outfile->cd();
  h_N_pz_pos[i]->Write("hentries_posz");
  h_R_pz_pos[i]->Write("hIntDistortionR_posz");
  h_P_pz_pos[i]->Write("hIntDistortionP_posz");
  h_Z_pz_pos[i]->Write("hIntDistortionZ_posz");
  h_RP_pz_pos[i]->Write("hIntDistortionRP_posz");
  h_N_pz_neg[i]->Write("hentries_negz");
  h_R_pz_neg[i]->Write("hIntDistortionR_negz");
  h_P_pz_neg[i]->Write("hIntDistortionP_negz");
  h_Z_pz_neg[i]->Write("hIntDistortionZ_negz");
  h_RP_pz_neg[i]->Write("hIntDistortionRP_negz");
  outfile->Write();
  outfile->Close();
}
}

}
