#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void plot1D_Zbin(TH3* h3, TH1* h1, float y, bool convert_P_2_RP)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int ybin = h3->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h3->GetNbinsZ(); i++)
  {
      double value = h3->GetBinContent(xbin, ybin, i);
      double error = h3->GetBinError(xbin, ybin, i);
      if (convert_P_2_RP)
      {
        value *= y;
        error *= y;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void draw1D_z_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::pair<double, double> phirange_NCO = {0.609928,0.917104};
std::pair<double, double> phirange_NCI = {0.0276558,0.566484};
std::pair<double, double> phirange_SCI = {-0.568656,-0.0325273};
std::pair<double, double> phirange_SCO = {-0.788212,-0.613404};

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
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun];

TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun];
TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];
TLine *l_nco_o_RP[nrun], *l_nco_i_RP[nrun], *l_nci_o_RP[nrun], *l_nci_i_RP[nrun];
TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];
TLine *l_sci_i_RP[nrun], *l_sci_o_RP[nrun], *l_sco_i_RP[nrun], *l_sco_o_RP[nrun];
double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
for (int i = 0; i < nrun; i++)
{
  file_3D_map[i] = new TFile(Form("Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hentries_posz");
  h_R_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hentries_negz");

  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos[i]->GetZaxis()->GetNbins(),h_N_prz_pos[i]->GetZaxis()->GetXmin(),h_N_prz_pos[i]->GetZaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos[i]->GetZaxis()->GetNbins(),h_R_prz_pos[i]->GetZaxis()->GetXmin(),h_R_prz_pos[i]->GetZaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (rad)",(int)selectR),h_P_prz_pos[i]->GetZaxis()->GetNbins(),h_P_prz_pos[i]->GetZaxis()->GetXmin(),h_P_prz_pos[i]->GetZaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos[i]->GetZaxis()->GetNbins(),h_Z_prz_pos[i]->GetZaxis()->GetXmin(),h_Z_prz_pos[i]->GetZaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg[i]->GetZaxis()->GetNbins(),h_N_prz_neg[i]->GetZaxis()->GetXmin(),h_N_prz_neg[i]->GetZaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg[i]->GetZaxis()->GetNbins(),h_R_prz_neg[i]->GetZaxis()->GetXmin(),h_R_prz_neg[i]->GetZaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (rad)",(int)selectR),h_P_prz_neg[i]->GetZaxis()->GetNbins(),h_P_prz_neg[i]->GetZaxis()->GetXmin(),h_P_prz_neg[i]->GetZaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg[i]->GetZaxis()->GetNbins(),h_Z_prz_neg[i]->GetZaxis()->GetXmin(),h_Z_prz_neg[i]->GetZaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos[i],h_N_pos[i],selectR,false);
  plot1D_Zbin(h_R_prz_pos[i],h_R_pos[i],selectR,false);
  plot1D_Zbin(h_P_prz_pos[i],h_P_pos[i],selectR,true);
  plot1D_Zbin(h_Z_prz_pos[i],h_Z_pos[i],selectR,false);
  plot1D_Zbin(h_N_prz_neg[i],h_N_neg[i],selectR,false);
  plot1D_Zbin(h_R_prz_neg[i],h_R_neg[i],selectR,false);
  plot1D_Zbin(h_P_prz_neg[i],h_P_neg[i],selectR,true);
  plot1D_Zbin(h_Z_prz_neg[i],h_Z_neg[i],selectR,false);
  h_N_pos[i]->SetLineColor(i+2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(i+2);
  h_R_pos[i]->SetLineColor(i+2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(i+2);
  h_P_pos[i]->SetLineColor(i+2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(i+2);
  h_Z_pos[i]->SetLineColor(i+2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(i+2);
  h_N_neg[i]->SetLineColor(i+2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_R_neg[i]->SetLineColor(i+2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_P_neg[i]->SetLineColor(i+2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(i+2);
  h_Z_neg[i]->SetLineColor(i+2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(i+2);

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg[i], h_R_pos[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg[i], h_P_pos[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg[i], h_Z_pos[i]});

  l_nco_o_R[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_R.first,selectR*tan(phirange_NCO.second),yrange_R.second);
  l_nco_o_P[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_P.first,selectR*tan(phirange_NCO.second),yrange_P.second);
  l_nco_o_Z[i] = new TLine(selectR*tan(phirange_NCO.second),yrange_Z.first,selectR*tan(phirange_NCO.second),yrange_Z.second);

  l_nco_i_R[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_R.first,selectR*tan(phirange_NCO.first),yrange_R.second);
  l_nco_i_P[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_P.first,selectR*tan(phirange_NCO.first),yrange_P.second);
  l_nco_i_Z[i] = new TLine(selectR*tan(phirange_NCO.first),yrange_Z.first,selectR*tan(phirange_NCO.first),yrange_Z.second);

  l_nci_o_R[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_R.first,selectR*tan(phirange_NCI.second),yrange_R.second);
  l_nci_o_P[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_P.first,selectR*tan(phirange_NCI.second),yrange_P.second);
  l_nci_o_Z[i] = new TLine(selectR*tan(phirange_NCI.second),yrange_Z.first,selectR*tan(phirange_NCI.second),yrange_Z.second);

  l_nci_i_R[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_R.first,selectR*tan(phirange_NCI.first),yrange_R.second);
  l_nci_i_P[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_P.first,selectR*tan(phirange_NCI.first),yrange_P.second);
  l_nci_i_Z[i] = new TLine(selectR*tan(phirange_NCI.first),yrange_Z.first,selectR*tan(phirange_NCI.first),yrange_Z.second);

  l_sci_i_R[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_R.first,selectR*tan(phirange_SCI.second),yrange_R.second);
  l_sci_i_P[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_P.first,selectR*tan(phirange_SCI.second),yrange_P.second);
  l_sci_i_Z[i] = new TLine(selectR*tan(phirange_SCI.second),yrange_Z.first,selectR*tan(phirange_SCI.second),yrange_Z.second);

  l_sci_o_R[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_R.first,selectR*tan(phirange_SCI.first),yrange_R.second);
  l_sci_o_P[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_P.first,selectR*tan(phirange_SCI.first),yrange_P.second);
  l_sci_o_Z[i] = new TLine(selectR*tan(phirange_SCI.first),yrange_Z.first,selectR*tan(phirange_SCI.first),yrange_Z.second);

  l_sco_i_R[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_R.first,selectR*tan(phirange_SCO.second),yrange_R.second);
  l_sco_i_P[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_P.first,selectR*tan(phirange_SCO.second),yrange_P.second);
  l_sco_i_Z[i] = new TLine(selectR*tan(phirange_SCO.second),yrange_Z.first,selectR*tan(phirange_SCO.second),yrange_Z.second);

  l_sco_o_R[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_R.first,selectR*tan(phirange_SCO.first),yrange_R.second);
  l_sco_o_P[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_P.first,selectR*tan(phirange_SCO.first),yrange_P.second);
  l_sco_o_Z[i] = new TLine(selectR*tan(phirange_SCO.first),yrange_Z.first,selectR*tan(phirange_SCO.first),yrange_Z.second);

}

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos[i]->Write();
  h_R_neg[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_rdphi_1Dz_R%d.root",(int)selectR),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos[i]->Write();
  h_P_neg[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos[i]->Write();
  h_Z_neg[i]->Write();
}
ofile3->Write();

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_P_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_nco_i_P[i]->Draw();
  l_nco_o_P[i]->Draw();
  l_nci_i_P[i]->Draw();
  l_nci_o_P[i]->Draw();
}
can->cd(2);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_R_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_nco_i_R[i]->Draw();
  l_nco_o_R[i]->Draw();
  l_nci_i_R[i]->Draw();
  l_nci_o_R[i]->Draw();
}
can->cd(3);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_Z_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_nco_i_Z[i]->Draw();
  l_nco_o_Z[i]->Draw();
  l_nci_i_Z[i]->Draw();
  l_nci_o_Z[i]->Draw();
}
can->cd(4);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_pos[i]->Draw("hist,e,same");
can->cd(5);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_P_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_sco_i_P[i]->Draw();
  l_sco_o_P[i]->Draw();
  l_sci_i_P[i]->Draw();
  l_sci_o_P[i]->Draw();
}
can->cd(6);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_R_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_sco_i_R[i]->Draw();
  l_sco_o_R[i]->Draw();
  l_sci_i_R[i]->Draw();
  l_sci_o_R[i]->Draw();
}
can->cd(7);
gPad->SetLogy(0);
for (int i=0; i<nrun; i++) h_Z_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++)
{
  l_sco_i_Z[i]->Draw();
  l_sco_o_Z[i]->Draw();
  l_sci_i_Z[i]->Draw();
  l_sci_o_Z[i]->Draw();
}
can->cd(8);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_neg[i]->Draw("hist,e,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d.pdf",(int)selectR));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos[i], Form("Run %d, MBD rate %d kHz",runs[i],mbdrates[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsZ_from3D_leg.pdf"));

/*
for (int i=0; i<nrun; i++)
{
  TFile* outfile = new TFile(Form("Rootfiles/Distortions_vsR_atZ%d_%d.root",(int)selectR,runs[i]),"recreate");
  outfile->cd();
  h_N_pos[i]->Write("hentries_posz");
  h_R_pos[i]->Write("hIntDistortionR_posz");
  h_P_pos[i]->Write("hIntDistortionP_posz");
  h_Z_pos[i]->Write("hIntDistortionZ_posz");
  h_N_neg[i]->Write("hentries_negz");
  h_R_neg[i]->Write("hIntDistortionR_negz");
  h_P_neg[i]->Write("hIntDistortionP_negz");
  h_Z_neg[i]->Write("hIntDistortionZ_negz");
  outfile->Write();
  outfile->Close();
}
*/

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

    if (yMin<-20) yMin=-20;
    if (yMax>20) yMax=20;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
