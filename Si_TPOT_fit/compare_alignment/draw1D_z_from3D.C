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

const int nrun = 2;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nrun] = {79507,79516};

TFile* file_3D_map_method1[nrun];
TH3 *h_N_prz_pos_method1[nrun], *h_R_prz_pos_method1[nrun], *h_P_prz_pos_method1[nrun], *h_Z_prz_pos_method1[nrun];
TH3 *h_N_prz_neg_method1[nrun], *h_R_prz_neg_method1[nrun], *h_P_prz_neg_method1[nrun], *h_Z_prz_neg_method1[nrun];
TH1 *h_N_pos_method1[nrun], *h_R_pos_method1[nrun], *h_P_pos_method1[nrun], *h_Z_pos_method1[nrun];
TH1 *h_N_neg_method1[nrun], *h_R_neg_method1[nrun], *h_P_neg_method1[nrun], *h_Z_neg_method1[nrun];

TFile* file_3D_map_method2[nrun];
TH3 *h_N_prz_pos_method2[nrun], *h_R_prz_pos_method2[nrun], *h_P_prz_pos_method2[nrun], *h_Z_prz_pos_method2[nrun];
TH3 *h_N_prz_neg_method2[nrun], *h_R_prz_neg_method2[nrun], *h_P_prz_neg_method2[nrun], *h_Z_prz_neg_method2[nrun];
TH1 *h_N_pos_method2[nrun], *h_R_pos_method2[nrun], *h_P_pos_method2[nrun], *h_Z_pos_method2[nrun];
TH1 *h_N_neg_method2[nrun], *h_R_neg_method2[nrun], *h_P_neg_method2[nrun], *h_Z_neg_method2[nrun];

TFile* file_3D_map_method3[nrun];
TH3 *h_N_prz_pos_method3[nrun], *h_R_prz_pos_method3[nrun], *h_P_prz_pos_method3[nrun], *h_Z_prz_pos_method3[nrun];
TH3 *h_N_prz_neg_method3[nrun], *h_R_prz_neg_method3[nrun], *h_P_prz_neg_method3[nrun], *h_Z_prz_neg_method3[nrun];
TH1 *h_N_pos_method3[nrun], *h_R_pos_method3[nrun], *h_P_pos_method3[nrun], *h_Z_pos_method3[nrun];
TH1 *h_N_neg_method3[nrun], *h_R_neg_method3[nrun], *h_P_neg_method3[nrun], *h_Z_neg_method3[nrun];

TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];
TLine *l_nco_o_RP[nrun], *l_nco_i_RP[nrun], *l_nci_o_RP[nrun], *l_nci_i_RP[nrun];
TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];
TLine *l_sci_i_RP[nrun], *l_sci_o_RP[nrun], *l_sco_i_RP[nrun], *l_sco_o_RP[nrun];

for (int i = 0; i < nrun; i++)
{
  file_3D_map_method1[i] = new TFile(Form("../run3pp_newAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hentries_posz");
  h_R_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hentries_negz");

  h_N_pos_method1[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos_method1[i]->GetZaxis()->GetNbins(),h_N_prz_pos_method1[i]->GetZaxis()->GetXmin(),h_N_prz_pos_method1[i]->GetZaxis()->GetXmax());
  h_R_pos_method1[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos_method1[i]->GetZaxis()->GetNbins(),h_R_prz_pos_method1[i]->GetZaxis()->GetXmin(),h_R_prz_pos_method1[i]->GetZaxis()->GetXmax());
  h_P_pos_method1[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos_method1[i]->GetZaxis()->GetNbins(),h_P_prz_pos_method1[i]->GetZaxis()->GetXmin(),h_P_prz_pos_method1[i]->GetZaxis()->GetXmax());
  h_Z_pos_method1[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos_method1[i]->GetZaxis()->GetNbins(),h_Z_prz_pos_method1[i]->GetZaxis()->GetXmin(),h_Z_prz_pos_method1[i]->GetZaxis()->GetXmax());
  h_N_neg_method1[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg_method1[i]->GetZaxis()->GetNbins(),h_N_prz_neg_method1[i]->GetZaxis()->GetXmin(),h_N_prz_neg_method1[i]->GetZaxis()->GetXmax());
  h_R_neg_method1[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg_method1[i]->GetZaxis()->GetNbins(),h_R_prz_neg_method1[i]->GetZaxis()->GetXmin(),h_R_prz_neg_method1[i]->GetZaxis()->GetXmax());
  h_P_neg_method1[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg_method1[i]->GetZaxis()->GetNbins(),h_P_prz_neg_method1[i]->GetZaxis()->GetXmin(),h_P_prz_neg_method1[i]->GetZaxis()->GetXmax());
  h_Z_neg_method1[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg_method1[i]->GetZaxis()->GetNbins(),h_Z_prz_neg_method1[i]->GetZaxis()->GetXmin(),h_Z_prz_neg_method1[i]->GetZaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method1[i],h_N_pos_method1[i],selectR,false);
  plot1D_Zbin(h_R_prz_pos_method1[i],h_R_pos_method1[i],selectR,false);
  plot1D_Zbin(h_P_prz_pos_method1[i],h_P_pos_method1[i],selectR,false);
  plot1D_Zbin(h_Z_prz_pos_method1[i],h_Z_pos_method1[i],selectR,false);
  plot1D_Zbin(h_N_prz_neg_method1[i],h_N_neg_method1[i],selectR,false);
  plot1D_Zbin(h_R_prz_neg_method1[i],h_R_neg_method1[i],selectR,false);
  plot1D_Zbin(h_P_prz_neg_method1[i],h_P_neg_method1[i],selectR,false);
  plot1D_Zbin(h_Z_prz_neg_method1[i],h_Z_neg_method1[i],selectR,false);
  h_N_pos_method1[i]->SetLineColor(2); h_N_pos_method1[i]->SetLineWidth(1); h_N_pos_method1[i]->SetFillColor(0); h_N_pos_method1[i]->SetMarkerColor(2);
  h_R_pos_method1[i]->SetLineColor(2); h_R_pos_method1[i]->SetLineWidth(1); h_R_pos_method1[i]->SetFillColor(0); h_R_pos_method1[i]->SetMarkerColor(2);
  h_P_pos_method1[i]->SetLineColor(2); h_P_pos_method1[i]->SetLineWidth(1); h_P_pos_method1[i]->SetFillColor(0); h_P_pos_method1[i]->SetMarkerColor(2);
  h_Z_pos_method1[i]->SetLineColor(2); h_Z_pos_method1[i]->SetLineWidth(1); h_Z_pos_method1[i]->SetFillColor(0); h_Z_pos_method1[i]->SetMarkerColor(2);
  h_N_neg_method1[i]->SetLineColor(2); h_N_neg_method1[i]->SetLineWidth(1); h_N_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_R_neg_method1[i]->SetLineColor(2); h_R_neg_method1[i]->SetLineWidth(1); h_R_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_P_neg_method1[i]->SetLineColor(2); h_P_neg_method1[i]->SetLineWidth(1); h_P_neg_method1[i]->SetFillColor(0); h_P_neg_method1[i]->SetMarkerColor(2);
  h_Z_neg_method1[i]->SetLineColor(2); h_Z_neg_method1[i]->SetLineWidth(1); h_Z_neg_method1[i]->SetFillColor(0); h_Z_neg_method1[i]->SetMarkerColor(2);

  file_3D_map_method2[i] = new TFile(Form("../run3pp_oldAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_posz");
  h_R_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_negz");

  h_N_pos_method2[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos_method2[i]->GetZaxis()->GetNbins(),h_N_prz_pos_method2[i]->GetZaxis()->GetXmin(),h_N_prz_pos_method2[i]->GetZaxis()->GetXmax());
  h_R_pos_method2[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos_method2[i]->GetZaxis()->GetNbins(),h_R_prz_pos_method2[i]->GetZaxis()->GetXmin(),h_R_prz_pos_method2[i]->GetZaxis()->GetXmax());
  h_P_pos_method2[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos_method2[i]->GetZaxis()->GetNbins(),h_P_prz_pos_method2[i]->GetZaxis()->GetXmin(),h_P_prz_pos_method2[i]->GetZaxis()->GetXmax());
  h_Z_pos_method2[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos_method2[i]->GetZaxis()->GetNbins(),h_Z_prz_pos_method2[i]->GetZaxis()->GetXmin(),h_Z_prz_pos_method2[i]->GetZaxis()->GetXmax());
  h_N_neg_method2[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg_method2[i]->GetZaxis()->GetNbins(),h_N_prz_neg_method2[i]->GetZaxis()->GetXmin(),h_N_prz_neg_method2[i]->GetZaxis()->GetXmax());
  h_R_neg_method2[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg_method2[i]->GetZaxis()->GetNbins(),h_R_prz_neg_method2[i]->GetZaxis()->GetXmin(),h_R_prz_neg_method2[i]->GetZaxis()->GetXmax());
  h_P_neg_method2[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg_method2[i]->GetZaxis()->GetNbins(),h_P_prz_neg_method2[i]->GetZaxis()->GetXmin(),h_P_prz_neg_method2[i]->GetZaxis()->GetXmax());
  h_Z_neg_method2[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg_method2[i]->GetZaxis()->GetNbins(),h_Z_prz_neg_method2[i]->GetZaxis()->GetXmin(),h_Z_prz_neg_method2[i]->GetZaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method2[i],h_N_pos_method2[i],selectR,false);
  plot1D_Zbin(h_R_prz_pos_method2[i],h_R_pos_method2[i],selectR,false);
  plot1D_Zbin(h_P_prz_pos_method2[i],h_P_pos_method2[i],selectR,false);
  plot1D_Zbin(h_Z_prz_pos_method2[i],h_Z_pos_method2[i],selectR,false);
  plot1D_Zbin(h_N_prz_neg_method2[i],h_N_neg_method2[i],selectR,false);
  plot1D_Zbin(h_R_prz_neg_method2[i],h_R_neg_method2[i],selectR,false);
  plot1D_Zbin(h_P_prz_neg_method2[i],h_P_neg_method2[i],selectR,false);
  plot1D_Zbin(h_Z_prz_neg_method2[i],h_Z_neg_method2[i],selectR,false);
  h_N_pos_method2[i]->SetLineColor(4); h_N_pos_method2[i]->SetLineWidth(1); h_N_pos_method2[i]->SetFillColor(0); h_N_pos_method2[i]->SetMarkerColor(4);
  h_R_pos_method2[i]->SetLineColor(4); h_R_pos_method2[i]->SetLineWidth(1); h_R_pos_method2[i]->SetFillColor(0); h_R_pos_method2[i]->SetMarkerColor(4);
  h_P_pos_method2[i]->SetLineColor(4); h_P_pos_method2[i]->SetLineWidth(1); h_P_pos_method2[i]->SetFillColor(0); h_P_pos_method2[i]->SetMarkerColor(4);
  h_Z_pos_method2[i]->SetLineColor(4); h_Z_pos_method2[i]->SetLineWidth(1); h_Z_pos_method2[i]->SetFillColor(0); h_Z_pos_method2[i]->SetMarkerColor(4);
  h_N_neg_method2[i]->SetLineColor(4); h_N_neg_method2[i]->SetLineWidth(1); h_N_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor(4);
  h_R_neg_method2[i]->SetLineColor(4); h_R_neg_method2[i]->SetLineWidth(1); h_R_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor(4);
  h_P_neg_method2[i]->SetLineColor(4); h_P_neg_method2[i]->SetLineWidth(1); h_P_neg_method2[i]->SetFillColor(0); h_P_neg_method2[i]->SetMarkerColor(4);
  h_Z_neg_method2[i]->SetLineColor(4); h_Z_neg_method2[i]->SetLineWidth(1); h_Z_neg_method2[i]->SetFillColor(0); h_Z_neg_method2[i]->SetMarkerColor(4);

  file_3D_map_method3[i] = new TFile(Form("../run3pp_newSiFieldonAlignment_newTPOTzfAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hentries_posz");
  h_R_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hentries_negz");

  h_N_pos_method3[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos_method3[i]->GetZaxis()->GetNbins(),h_N_prz_pos_method3[i]->GetZaxis()->GetXmin(),h_N_prz_pos_method3[i]->GetZaxis()->GetXmax());
  h_R_pos_method3[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos_method3[i]->GetZaxis()->GetNbins(),h_R_prz_pos_method3[i]->GetZaxis()->GetXmin(),h_R_prz_pos_method3[i]->GetZaxis()->GetXmax());
  h_P_pos_method3[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos_method3[i]->GetZaxis()->GetNbins(),h_P_prz_pos_method3[i]->GetZaxis()->GetXmin(),h_P_prz_pos_method3[i]->GetZaxis()->GetXmax());
  h_Z_pos_method3[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos_method3[i]->GetZaxis()->GetNbins(),h_Z_prz_pos_method3[i]->GetZaxis()->GetXmin(),h_Z_prz_pos_method3[i]->GetZaxis()->GetXmax());
  h_N_neg_method3[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg_method3[i]->GetZaxis()->GetNbins(),h_N_prz_neg_method3[i]->GetZaxis()->GetXmin(),h_N_prz_neg_method3[i]->GetZaxis()->GetXmax());
  h_R_neg_method3[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg_method3[i]->GetZaxis()->GetNbins(),h_R_prz_neg_method3[i]->GetZaxis()->GetXmin(),h_R_prz_neg_method3[i]->GetZaxis()->GetXmax());
  h_P_neg_method3[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg_method3[i]->GetZaxis()->GetNbins(),h_P_prz_neg_method3[i]->GetZaxis()->GetXmin(),h_P_prz_neg_method3[i]->GetZaxis()->GetXmax());
  h_Z_neg_method3[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg_method3[i]->GetZaxis()->GetNbins(),h_Z_prz_neg_method3[i]->GetZaxis()->GetXmin(),h_Z_prz_neg_method3[i]->GetZaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method3[i],h_N_pos_method3[i],selectR,false);
  plot1D_Zbin(h_R_prz_pos_method3[i],h_R_pos_method3[i],selectR,false);
  plot1D_Zbin(h_P_prz_pos_method3[i],h_P_pos_method3[i],selectR,false);
  plot1D_Zbin(h_Z_prz_pos_method3[i],h_Z_pos_method3[i],selectR,false);
  plot1D_Zbin(h_N_prz_neg_method3[i],h_N_neg_method3[i],selectR,false);
  plot1D_Zbin(h_R_prz_neg_method3[i],h_R_neg_method3[i],selectR,false);
  plot1D_Zbin(h_P_prz_neg_method3[i],h_P_neg_method3[i],selectR,false);
  plot1D_Zbin(h_Z_prz_neg_method3[i],h_Z_neg_method3[i],selectR,false);
  h_N_pos_method3[i]->SetLineColor(8); h_N_pos_method3[i]->SetLineWidth(1); h_N_pos_method3[i]->SetFillColor(0); h_N_pos_method3[i]->SetMarkerColor(8);
  h_R_pos_method3[i]->SetLineColor(8); h_R_pos_method3[i]->SetLineWidth(1); h_R_pos_method3[i]->SetFillColor(0); h_R_pos_method3[i]->SetMarkerColor(8);
  h_P_pos_method3[i]->SetLineColor(8); h_P_pos_method3[i]->SetLineWidth(1); h_P_pos_method3[i]->SetFillColor(0); h_P_pos_method3[i]->SetMarkerColor(8);
  h_Z_pos_method3[i]->SetLineColor(8); h_Z_pos_method3[i]->SetLineWidth(1); h_Z_pos_method3[i]->SetFillColor(0); h_Z_pos_method3[i]->SetMarkerColor(8);
  h_N_neg_method3[i]->SetLineColor(8); h_N_neg_method3[i]->SetLineWidth(1); h_N_neg_method3[i]->SetFillColor(0); h_R_neg_method3[i]->SetMarkerColor(8);
  h_R_neg_method3[i]->SetLineColor(8); h_R_neg_method3[i]->SetLineWidth(1); h_R_neg_method3[i]->SetFillColor(0); h_R_neg_method3[i]->SetMarkerColor(8);
  h_P_neg_method3[i]->SetLineColor(8); h_P_neg_method3[i]->SetLineWidth(1); h_P_neg_method3[i]->SetFillColor(0); h_P_neg_method3[i]->SetMarkerColor(8);
  h_Z_neg_method3[i]->SetLineColor(8); h_Z_neg_method3[i]->SetLineWidth(1); h_Z_neg_method3[i]->SetFillColor(0); h_Z_neg_method3[i]->SetMarkerColor(8);

  /*
  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_method1[i], h_R_pos_method1[i], h_R_neg_method2[i], h_R_pos_method2[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_method1[i], h_P_pos_method1[i], h_P_neg_method2[i], h_P_pos_method2[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method2[i], h_Z_pos_method2[i]});
  */
}

std::vector<TH1*> hists_P; hists_P.clear();
std::vector<TH1*> hists_R; hists_R.clear();
std::vector<TH1*> hists_Z; hists_Z.clear();
for (int i=0; i<nrun; i++)
{
	hists_P.push_back(h_P_neg_method1[i]);
	hists_P.push_back(h_P_pos_method1[i]);
	hists_P.push_back(h_P_neg_method2[i]);
	hists_P.push_back(h_P_pos_method2[i]);
	hists_P.push_back(h_P_neg_method3[i]);
	hists_P.push_back(h_P_pos_method3[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_R.push_back(h_R_neg_method1[i]);
	hists_R.push_back(h_R_pos_method1[i]);
	hists_R.push_back(h_R_neg_method2[i]);
	hists_R.push_back(h_R_pos_method2[i]);
	hists_R.push_back(h_R_neg_method3[i]);
	hists_R.push_back(h_R_pos_method3[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_Z.push_back(h_Z_neg_method1[i]);
	hists_Z.push_back(h_Z_pos_method1[i]);
	hists_Z.push_back(h_Z_neg_method2[i]);
	hists_Z.push_back(h_Z_pos_method2[i]);
	hists_Z.push_back(h_Z_neg_method3[i]);
	hists_Z.push_back(h_Z_pos_method3[i]);
}
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

for (int i=0; i<nrun; i++)
{
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
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_pos_method3[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
  h_R_neg_method3[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_pos_method3[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
  h_P_neg_method3[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_pos_method3[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
  h_Z_neg_method3[i]->Write();
}
ofile3->Write();

for (int i = 0; i < nrun; i++)
{
TCanvas* can = new TCanvas("can","",2400,1200);
can->Divide(3,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->Draw("hist,e,same");
h_P_pos_method2[i]->Draw("hist,e,same");
h_P_pos_method3[i]->Draw("hist,e,same");
l_nco_i_P[i]->Draw();
l_nco_o_P[i]->Draw();
l_nci_i_P[i]->Draw();
l_nci_o_P[i]->Draw();
can->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->Draw("hist,e,same");
h_R_pos_method2[i]->Draw("hist,e,same");
h_R_pos_method3[i]->Draw("hist,e,same");
l_nco_i_R[i]->Draw();
l_nco_o_R[i]->Draw();
l_nci_i_R[i]->Draw();
l_nci_o_R[i]->Draw();
can->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
h_Z_pos_method2[i]->Draw("hist,e,same");
h_Z_pos_method3[i]->Draw("hist,e,same");
l_nco_i_Z[i]->Draw();
l_nco_o_Z[i]->Draw();
l_nci_i_Z[i]->Draw();
l_nci_o_Z[i]->Draw();
can->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->Draw("hist,e,same");
h_P_neg_method2[i]->Draw("hist,e,same");
h_P_neg_method3[i]->Draw("hist,e,same");
l_sco_i_P[i]->Draw();
l_sco_o_P[i]->Draw();
l_sci_i_P[i]->Draw();
l_sci_o_P[i]->Draw();
can->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->Draw("hist,e,same");
h_R_neg_method2[i]->Draw("hist,e,same");
h_R_neg_method3[i]->Draw("hist,e,same");
l_sco_i_R[i]->Draw();
l_sco_o_R[i]->Draw();
l_sci_i_R[i]->Draw();
l_sci_o_R[i]->Draw();
can->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
h_Z_neg_method2[i]->Draw("hist,e,same");
h_Z_neg_method3[i]->Draw("hist,e,same");
l_sco_i_Z[i]->Draw();
l_sco_o_Z[i]->Draw();
l_sci_i_Z[i]->Draw();
l_sci_o_Z[i]->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d_%d.pdf",(int)selectR,runs[i]));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d",runs[i]));
legend->AddEntry(h_P_pos_method2[i], Form("TPOT-based correction, old alignment"), "l");
legend->AddEntry(h_P_pos_method1[i], Form("TPOT-based correction, new Si field-on & TPOT old alignment"), "l");
legend->AddEntry(h_P_pos_method3[i], Form("TPOT-based correction, new Si field-on & TPOT new ZF alignment"), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsZ_from3D_leg_%d.pdf",runs[i]));

delete can;
delete can_leg;
}

}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
	//recover default max and min
	h->SetMinimum();
	h->SetMaximum();
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin > 0) { yMin *= 0.8; }
    else if (yMin < 0) { yMin *= 1.2; }

    if (yMax > 0) { yMax *= 1.2; }
    else if (yMax < 0) { yMax *= 0.8; }

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
