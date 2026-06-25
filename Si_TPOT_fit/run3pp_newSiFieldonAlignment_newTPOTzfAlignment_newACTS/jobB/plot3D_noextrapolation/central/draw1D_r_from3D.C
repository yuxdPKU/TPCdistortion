#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>
#include "utilities.h"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

TString phiType = "central";

void draw1D_r_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::vector<float> selectZs={2, 5, 10, 15, 20, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

const int nrun = 2;
int runs[nrun] = {79507,79516};

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun];

//get map from si-tpot fit matrix inversion
TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun];

//get map from si-tpot fit matrix inversion -- phi
TFile* file_3D_map_phi[nrun];
TH3 *h_N_prz_pos_phi[nrun], *h_R_prz_pos_phi[nrun], *h_P_prz_pos_phi[nrun], *h_Z_prz_pos_phi[nrun];
TH3 *h_N_prz_neg_phi[nrun], *h_R_prz_neg_phi[nrun], *h_P_prz_neg_phi[nrun], *h_Z_prz_neg_phi[nrun];
TH1 *h_N_pos_phi[nrun], *h_R_pos_phi[nrun], *h_P_pos_phi[nrun], *h_Z_pos_phi[nrun];
TH1 *h_N_neg_phi[nrun], *h_R_neg_phi[nrun], *h_P_neg_phi[nrun], *h_Z_neg_phi[nrun];

//get map from si-tpot fit matrix inversion -- z
TFile* file_3D_map_z[nrun];
TH3 *h_N_prz_pos_z[nrun], *h_R_prz_pos_z[nrun], *h_P_prz_pos_z[nrun], *h_Z_prz_pos_z[nrun];
TH3 *h_N_prz_neg_z[nrun], *h_R_prz_neg_z[nrun], *h_P_prz_neg_z[nrun], *h_Z_prz_neg_z[nrun];
TH1 *h_N_pos_z[nrun], *h_R_pos_z[nrun], *h_P_pos_z[nrun], *h_Z_pos_z[nrun];
TH1 *h_N_neg_z[nrun], *h_R_neg_z[nrun], *h_P_neg_z[nrun], *h_Z_neg_z[nrun];

for (int i = 0; i < nrun; i++)
{
  //get TH3 from matrix inversion root file
  file_3D_map[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hentries_posz");
  h_R_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hentries_negz");

  //get TH3 from matrix inversion root file -- phi
  file_3D_map_phi[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d_phi.root",runs[i]),"");
  h_R_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hentries_posz");
  h_R_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hentries_negz");

  //get TH3 from matrix inversion root file -- z
  file_3D_map_z[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d_z.root",runs[i]),"");
  h_R_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hentries_posz");
  h_R_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hentries_negz");

  //define 1D hist
  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());

  //define 1D hist -- phi
  h_N_pos_phi[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos_phi[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos_phi[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos_phi[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg_phi[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg_phi[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg_phi[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg_phi[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());

  //define 1D hist -- z
  h_N_pos_z[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos_z[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos_z[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos_z[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg_z[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg_z[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg_z[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg_z[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());

  //plot 1D hist
  plot1D_Zbin(h_N_prz_pos[i],h_N_pos[i],selectZ, phiType);
  plot1D_Zbin(h_R_prz_pos[i],h_R_pos[i],selectZ, phiType);
  plot1D_Zbin(h_P_prz_pos[i],h_P_pos[i],selectZ, phiType);
  plot1D_Zbin(h_Z_prz_pos[i],h_Z_pos[i],selectZ, phiType);
  plot1D_Zbin(h_N_prz_neg[i],h_N_neg[i],-selectZ, phiType);
  plot1D_Zbin(h_R_prz_neg[i],h_R_neg[i],-selectZ, phiType);
  plot1D_Zbin(h_P_prz_neg[i],h_P_neg[i],-selectZ, phiType);
  plot1D_Zbin(h_Z_prz_neg[i],h_Z_neg[i],-selectZ, phiType);

  //plot 1D hist -- hist
  plot1D_Zbin(h_N_prz_pos_phi[i],h_N_pos_phi[i],selectZ, phiType);
  plot1D_Zbin(h_R_prz_pos_phi[i],h_R_pos_phi[i],selectZ, phiType);
  plot1D_Zbin(h_P_prz_pos_phi[i],h_P_pos_phi[i],selectZ, phiType);
  plot1D_Zbin(h_Z_prz_pos_phi[i],h_Z_pos_phi[i],selectZ, phiType);
  plot1D_Zbin(h_N_prz_neg_phi[i],h_N_neg_phi[i],-selectZ, phiType);
  plot1D_Zbin(h_R_prz_neg_phi[i],h_R_neg_phi[i],-selectZ, phiType);
  plot1D_Zbin(h_P_prz_neg_phi[i],h_P_neg_phi[i],-selectZ, phiType);
  plot1D_Zbin(h_Z_prz_neg_phi[i],h_Z_neg_phi[i],-selectZ, phiType);

  //plot 1D hist -- hist
  plot1D_Zbin(h_N_prz_pos_z[i],h_N_pos_z[i],selectZ, phiType);
  plot1D_Zbin(h_R_prz_pos_z[i],h_R_pos_z[i],selectZ, phiType);
  plot1D_Zbin(h_P_prz_pos_z[i],h_P_pos_z[i],selectZ, phiType);
  plot1D_Zbin(h_Z_prz_pos_z[i],h_Z_pos_z[i],selectZ, phiType);
  plot1D_Zbin(h_N_prz_neg_z[i],h_N_neg_z[i],-selectZ, phiType);
  plot1D_Zbin(h_R_prz_neg_z[i],h_R_neg_z[i],-selectZ, phiType);
  plot1D_Zbin(h_P_prz_neg_z[i],h_P_neg_z[i],-selectZ, phiType);
  plot1D_Zbin(h_Z_prz_neg_z[i],h_Z_neg_z[i],-selectZ, phiType);

  //set 1D hist style
  h_N_pos[i]->SetLineColor(2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(2);
  h_R_pos[i]->SetLineColor(2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(2);
  h_P_pos[i]->SetLineColor(2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(2);
  h_Z_pos[i]->SetLineColor(2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(2);
  h_N_neg[i]->SetLineColor(2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(2);
  h_R_neg[i]->SetLineColor(2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(2);
  h_P_neg[i]->SetLineColor(2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(2);
  h_Z_neg[i]->SetLineColor(2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(2);

  //set 1D hist style -- phi
  h_N_pos_phi[i]->SetLineColor(4); h_N_pos_phi[i]->SetLineWidth(1); h_N_pos_phi[i]->SetFillColor(0); h_N_pos_phi[i]->SetMarkerColor(4);
  h_R_pos_phi[i]->SetLineColor(4); h_R_pos_phi[i]->SetLineWidth(1); h_R_pos_phi[i]->SetFillColor(0); h_R_pos_phi[i]->SetMarkerColor(4);
  h_P_pos_phi[i]->SetLineColor(4); h_P_pos_phi[i]->SetLineWidth(1); h_P_pos_phi[i]->SetFillColor(0); h_P_pos_phi[i]->SetMarkerColor(4);
  h_Z_pos_phi[i]->SetLineColor(4); h_Z_pos_phi[i]->SetLineWidth(1); h_Z_pos_phi[i]->SetFillColor(0); h_Z_pos_phi[i]->SetMarkerColor(4);
  h_N_neg_phi[i]->SetLineColor(4); h_N_neg_phi[i]->SetLineWidth(1); h_N_neg_phi[i]->SetFillColor(0); h_R_neg_phi[i]->SetMarkerColor(4);
  h_R_neg_phi[i]->SetLineColor(4); h_R_neg_phi[i]->SetLineWidth(1); h_R_neg_phi[i]->SetFillColor(0); h_R_neg_phi[i]->SetMarkerColor(4);
  h_P_neg_phi[i]->SetLineColor(4); h_P_neg_phi[i]->SetLineWidth(1); h_P_neg_phi[i]->SetFillColor(0); h_P_neg_phi[i]->SetMarkerColor(4);
  h_Z_neg_phi[i]->SetLineColor(4); h_Z_neg_phi[i]->SetLineWidth(1); h_Z_neg_phi[i]->SetFillColor(0); h_Z_neg_phi[i]->SetMarkerColor(4);

  //set 1D hist style -- z
  h_N_pos_z[i]->SetLineColor(3); h_N_pos_z[i]->SetLineWidth(1); h_N_pos_z[i]->SetFillColor(0); h_N_pos_z[i]->SetMarkerColor(3);
  h_R_pos_z[i]->SetLineColor(3); h_R_pos_z[i]->SetLineWidth(1); h_R_pos_z[i]->SetFillColor(0); h_R_pos_z[i]->SetMarkerColor(3);
  h_P_pos_z[i]->SetLineColor(3); h_P_pos_z[i]->SetLineWidth(1); h_P_pos_z[i]->SetFillColor(0); h_P_pos_z[i]->SetMarkerColor(3);
  h_Z_pos_z[i]->SetLineColor(3); h_Z_pos_z[i]->SetLineWidth(1); h_Z_pos_z[i]->SetFillColor(0); h_Z_pos_z[i]->SetMarkerColor(3);
  h_N_neg_z[i]->SetLineColor(3); h_N_neg_z[i]->SetLineWidth(1); h_N_neg_z[i]->SetFillColor(0); h_R_neg_z[i]->SetMarkerColor(3);
  h_R_neg_z[i]->SetLineColor(3); h_R_neg_z[i]->SetLineWidth(1); h_R_neg_z[i]->SetFillColor(0); h_R_neg_z[i]->SetMarkerColor(3);
  h_P_neg_z[i]->SetLineColor(3); h_P_neg_z[i]->SetLineWidth(1); h_P_neg_z[i]->SetFillColor(0); h_P_neg_z[i]->SetMarkerColor(3);
  h_Z_neg_z[i]->SetLineColor(3); h_Z_neg_z[i]->SetLineWidth(1); h_Z_neg_z[i]->SetFillColor(0); h_Z_neg_z[i]->SetMarkerColor(3);

  std::pair<double,double> yrange_N = SetCommonYRange({h_N_neg[i],h_N_pos[i],h_N_neg_phi[i],h_N_pos_phi[i],h_N_neg_z[i],h_N_pos_z[i]});
  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg[i],h_R_pos[i],h_R_neg_phi[i],h_R_pos_phi[i],h_R_neg_z[i],h_R_pos_z[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg[i],h_P_pos[i],h_P_neg_phi[i],h_P_pos_phi[i],h_P_neg_z[i],h_P_pos_z[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg[i],h_Z_pos[i],h_Z_neg_phi[i],h_Z_pos_phi[i],h_Z_neg_z[i],h_Z_pos_z[i]});

  TLine *l_zero_pos = new TLine(h_N_prz_pos[0]->GetYaxis()->GetXmin(),0,h_N_prz_pos[0]->GetYaxis()->GetXmax(),0);
  TLine *l_zero_neg = new TLine(h_N_prz_neg[0]->GetYaxis()->GetXmin(),0,h_N_prz_neg[0]->GetYaxis()->GetXmax(),0);

  l_zero_pos->SetLineStyle(2); l_zero_pos->SetLineColor(kBlack); l_zero_pos->SetLineWidth(1);
  l_zero_neg->SetLineStyle(2); l_zero_neg->SetLineColor(kBlack); l_zero_neg->SetLineWidth(1);

//plot
TCanvas* can = new TCanvas("can","",2400,1200);
can->Divide(3,2);
can->cd(1);
h_P_pos[i]->Draw("hist,e,same");
h_P_pos_phi[i]->Draw("hist,e,same");
l_zero_pos->Draw();
can->cd(2);
h_R_pos[i]->Draw("hist,e,same");
h_R_pos_phi[i]->Draw("hist,e,same");
h_R_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
can->cd(3);
h_Z_pos[i]->Draw("hist,e,same");
h_Z_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
can->cd(4);
h_P_neg[i]->Draw("hist,e,same");
h_P_neg_phi[i]->Draw("hist,e,same");
l_zero_neg->Draw();
can->cd(5);
h_R_neg[i]->Draw("hist,e,same");
h_R_neg_phi[i]->Draw("hist,e,same");
h_R_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();
can->cd(6);
h_Z_neg[i]->Draw("hist,e,same");
h_Z_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_run%d.pdf",(int)selectZ,runs[i]));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(phiType);
legend->AddEntry(h_P_pos[i], Form("Run %d",runs[i]), "l");
legend->AddEntry(h_P_pos_phi[i], Form("Run %d (phi)",runs[i]), "l");
legend->AddEntry(h_P_pos_z[i], Form("Run %d (z)",runs[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsR_from3D_run%d_leg.pdf",runs[i]));

delete can;
delete can_leg;

}

//write 1D hist to root file
TFile* ofile_R = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile_R->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos[i]->Write();
  h_R_neg[i]->Write();
}
ofile_R->Write();

TFile* ofile_P = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile_P->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos[i]->Write();
  h_P_neg[i]->Write();
}
ofile_P->Write();

TFile* ofile_Z = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile_Z->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos[i]->Write();
  h_Z_neg[i]->Write();
}
ofile_Z->Write();

//write 1D hist to root file -- phi
TFile* ofile_R_phi = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d_phi.root",(int)selectZ),"recreate");
ofile_R_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_phi[i]->Write();
  h_R_neg_phi[i]->Write();
}
ofile_R_phi->Write();

TFile* ofile_P_phi = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d_phi.root",(int)selectZ),"recreate");
ofile_P_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_phi[i]->Write();
  h_P_neg_phi[i]->Write();
}
ofile_P_phi->Write();

TFile* ofile_Z_phi = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d_phi.root",(int)selectZ),"recreate");
ofile_Z_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_phi[i]->Write();
  h_Z_neg_phi[i]->Write();
}
ofile_Z_phi->Write();

//write 1D hist to root file -- z
TFile* ofile_R_z = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d_z.root",(int)selectZ),"recreate");
ofile_R_z->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_z[i]->Write();
  h_R_neg_z[i]->Write();
}
ofile_R_z->Write();

TFile* ofile_P_z = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d_z.root",(int)selectZ),"recreate");
ofile_P_z->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_z[i]->Write();
  h_P_neg_z[i]->Write();
}
ofile_P_z->Write();

TFile* ofile_Z_z = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d_z.root",(int)selectZ),"recreate");
ofile_Z_z->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_z[i]->Write();
  h_Z_neg_z[i]->Write();
}
ofile_Z_z->Write();

}

}
