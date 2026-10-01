#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

void plot1D_Ybin(TH3* h3, TH1* h1, float y)
{
//  TpcSpaceChargeReconstructionHelper::set_phi_range_central( {-1.73246,-1.43608} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_east( {-2.26272,-1.96089} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_west( {-1.21241,-0.909953} );

  //int xbin = h3->GetXaxis()->FindBin(((-1.73246)+(-1.43608))/2.+2*TMath::Pi());//central
  //int xbin = h3->GetXaxis()->FindBin(((-2.26272)+(-1.96089))/2.+2*TMath::Pi());//east
  int xbin = h3->GetXaxis()->FindBin(((-1.21241)+(-0.909953))/2.+2*TMath::Pi());//west
  int ybin = h3->GetYaxis()->FindBin(y);
  for (int i = 1; i <= h3->GetNbinsZ(); i++)
  {
      h1->SetBinContent(i, h3->GetBinContent(xbin, ybin, i));
      h1->SetBinError(i, h3->GetBinError(xbin, ybin, i));
  }
}

void draw1D_z_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::pair<double, double> thetarange_NCO = {0.609928,0.917104};
std::pair<double, double> thetarange_NCI = {0.0276558,0.566484};
std::pair<double, double> thetarange_SCI = {-0.568656,-0.0325273};
std::pair<double, double> thetarange_SCO = {-0.788212,-0.613404};

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
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun], *hcdb_RP_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun], *hcdb_RP_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun], *hcdb_RP_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun], *hcdb_RP_neg[nrun];

//get map from si-tpot fit matrix inversion
TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun], *h_RP_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun], *h_RP_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun], *h_RP_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun], *h_RP_neg[nrun];

//get map from si-tpot fit matrix inversion -- phi
TFile* file_3D_map_phi[nrun];
TH3 *h_N_prz_pos_phi[nrun], *h_R_prz_pos_phi[nrun], *h_P_prz_pos_phi[nrun], *h_Z_prz_pos_phi[nrun], *h_RP_prz_pos_phi[nrun];
TH3 *h_N_prz_neg_phi[nrun], *h_R_prz_neg_phi[nrun], *h_P_prz_neg_phi[nrun], *h_Z_prz_neg_phi[nrun], *h_RP_prz_neg_phi[nrun];
TH1 *h_N_pos_phi[nrun], *h_R_pos_phi[nrun], *h_P_pos_phi[nrun], *h_Z_pos_phi[nrun], *h_RP_pos_phi[nrun];
TH1 *h_N_neg_phi[nrun], *h_R_neg_phi[nrun], *h_P_neg_phi[nrun], *h_Z_neg_phi[nrun], *h_RP_neg_phi[nrun];

//get map from si-tpot fit matrix inversion -- z
TFile* file_3D_map_z[nrun];
TH3 *h_N_prz_pos_z[nrun], *h_R_prz_pos_z[nrun], *h_P_prz_pos_z[nrun], *h_Z_prz_pos_z[nrun], *h_RP_prz_pos_z[nrun];
TH3 *h_N_prz_neg_z[nrun], *h_R_prz_neg_z[nrun], *h_P_prz_neg_z[nrun], *h_Z_prz_neg_z[nrun], *h_RP_prz_neg_z[nrun];
TH1 *h_N_pos_z[nrun], *h_R_pos_z[nrun], *h_P_pos_z[nrun], *h_Z_pos_z[nrun], *h_RP_pos_z[nrun];
TH1 *h_N_neg_z[nrun], *h_R_neg_z[nrun], *h_P_neg_z[nrun], *h_Z_neg_z[nrun], *h_RP_neg_z[nrun];

//boundary of nco and nci
TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];
TLine *l_nco_o_RP[nrun], *l_nco_i_RP[nrun], *l_nci_o_RP[nrun], *l_nci_i_RP[nrun];

//boundary of sco and sci
TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];
TLine *l_sci_i_RP[nrun], *l_sci_o_RP[nrun], *l_sco_i_RP[nrun], *l_sco_o_RP[nrun];

double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
double ymax_RP=-100, ymin_RP=100;
for (int i = 0; i < nrun; i++)
{
  //get TH3 from matrix inversion root file
  file_3D_map[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_posz");
  //h_RP_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_posz");
  h_N_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hentries_posz");
  h_R_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_negz");
  //h_RP_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_negz");
  h_N_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hentries_negz");

  //get TH3 from matrix inversion root file -- phi
  file_3D_map_phi[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d_phi.root",runs[i]),"");
  h_R_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionZ_posz");
  //h_RP_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionRP_posz");
  h_N_prz_pos_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hentries_posz");
  h_R_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionZ_negz");
  //h_RP_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hIntDistortionRP_negz");
  h_N_prz_neg_phi[i] = (TH3*) file_3D_map_phi[i]->Get("hentries_negz");

  //get TH3 from matrix inversion root file -- z
  file_3D_map_z[i] = new TFile(Form("../../Rootfiles/Distortions_full_mm_%d_z.root",runs[i]),"");
  h_R_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionZ_posz");
  //h_RP_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionRP_posz");
  h_N_prz_pos_z[i] = (TH3*) file_3D_map_z[i]->Get("hentries_posz");
  h_R_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionZ_negz");
  //h_RP_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hIntDistortionRP_negz");
  h_N_prz_neg_z[i] = (TH3*) file_3D_map_z[i]->Get("hentries_negz");

  //define 1D hist
  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos[i]->GetZaxis()->GetNbins(),h_N_prz_pos[i]->GetZaxis()->GetXmin(),h_N_prz_pos[i]->GetZaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos[i]->GetZaxis()->GetNbins(),h_R_prz_pos[i]->GetZaxis()->GetXmin(),h_R_prz_pos[i]->GetZaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos[i]->GetZaxis()->GetNbins(),h_P_prz_pos[i]->GetZaxis()->GetXmin(),h_P_prz_pos[i]->GetZaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos[i]->GetZaxis()->GetNbins(),h_Z_prz_pos[i]->GetZaxis()->GetXmin(),h_Z_prz_pos[i]->GetZaxis()->GetXmax());
  //h_RP_pos[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_pos[i]->GetZaxis()->GetNbins(),h_RP_prz_pos[i]->GetZaxis()->GetXmin(),h_RP_prz_pos[i]->GetZaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg[i]->GetZaxis()->GetNbins(),h_N_prz_neg[i]->GetZaxis()->GetXmin(),h_N_prz_neg[i]->GetZaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg[i]->GetZaxis()->GetNbins(),h_R_prz_neg[i]->GetZaxis()->GetXmin(),h_R_prz_neg[i]->GetZaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg[i]->GetZaxis()->GetNbins(),h_P_prz_neg[i]->GetZaxis()->GetXmin(),h_P_prz_neg[i]->GetZaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg[i]->GetZaxis()->GetNbins(),h_Z_prz_neg[i]->GetZaxis()->GetXmin(),h_Z_prz_neg[i]->GetZaxis()->GetXmax());
  //h_RP_neg[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_neg[i]->GetZaxis()->GetNbins(),h_RP_prz_neg[i]->GetZaxis()->GetXmin(),h_RP_prz_neg[i]->GetZaxis()->GetXmax());

  //define 1D hist -- phi
  h_N_pos_phi[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos[i]->GetZaxis()->GetNbins(),h_N_prz_pos[i]->GetZaxis()->GetXmin(),h_N_prz_pos[i]->GetZaxis()->GetXmax());
  h_R_pos_phi[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos[i]->GetZaxis()->GetNbins(),h_R_prz_pos[i]->GetZaxis()->GetXmin(),h_R_prz_pos[i]->GetZaxis()->GetXmax());
  h_P_pos_phi[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos[i]->GetZaxis()->GetNbins(),h_P_prz_pos[i]->GetZaxis()->GetXmin(),h_P_prz_pos[i]->GetZaxis()->GetXmax());
  h_Z_pos_phi[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos[i]->GetZaxis()->GetNbins(),h_Z_prz_pos[i]->GetZaxis()->GetXmin(),h_Z_prz_pos[i]->GetZaxis()->GetXmax());
  //h_RP_pos_phi[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_pos[i]->GetZaxis()->GetNbins(),h_RP_prz_pos[i]->GetZaxis()->GetXmin(),h_RP_prz_pos[i]->GetZaxis()->GetXmax());
  h_N_neg_phi[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg[i]->GetZaxis()->GetNbins(),h_N_prz_neg[i]->GetZaxis()->GetXmin(),h_N_prz_neg[i]->GetZaxis()->GetXmax());
  h_R_neg_phi[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg[i]->GetZaxis()->GetNbins(),h_R_prz_neg[i]->GetZaxis()->GetXmin(),h_R_prz_neg[i]->GetZaxis()->GetXmax());
  h_P_neg_phi[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg[i]->GetZaxis()->GetNbins(),h_P_prz_neg[i]->GetZaxis()->GetXmin(),h_P_prz_neg[i]->GetZaxis()->GetXmax());
  h_Z_neg_phi[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg[i]->GetZaxis()->GetNbins(),h_Z_prz_neg[i]->GetZaxis()->GetXmin(),h_Z_prz_neg[i]->GetZaxis()->GetXmax());
  //h_RP_neg_phi[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_neg[i]->GetZaxis()->GetNbins(),h_RP_prz_neg[i]->GetZaxis()->GetXmin(),h_RP_prz_neg[i]->GetZaxis()->GetXmax());

  //define 1D hist -- z
  h_N_pos_z[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_pos[i]->GetZaxis()->GetNbins(),h_N_prz_pos[i]->GetZaxis()->GetXmin(),h_N_prz_pos[i]->GetZaxis()->GetXmax());
  h_R_pos_z[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_pos[i]->GetZaxis()->GetNbins(),h_R_prz_pos[i]->GetZaxis()->GetXmin(),h_R_prz_pos[i]->GetZaxis()->GetXmax());
  h_P_pos_z[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_pos[i]->GetZaxis()->GetNbins(),h_P_prz_pos[i]->GetZaxis()->GetXmin(),h_P_prz_pos[i]->GetZaxis()->GetXmax());
  h_Z_pos_z[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_pos[i]->GetZaxis()->GetNbins(),h_Z_prz_pos[i]->GetZaxis()->GetXmin(),h_Z_prz_pos[i]->GetZaxis()->GetXmax());
  //h_RP_pos_z[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_pos[i]->GetZaxis()->GetNbins(),h_RP_prz_pos[i]->GetZaxis()->GetXmin(),h_RP_prz_pos[i]->GetZaxis()->GetXmax());
  h_N_neg_z[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ R=%d cm;Z (cm);N",(int)selectR),h_N_prz_neg[i]->GetZaxis()->GetNbins(),h_N_prz_neg[i]->GetZaxis()->GetXmin(),h_N_prz_neg[i]->GetZaxis()->GetXmax());
  h_R_neg_z[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ R=%d cm;Z (cm);dR (cm)",(int)selectR),h_R_prz_neg[i]->GetZaxis()->GetNbins(),h_R_prz_neg[i]->GetZaxis()->GetXmin(),h_R_prz_neg[i]->GetZaxis()->GetXmax());
  h_P_neg_z[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ R=%d cm;Z (cm);dphi (rad)",(int)selectR),h_P_prz_neg[i]->GetZaxis()->GetNbins(),h_P_prz_neg[i]->GetZaxis()->GetXmin(),h_P_prz_neg[i]->GetZaxis()->GetXmax());
  h_Z_neg_z[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ R=%d cm;Z (cm);dz (cm)",(int)selectR),h_Z_prz_neg[i]->GetZaxis()->GetNbins(),h_Z_prz_neg[i]->GetZaxis()->GetXmin(),h_Z_prz_neg[i]->GetZaxis()->GetXmax());
  //h_RP_neg_z[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rdphi @ R=%d cm;Z (cm);Rdphi (cm)",(int)selectR),h_RP_prz_neg[i]->GetZaxis()->GetNbins(),h_RP_prz_neg[i]->GetZaxis()->GetXmin(),h_RP_prz_neg[i]->GetZaxis()->GetXmax());

  //plot 1D hist
  plot1D_Ybin(h_N_prz_pos[i],h_N_pos[i],selectR);
  plot1D_Ybin(h_R_prz_pos[i],h_R_pos[i],selectR);
  plot1D_Ybin(h_P_prz_pos[i],h_P_pos[i],selectR);
  plot1D_Ybin(h_Z_prz_pos[i],h_Z_pos[i],selectR);
  //plot1D_Ybin(h_RP_prz_pos[i],h_RP_pos[i],selectR);
  plot1D_Ybin(h_N_prz_neg[i],h_N_neg[i],selectR);
  plot1D_Ybin(h_R_prz_neg[i],h_R_neg[i],selectR);
  plot1D_Ybin(h_P_prz_neg[i],h_P_neg[i],selectR);
  plot1D_Ybin(h_Z_prz_neg[i],h_Z_neg[i],selectR);
  //plot1D_Ybin(h_RP_prz_neg[i],h_RP_neg[i],selectR);

  //plot 1D hist -- hist
  plot1D_Ybin(h_N_prz_pos_phi[i],h_N_pos_phi[i],selectR);
  plot1D_Ybin(h_R_prz_pos_phi[i],h_R_pos_phi[i],selectR);
  plot1D_Ybin(h_P_prz_pos_phi[i],h_P_pos_phi[i],selectR);
  plot1D_Ybin(h_Z_prz_pos_phi[i],h_Z_pos_phi[i],selectR);
  //plot1D_Ybin(h_RP_prz_pos_phi[i],h_RP_pos_phi[i],selectR);
  plot1D_Ybin(h_N_prz_neg_phi[i],h_N_neg_phi[i],selectR);
  plot1D_Ybin(h_R_prz_neg_phi[i],h_R_neg_phi[i],selectR);
  plot1D_Ybin(h_P_prz_neg_phi[i],h_P_neg_phi[i],selectR);
  plot1D_Ybin(h_Z_prz_neg_phi[i],h_Z_neg_phi[i],selectR);
  //plot1D_Ybin(h_RP_prz_neg_phi[i],h_RP_neg_phi[i],selectR);

  //plot 1D hist -- hist
  plot1D_Ybin(h_N_prz_pos_z[i],h_N_pos_z[i],selectR);
  plot1D_Ybin(h_R_prz_pos_z[i],h_R_pos_z[i],selectR);
  plot1D_Ybin(h_P_prz_pos_z[i],h_P_pos_z[i],selectR);
  plot1D_Ybin(h_Z_prz_pos_z[i],h_Z_pos_z[i],selectR);
  //plot1D_Ybin(h_RP_prz_pos_z[i],h_RP_pos_z[i],selectR);
  plot1D_Ybin(h_N_prz_neg_z[i],h_N_neg_z[i],selectR);
  plot1D_Ybin(h_R_prz_neg_z[i],h_R_neg_z[i],selectR);
  plot1D_Ybin(h_P_prz_neg_z[i],h_P_neg_z[i],selectR);
  plot1D_Ybin(h_Z_prz_neg_z[i],h_Z_neg_z[i],selectR);
  //plot1D_Ybin(h_RP_prz_neg_z[i],h_RP_neg_z[i],selectR);

  //set 1D hist style
  h_N_pos[i]->SetLineColor(i+2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(i+2);
  h_R_pos[i]->SetLineColor(i+2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(i+2);
  h_P_pos[i]->SetLineColor(i+2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(i+2);
  h_Z_pos[i]->SetLineColor(i+2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(i+2);
  //h_RP_pos[i]->SetLineColor(i+2); h_RP_pos[i]->SetLineWidth(1); h_RP_pos[i]->SetFillColor(0); h_RP_pos[i]->SetMarkerColor(i+2);
  h_N_neg[i]->SetLineColor(i+2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_R_neg[i]->SetLineColor(i+2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_P_neg[i]->SetLineColor(i+2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(i+2);
  h_Z_neg[i]->SetLineColor(i+2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(i+2);
  //h_RP_neg[i]->SetLineColor(i+2); h_RP_neg[i]->SetLineWidth(1); h_RP_neg[i]->SetFillColor(0); h_RP_neg[i]->SetMarkerColor(i+2);

  //set 1D hist style -- phi
  h_N_pos_phi[i]->SetLineColor(i+4); h_N_pos_phi[i]->SetLineWidth(1); h_N_pos_phi[i]->SetFillColor(0); h_N_pos_phi[i]->SetMarkerColor(i+4);
  h_R_pos_phi[i]->SetLineColor(i+4); h_R_pos_phi[i]->SetLineWidth(1); h_R_pos_phi[i]->SetFillColor(0); h_R_pos_phi[i]->SetMarkerColor(i+4);
  h_P_pos_phi[i]->SetLineColor(i+4); h_P_pos_phi[i]->SetLineWidth(1); h_P_pos_phi[i]->SetFillColor(0); h_P_pos_phi[i]->SetMarkerColor(i+4);
  h_Z_pos_phi[i]->SetLineColor(i+4); h_Z_pos_phi[i]->SetLineWidth(1); h_Z_pos_phi[i]->SetFillColor(0); h_Z_pos_phi[i]->SetMarkerColor(i+4);
  //h_RP_pos_phi[i]->SetLineColor(i+4); h_RP_pos_phi[i]->SetLineWidth(1); h_RP_pos_phi[i]->SetFillColor(0); h_RP_pos_phi[i]->SetMarkerColor(i+4);
  h_N_neg_phi[i]->SetLineColor(i+4); h_N_neg_phi[i]->SetLineWidth(1); h_N_neg_phi[i]->SetFillColor(0); h_R_neg_phi[i]->SetMarkerColor(i+4);
  h_R_neg_phi[i]->SetLineColor(i+4); h_R_neg_phi[i]->SetLineWidth(1); h_R_neg_phi[i]->SetFillColor(0); h_R_neg_phi[i]->SetMarkerColor(i+4);
  h_P_neg_phi[i]->SetLineColor(i+4); h_P_neg_phi[i]->SetLineWidth(1); h_P_neg_phi[i]->SetFillColor(0); h_P_neg_phi[i]->SetMarkerColor(i+4);
  h_Z_neg_phi[i]->SetLineColor(i+4); h_Z_neg_phi[i]->SetLineWidth(1); h_Z_neg_phi[i]->SetFillColor(0); h_Z_neg_phi[i]->SetMarkerColor(i+4);
  //h_RP_neg_phi[i]->SetLineColor(i+4); h_RP_neg_phi[i]->SetLineWidth(1); h_RP_neg_phi[i]->SetFillColor(0); h_RP_neg_phi[i]->SetMarkerColor(i+4);

  //set 1D hist style -- z
  h_N_pos_z[i]->SetLineColor(i+3); h_N_pos_z[i]->SetLineWidth(1); h_N_pos_z[i]->SetFillColor(0); h_N_pos_z[i]->SetMarkerColor(i+3);
  h_R_pos_z[i]->SetLineColor(i+3); h_R_pos_z[i]->SetLineWidth(1); h_R_pos_z[i]->SetFillColor(0); h_R_pos_z[i]->SetMarkerColor(i+3);
  h_P_pos_z[i]->SetLineColor(i+3); h_P_pos_z[i]->SetLineWidth(1); h_P_pos_z[i]->SetFillColor(0); h_P_pos_z[i]->SetMarkerColor(i+3);
  h_Z_pos_z[i]->SetLineColor(i+3); h_Z_pos_z[i]->SetLineWidth(1); h_Z_pos_z[i]->SetFillColor(0); h_Z_pos_z[i]->SetMarkerColor(i+3);
  //h_RP_pos_z[i]->SetLineColor(i+3); h_RP_pos_z[i]->SetLineWidth(1); h_RP_pos_z[i]->SetFillColor(0); h_RP_pos_z[i]->SetMarkerColor(i+3);
  h_N_neg_z[i]->SetLineColor(i+3); h_N_neg_z[i]->SetLineWidth(1); h_N_neg_z[i]->SetFillColor(0); h_R_neg_z[i]->SetMarkerColor(i+3);
  h_R_neg_z[i]->SetLineColor(i+3); h_R_neg_z[i]->SetLineWidth(1); h_R_neg_z[i]->SetFillColor(0); h_R_neg_z[i]->SetMarkerColor(i+3);
  h_P_neg_z[i]->SetLineColor(i+3); h_P_neg_z[i]->SetLineWidth(1); h_P_neg_z[i]->SetFillColor(0); h_P_neg_z[i]->SetMarkerColor(i+3);
  h_Z_neg_z[i]->SetLineColor(i+3); h_Z_neg_z[i]->SetLineWidth(1); h_Z_neg_z[i]->SetFillColor(0); h_Z_neg_z[i]->SetMarkerColor(i+3);
  //h_RP_neg_z[i]->SetLineColor(i+3); h_RP_neg_z[i]->SetLineWidth(1); h_RP_neg_z[i]->SetFillColor(0); h_RP_neg_z[i]->SetMarkerColor(i+3);

  //set ymax, ymin
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

  //if (h_RP_neg[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_neg[i]->GetMaximum();
  //if (h_RP_pos[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_pos[i]->GetMaximum();
  //if (h_RP_neg[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_neg[i]->GetMinimum();
  //if (h_RP_pos[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_pos[i]->GetMinimum();

  if (h_N_neg_phi[i]->GetMaximum()>ymax_N) ymax_N = h_N_neg_phi[i]->GetMaximum();
  if (h_N_pos_phi[i]->GetMaximum()>ymax_N) ymax_N = h_N_pos_phi[i]->GetMaximum();
  if (h_N_neg_phi[i]->GetMinimum()<ymin_N) ymin_N = h_N_neg_phi[i]->GetMinimum();
  if (h_N_pos_phi[i]->GetMinimum()<ymin_N) ymin_N = h_N_pos_phi[i]->GetMinimum();

  if (h_R_neg_phi[i]->GetMaximum()>ymax_R) ymax_R = h_R_neg_phi[i]->GetMaximum();
  if (h_R_pos_phi[i]->GetMaximum()>ymax_R) ymax_R = h_R_pos_phi[i]->GetMaximum();
  if (h_R_neg_phi[i]->GetMinimum()<ymin_R) ymin_R = h_R_neg_phi[i]->GetMinimum();
  if (h_R_pos_phi[i]->GetMinimum()<ymin_R) ymin_R = h_R_pos_phi[i]->GetMinimum();

  if (h_P_neg_phi[i]->GetMaximum()>ymax_P) ymax_P = h_P_neg_phi[i]->GetMaximum();
  if (h_P_pos_phi[i]->GetMaximum()>ymax_P) ymax_P = h_P_pos_phi[i]->GetMaximum();
  if (h_P_neg_phi[i]->GetMinimum()<ymin_P) ymin_P = h_P_neg_phi[i]->GetMinimum();
  if (h_P_pos_phi[i]->GetMinimum()<ymin_P) ymin_P = h_P_pos_phi[i]->GetMinimum();

  if (h_Z_neg_phi[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_neg_phi[i]->GetMaximum();
  if (h_Z_pos_phi[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_pos_phi[i]->GetMaximum();
  if (h_Z_neg_phi[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_neg_phi[i]->GetMinimum();
  if (h_Z_pos_phi[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_pos_phi[i]->GetMinimum();

  //if (h_RP_neg_phi[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_neg_phi[i]->GetMaximum();
  //if (h_RP_pos_phi[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_pos_phi[i]->GetMaximum();
  //if (h_RP_neg_phi[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_neg_phi[i]->GetMinimum();
  //if (h_RP_pos_phi[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_pos_phi[i]->GetMinimum();

  if (h_N_neg_z[i]->GetMaximum()>ymax_N) ymax_N = h_N_neg_z[i]->GetMaximum();
  if (h_N_pos_z[i]->GetMaximum()>ymax_N) ymax_N = h_N_pos_z[i]->GetMaximum();
  if (h_N_neg_z[i]->GetMinimum()<ymin_N) ymin_N = h_N_neg_z[i]->GetMinimum();
  if (h_N_pos_z[i]->GetMinimum()<ymin_N) ymin_N = h_N_pos_z[i]->GetMinimum();

  if (h_R_neg_z[i]->GetMaximum()>ymax_R) ymax_R = h_R_neg_z[i]->GetMaximum();
  if (h_R_pos_z[i]->GetMaximum()>ymax_R) ymax_R = h_R_pos_z[i]->GetMaximum();
  if (h_R_neg_z[i]->GetMinimum()<ymin_R) ymin_R = h_R_neg_z[i]->GetMinimum();
  if (h_R_pos_z[i]->GetMinimum()<ymin_R) ymin_R = h_R_pos_z[i]->GetMinimum();

  if (h_P_neg_z[i]->GetMaximum()>ymax_P) ymax_P = h_P_neg_z[i]->GetMaximum();
  if (h_P_pos_z[i]->GetMaximum()>ymax_P) ymax_P = h_P_pos_z[i]->GetMaximum();
  if (h_P_neg_z[i]->GetMinimum()<ymin_P) ymin_P = h_P_neg_z[i]->GetMinimum();
  if (h_P_pos_z[i]->GetMinimum()<ymin_P) ymin_P = h_P_pos_z[i]->GetMinimum();

  if (h_Z_neg_z[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_neg_z[i]->GetMaximum();
  if (h_Z_pos_z[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_pos_z[i]->GetMaximum();
  if (h_Z_neg_z[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_neg_z[i]->GetMinimum();
  if (h_Z_pos_z[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_pos_z[i]->GetMinimum();

  //if (h_RP_neg_z[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_neg_z[i]->GetMaximum();
  //if (h_RP_pos_z[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_pos_z[i]->GetMaximum();
  //if (h_RP_neg_z[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_neg_z[i]->GetMinimum();
  //if (h_RP_pos_z[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_pos_z[i]->GetMinimum();

}

if (ymax_N>0) ymax_N *= 1.1; else ymax_N = 1;
if (ymax_R>0) ymax_R *= 1.1; else ymax_R = 0;
if (ymax_P>0) ymax_P *= 1.1; else ymax_P = 0;
if (ymax_Z>0) ymax_Z *= 1.1; else ymax_Z = 0;
//if (ymax_RP>0) ymax_RP *= 1.1; else ymax_RP = 0;

if (ymin_N<0) ymin_N *= 1.1; else ymin_N = 0.9;
if (ymin_R<0) ymin_R *= 1.1; else ymin_R = 0.0;
if (ymin_P<0) ymin_P *= 1.1; else ymin_P = 0.0;
if (ymin_Z<0) ymin_Z *= 1.1; else ymin_Z = 0.0;
//if (ymin_RP<0) ymin_RP *= 1.1; else ymin_RP = 0.0;

ymax_R = 2.5;
ymin_R = -2.5;
ymax_Z = 2;
ymin_Z = -2;

for (int i = 0; i < nrun; i++)
{
  // set ymin, ymax
  h_N_pos[i]->SetMinimum(ymin_N); h_N_pos[i]->SetMaximum(ymax_N);
  h_R_pos[i]->SetMinimum(ymin_R); h_R_pos[i]->SetMaximum(ymax_R);
  h_P_pos[i]->SetMinimum(ymin_P); h_P_pos[i]->SetMaximum(ymax_P);
  h_Z_pos[i]->SetMinimum(ymin_Z); h_Z_pos[i]->SetMaximum(ymax_Z);
  //h_RP_pos[i]->SetMinimum(ymin_RP); h_RP_pos[i]->SetMaximum(ymax_RP);

  h_N_neg[i]->SetMinimum(ymin_N); h_N_neg[i]->SetMaximum(ymax_N);
  h_R_neg[i]->SetMinimum(ymin_R); h_R_neg[i]->SetMaximum(ymax_R);
  h_P_neg[i]->SetMinimum(ymin_P); h_P_neg[i]->SetMaximum(ymax_P);
  h_Z_neg[i]->SetMinimum(ymin_Z); h_Z_neg[i]->SetMaximum(ymax_Z);
  //h_RP_neg[i]->SetMinimum(ymin_RP); h_RP_neg[i]->SetMaximum(ymax_RP);

  // set ymin, ymax -- phi
  h_N_pos_phi[i]->SetMinimum(ymin_N); h_N_pos_phi[i]->SetMaximum(ymax_N);
  h_R_pos_phi[i]->SetMinimum(ymin_R); h_R_pos_phi[i]->SetMaximum(ymax_R);
  h_P_pos_phi[i]->SetMinimum(ymin_P); h_P_pos_phi[i]->SetMaximum(ymax_P);
  h_Z_pos_phi[i]->SetMinimum(ymin_Z); h_Z_pos_phi[i]->SetMaximum(ymax_Z);
  //h_RP_pos_phi[i]->SetMinimum(ymin_RP); h_RP_pos_phi[i]->SetMaximum(ymax_RP);

  h_N_neg_phi[i]->SetMinimum(ymin_N); h_N_neg_phi[i]->SetMaximum(ymax_N);
  h_R_neg_phi[i]->SetMinimum(ymin_R); h_R_neg_phi[i]->SetMaximum(ymax_R);
  h_P_neg_phi[i]->SetMinimum(ymin_P); h_P_neg_phi[i]->SetMaximum(ymax_P);
  h_Z_neg_phi[i]->SetMinimum(ymin_Z); h_Z_neg_phi[i]->SetMaximum(ymax_Z);
  //h_RP_neg_phi[i]->SetMinimum(ymin_RP); h_RP_neg_phi[i]->SetMaximum(ymax_RP);

  // set ymin, ymax -- z
  h_N_pos_z[i]->SetMinimum(ymin_N); h_N_pos_z[i]->SetMaximum(ymax_N);
  h_R_pos_z[i]->SetMinimum(ymin_R); h_R_pos_z[i]->SetMaximum(ymax_R);
  h_P_pos_z[i]->SetMinimum(ymin_P); h_P_pos_z[i]->SetMaximum(ymax_P);
  h_Z_pos_z[i]->SetMinimum(ymin_Z); h_Z_pos_z[i]->SetMaximum(ymax_Z);
  //h_RP_pos_z[i]->SetMinimum(ymin_RP); h_RP_pos_z[i]->SetMaximum(ymax_RP);

  h_N_neg_z[i]->SetMinimum(ymin_N); h_N_neg_z[i]->SetMaximum(ymax_N);
  h_R_neg_z[i]->SetMinimum(ymin_R); h_R_neg_z[i]->SetMaximum(ymax_R);
  h_P_neg_z[i]->SetMinimum(ymin_P); h_P_neg_z[i]->SetMaximum(ymax_P);
  h_Z_neg_z[i]->SetMinimum(ymin_Z); h_Z_neg_z[i]->SetMaximum(ymax_Z);
  //h_RP_neg_z[i]->SetMinimum(ymin_RP); h_RP_neg_z[i]->SetMaximum(ymax_RP);

  l_nco_o_R[i] = new TLine(selectR*tan(thetarange_NCO.second),ymin_R,selectR*tan(thetarange_NCO.second),ymax_R);
  l_nco_o_P[i] = new TLine(selectR*tan(thetarange_NCO.second),ymin_P,selectR*tan(thetarange_NCO.second),ymax_P);
  l_nco_o_Z[i] = new TLine(selectR*tan(thetarange_NCO.second),ymin_Z,selectR*tan(thetarange_NCO.second),ymax_Z);
  //l_nco_o_RP[i] = new TLine(selectR*tan(thetarange_NCO.second),ymin_RP,selectR*tan(thetarange_NCO.second),ymax_RP);

  l_nco_i_R[i] = new TLine(selectR*tan(thetarange_NCO.first),ymin_R,selectR*tan(thetarange_NCO.first),ymax_R);
  l_nco_i_P[i] = new TLine(selectR*tan(thetarange_NCO.first),ymin_P,selectR*tan(thetarange_NCO.first),ymax_P);
  l_nco_i_Z[i] = new TLine(selectR*tan(thetarange_NCO.first),ymin_Z,selectR*tan(thetarange_NCO.first),ymax_Z);
  //l_nco_i_RP[i] = new TLine(selectR*tan(thetarange_NCO.first),ymin_RP,selectR*tan(thetarange_NCO.first),ymax_RP);

  l_nci_o_R[i] = new TLine(selectR*tan(thetarange_NCI.second),ymin_R,selectR*tan(thetarange_NCI.second),ymax_R);
  l_nci_o_P[i] = new TLine(selectR*tan(thetarange_NCI.second),ymin_P,selectR*tan(thetarange_NCI.second),ymax_P);
  l_nci_o_Z[i] = new TLine(selectR*tan(thetarange_NCI.second),ymin_Z,selectR*tan(thetarange_NCI.second),ymax_Z);
  //l_nci_o_RP[i] = new TLine(selectR*tan(thetarange_NCI.second),ymin_RP,selectR*tan(thetarange_NCI.second),ymax_RP);

  l_nci_i_R[i] = new TLine(selectR*tan(thetarange_NCI.first),ymin_R,selectR*tan(thetarange_NCI.first),ymax_R);
  l_nci_i_P[i] = new TLine(selectR*tan(thetarange_NCI.first),ymin_P,selectR*tan(thetarange_NCI.first),ymax_P);
  l_nci_i_Z[i] = new TLine(selectR*tan(thetarange_NCI.first),ymin_Z,selectR*tan(thetarange_NCI.first),ymax_Z);
  //l_nci_i_RP[i] = new TLine(selectR*tan(thetarange_NCI.first),ymin_RP,selectR*tan(thetarange_NCI.first),ymax_RP);

  l_sci_i_R[i] = new TLine(selectR*tan(thetarange_SCI.second),ymin_R,selectR*tan(thetarange_SCI.second),ymax_R);
  l_sci_i_P[i] = new TLine(selectR*tan(thetarange_SCI.second),ymin_P,selectR*tan(thetarange_SCI.second),ymax_P);
  l_sci_i_Z[i] = new TLine(selectR*tan(thetarange_SCI.second),ymin_Z,selectR*tan(thetarange_SCI.second),ymax_Z);
  //l_sci_i_RP[i] = new TLine(selectR*tan(thetarange_SCI.second),ymin_RP,selectR*tan(thetarange_SCI.second),ymax_RP);

  l_sci_o_R[i] = new TLine(selectR*tan(thetarange_SCI.first),ymin_R,selectR*tan(thetarange_SCI.first),ymax_R);
  l_sci_o_P[i] = new TLine(selectR*tan(thetarange_SCI.first),ymin_P,selectR*tan(thetarange_SCI.first),ymax_P);
  l_sci_o_Z[i] = new TLine(selectR*tan(thetarange_SCI.first),ymin_Z,selectR*tan(thetarange_SCI.first),ymax_Z);
  //l_sci_o_RP[i] = new TLine(selectR*tan(thetarange_SCI.first),ymin_RP,selectR*tan(thetarange_SCI.first),ymax_RP);

  l_sco_i_R[i] = new TLine(selectR*tan(thetarange_SCO.second),ymin_R,selectR*tan(thetarange_SCO.second),ymax_R);
  l_sco_i_P[i] = new TLine(selectR*tan(thetarange_SCO.second),ymin_P,selectR*tan(thetarange_SCO.second),ymax_P);
  l_sco_i_Z[i] = new TLine(selectR*tan(thetarange_SCO.second),ymin_Z,selectR*tan(thetarange_SCO.second),ymax_Z);
  //l_sco_i_RP[i] = new TLine(selectR*tan(thetarange_SCO.second),ymin_RP,selectR*tan(thetarange_SCO.second),ymax_RP);

  l_sco_o_R[i] = new TLine(selectR*tan(thetarange_SCO.first),ymin_R,selectR*tan(thetarange_SCO.first),ymax_R);
  l_sco_o_P[i] = new TLine(selectR*tan(thetarange_SCO.first),ymin_P,selectR*tan(thetarange_SCO.first),ymax_P);
  l_sco_o_Z[i] = new TLine(selectR*tan(thetarange_SCO.first),ymin_Z,selectR*tan(thetarange_SCO.first),ymax_Z);
  //l_sco_o_RP[i] = new TLine(selectR*tan(thetarange_SCO.first),ymin_RP,selectR*tan(thetarange_SCO.first),ymax_RP);

}

//write 1D hist to root file
TFile* ofile_R = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"recreate");
ofile_R->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos[i]->Write();
  h_R_neg[i]->Write();
}
ofile_R->Write();

TFile* ofile_P = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"recreate");
ofile_P->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos[i]->Write();
  h_P_neg[i]->Write();
}
ofile_P->Write();

TFile* ofile_Z = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"recreate");
ofile_Z->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos[i]->Write();
  h_Z_neg[i]->Write();
}
ofile_Z->Write();

//write 1D hist to root file -- phi
TFile* ofile_R_phi = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_phi.root",(int)selectR),"recreate");
ofile_R_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_phi[i]->Write();
  h_R_neg_phi[i]->Write();
}
ofile_R_phi->Write();

TFile* ofile_P_phi = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d_phi.root",(int)selectR),"recreate");
ofile_P_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_phi[i]->Write();
  h_P_neg_phi[i]->Write();
}
ofile_P_phi->Write();

TFile* ofile_Z_phi = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d_phi.root",(int)selectR),"recreate");
ofile_Z_phi->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_phi[i]->Write();
  h_Z_neg_phi[i]->Write();
}
ofile_Z_phi->Write();

//write 1D hist to root file -- z
TFile* ofile_R_z = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_z.root",(int)selectR),"recreate");
ofile_R_z->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_z[i]->Write();
  h_R_neg_z[i]->Write();
}
ofile_R_z->Write();

TFile* ofile_P_z = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d_z.root",(int)selectR),"recreate");
ofile_P_z->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_z[i]->Write();
  h_P_neg_z[i]->Write();
}
ofile_P_z->Write();

TFile* ofile_Z_z = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d_z.root",(int)selectR),"recreate");
ofile_Z_z->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_z[i]->Write();
  h_Z_neg_z[i]->Write();
}
ofile_Z_z->Write();

TLine *l_zero_pos = new TLine(h_N_prz_pos[0]->GetZaxis()->GetXmin(),0,h_N_prz_pos[0]->GetZaxis()->GetXmax(),0);
TLine *l_zero_neg = new TLine(h_N_prz_neg[0]->GetZaxis()->GetXmin(),0,h_N_prz_neg[0]->GetZaxis()->GetXmax(),0);

l_zero_pos->SetLineStyle(2); l_zero_pos->SetLineColor(kBlack); l_zero_pos->SetLineWidth(1);
l_zero_neg->SetLineStyle(2); l_zero_neg->SetLineColor(kBlack); l_zero_neg->SetLineWidth(1);

//plot
TCanvas* can = new TCanvas("can","",4000,1200);
can->Divide(5,2);
can->cd(1);
//for (int i=0; i<nrun; i++) h_RP_pos[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_RP_pos_phi[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_RP_pos_z[i]->Draw("hist,e,same");
//l_zero_pos->Draw();
for (int i=0; i<nrun; i++)
{
//  l_nco_i_RP[i]->Draw();
//  l_nco_o_RP[i]->Draw();
//  l_nci_i_RP[i]->Draw();
//  l_nci_o_RP[i]->Draw();
}
can->cd(2);
for (int i=0; i<nrun; i++) h_P_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_P_pos_phi[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_P_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
for (int i=0; i<nrun; i++)
{
  l_nco_i_P[i]->Draw();
  l_nco_o_P[i]->Draw();
  l_nci_i_P[i]->Draw();
  l_nci_o_P[i]->Draw();
}
can->cd(3);
for (int i=0; i<nrun; i++) h_R_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_R_pos_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_R_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
for (int i=0; i<nrun; i++)
{
  l_nco_i_R[i]->Draw();
  l_nco_o_R[i]->Draw();
  l_nci_i_R[i]->Draw();
  l_nci_o_R[i]->Draw();
}
can->cd(4);
for (int i=0; i<nrun; i++) h_Z_pos[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_Z_pos_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_Z_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
for (int i=0; i<nrun; i++)
{
  l_nco_i_Z[i]->Draw();
  l_nco_o_Z[i]->Draw();
  l_nci_i_Z[i]->Draw();
  l_nci_o_Z[i]->Draw();
}
can->cd(5);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_N_pos_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_N_pos_z[i]->Draw("hist,e,same");
l_zero_pos->Draw();
can->cd(6);
//for (int i=0; i<nrun; i++) h_RP_neg[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_RP_neg_phi[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_RP_neg_z[i]->Draw("hist,e,same");
//l_zero_neg->Draw();
for (int i=0; i<nrun; i++)
{
//  l_sco_i_RP[i]->Draw();
//  l_sco_o_RP[i]->Draw();
//  l_sci_i_RP[i]->Draw();
//  l_sci_o_RP[i]->Draw();
}
can->cd(7);
for (int i=0; i<nrun; i++) h_P_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_P_neg_phi[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_P_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();
for (int i=0; i<nrun; i++)
{
  l_sco_i_P[i]->Draw();
  l_sco_o_P[i]->Draw();
  l_sci_i_P[i]->Draw();
  l_sci_o_P[i]->Draw();
}
can->cd(8);
for (int i=0; i<nrun; i++) h_R_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_R_neg_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_R_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();
for (int i=0; i<nrun; i++)
{
  l_sco_i_R[i]->Draw();
  l_sco_o_R[i]->Draw();
  l_sci_i_R[i]->Draw();
  l_sci_o_R[i]->Draw();
}
can->cd(9);
for (int i=0; i<nrun; i++) h_Z_neg[i]->Draw("hist,e,same");
//for (int i=0; i<nrun; i++) h_Z_neg_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_Z_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();
for (int i=0; i<nrun; i++)
{
  l_sco_i_Z[i]->Draw();
  l_sco_o_Z[i]->Draw();
  l_sci_i_Z[i]->Draw();
  l_sci_o_Z[i]->Draw();
}
can->cd(10);
gPad->SetLogy(1);
for (int i=0; i<nrun; i++) h_N_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_N_neg_phi[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) h_N_neg_z[i]->Draw("hist,e,same");
l_zero_neg->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d.pdf",(int)selectR));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos[i], Form("Run %d, MBD rate %d kHz",runs[i],mbdrates[i]), "l");
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos_phi[i], Form("Run %d, MBD rate %d kHz (phi)",runs[i],mbdrates[i]), "l");
for (int i=0; i<nrun; i++) legend->AddEntry(h_P_pos_z[i], Form("Run %d, MBD rate %d kHz (z)",runs[i],mbdrates[i]), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsZ_from3D_leg.pdf"));

delete can;
delete can_leg;
}

}
