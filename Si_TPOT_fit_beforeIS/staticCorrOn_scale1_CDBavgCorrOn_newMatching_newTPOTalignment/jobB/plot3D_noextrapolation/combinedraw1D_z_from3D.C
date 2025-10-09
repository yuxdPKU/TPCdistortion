#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

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

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}

void combinedraw1D_z_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::pair<double, double> phirange_NCO = {0.609928,0.917104};
std::pair<double, double> phirange_NCI = {0.0276558,0.566484};
std::pair<double, double> phirange_SCI = {-0.568656,-0.0325273};
//std::pair<double, double> phirange_SCO = {-0.788212,-0.613404};
std::pair<double, double> phirange_SCO = {-0.92,-0.613404};

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
{

const int nrun = 1;
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

TH1 *h_R_central_pos[nrun], *h_P_central_pos[nrun], *h_Z_central_pos[nrun], *h_RP_central_pos[nrun];
TH1 *h_R_west_pos[nrun], *h_P_west_pos[nrun], *h_Z_west_pos[nrun], *h_RP_west_pos[nrun];
TH1 *h_R_east_pos[nrun], *h_P_east_pos[nrun], *h_Z_east_pos[nrun], *h_RP_east_pos[nrun];
TH1 *hfit_R_rphi_central_pos[nrun], *hfit_R_zr_central_pos[nrun], *hfit_P_central_pos[nrun], *hfit_Z_central_pos[nrun];
TH1 *hfit_R_rphi_west_pos[nrun], *hfit_R_zr_west_pos[nrun], *hfit_P_west_pos[nrun], *hfit_Z_west_pos[nrun];
TH1 *hfit_R_rphi_east_pos[nrun], *hfit_R_zr_east_pos[nrun], *hfit_P_east_pos[nrun], *hfit_Z_east_pos[nrun];
TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];
TH1 *h_R_central_neg[nrun], *h_P_central_neg[nrun], *h_Z_central_neg[nrun], *h_RP_central_neg[nrun];
TH1 *h_R_west_neg[nrun], *h_P_west_neg[nrun], *h_Z_west_neg[nrun], *h_RP_west_neg[nrun];
TH1 *h_R_east_neg[nrun], *h_P_east_neg[nrun], *h_Z_east_neg[nrun], *h_RP_east_neg[nrun];
TH1 *hfit_R_rphi_central_neg[nrun], *hfit_R_zr_central_neg[nrun], *hfit_P_central_neg[nrun], *hfit_Z_central_neg[nrun];
TH1 *hfit_R_rphi_west_neg[nrun], *hfit_R_zr_west_neg[nrun], *hfit_P_west_neg[nrun], *hfit_Z_west_neg[nrun];
TH1 *hfit_R_rphi_east_neg[nrun], *hfit_R_zr_east_neg[nrun], *hfit_P_east_neg[nrun], *hfit_Z_east_neg[nrun];
TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];
double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
double ymax_RP=-100, ymin_RP=100;
for (int i = 0; i < nrun; i++)
{

  TFile* infile_dr_central = new TFile(Form("central/Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"");
  infile_dr_central->cd();
  h_R_central_pos[i] = (TH1*) infile_dr_central->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h_R_central_neg[i] = (TH1*) infile_dr_central->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  TFile* infile_dphi_central = new TFile(Form("central/Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"");
  infile_dphi_central->cd();
  h_P_central_pos[i] = (TH1*) infile_dphi_central->Get(Form("hIntDistortionP_posz_%d",runs[i]));
  h_P_central_neg[i] = (TH1*) infile_dphi_central->Get(Form("hIntDistortionP_negz_%d",runs[i]));
 
  TFile* infile_dz_central = new TFile(Form("central/Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"");
  infile_dz_central->cd();
  h_Z_central_pos[i] = (TH1*) infile_dz_central->Get(Form("hIntDistortionZ_posz_%d",runs[i]));
  h_Z_central_neg[i] = (TH1*) infile_dz_central->Get(Form("hIntDistortionZ_negz_%d",runs[i]));

  TFile* infile_dr_west = new TFile(Form("west/Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"");
  infile_dr_west->cd();
  h_R_west_pos[i] = (TH1*) infile_dr_west->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h_R_west_neg[i] = (TH1*) infile_dr_west->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  TFile* infile_dphi_west = new TFile(Form("west/Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"");
  infile_dphi_west->cd();
  h_P_west_pos[i] = (TH1*) infile_dphi_west->Get(Form("hIntDistortionP_posz_%d",runs[i]));
  h_P_west_neg[i] = (TH1*) infile_dphi_west->Get(Form("hIntDistortionP_negz_%d",runs[i]));

  TFile* infile_dz_west = new TFile(Form("west/Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"");
  infile_dz_west->cd();
  h_Z_west_pos[i] = (TH1*) infile_dz_west->Get(Form("hIntDistortionZ_posz_%d",runs[i]));
  h_Z_west_neg[i] = (TH1*) infile_dz_west->Get(Form("hIntDistortionZ_negz_%d",runs[i]));

  TFile* infile_dr_east = new TFile(Form("east/Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"");
  infile_dr_east->cd();
  h_R_east_pos[i] = (TH1*) infile_dr_east->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h_R_east_neg[i] = (TH1*) infile_dr_east->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  TFile* infile_dphi_east = new TFile(Form("east/Rootfiles/hist_dphi_1Dz_R%d.root",(int)selectR),"");
  infile_dphi_east->cd();
  h_P_east_pos[i] = (TH1*) infile_dphi_east->Get(Form("hIntDistortionP_posz_%d",runs[i]));
  h_P_east_neg[i] = (TH1*) infile_dphi_east->Get(Form("hIntDistortionP_negz_%d",runs[i]));
  
  TFile* infile_dz_east = new TFile(Form("east/Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"");
  infile_dz_east->cd();
  h_Z_east_pos[i] = (TH1*) infile_dz_east->Get(Form("hIntDistortionZ_posz_%d",runs[i]));
  h_Z_east_neg[i] = (TH1*) infile_dz_east->Get(Form("hIntDistortionZ_negz_%d",runs[i]));

  TFile* infile_fit_dr_central = new TFile(Form("central/Rootfiles/hist_dr_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dr_central->cd();
  hfit_R_rphi_central_pos[i] = (TH1*) infile_fit_dr_central->Get(Form("h_dr_rphi_posz"));
  hfit_R_rphi_central_neg[i] = (TH1*) infile_fit_dr_central->Get(Form("h_dr_rphi_negz"));
  hfit_R_zr_central_pos[i] = (TH1*) infile_fit_dr_central->Get(Form("h_dr_zr_posz"));
  hfit_R_zr_central_neg[i] = (TH1*) infile_fit_dr_central->Get(Form("h_dr_zr_negz"));

  TFile* infile_fit_dphi_central = new TFile(Form("central/Rootfiles/hist_dphi_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dphi_central->cd();
  hfit_P_central_pos[i] = (TH1*) infile_fit_dphi_central->Get(Form("h_dphi_rphi_posz"));
  hfit_P_central_neg[i] = (TH1*) infile_fit_dphi_central->Get(Form("h_dphi_rphi_negz"));
 
  TFile* infile_fit_dz_central = new TFile(Form("central/Rootfiles/hist_dz_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dz_central->cd();
  hfit_Z_central_pos[i] = (TH1*) infile_fit_dz_central->Get(Form("h_dz_zr_posz"));
  hfit_Z_central_neg[i] = (TH1*) infile_fit_dz_central->Get(Form("h_dz_zr_negz"));

  TFile* infile_fit_dr_west = new TFile(Form("west/Rootfiles/hist_dr_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dr_west->cd();
  hfit_R_rphi_west_pos[i] = (TH1*) infile_fit_dr_west->Get(Form("h_dr_rphi_posz"));
  hfit_R_rphi_west_neg[i] = (TH1*) infile_fit_dr_west->Get(Form("h_dr_rphi_negz"));
  hfit_R_zr_west_pos[i] = (TH1*) infile_fit_dr_west->Get(Form("h_dr_zr_posz"));
  hfit_R_zr_west_neg[i] = (TH1*) infile_fit_dr_west->Get(Form("h_dr_zr_negz"));

  TFile* infile_fit_dphi_west = new TFile(Form("west/Rootfiles/hist_dphi_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dphi_west->cd();
  hfit_P_west_pos[i] = (TH1*) infile_fit_dphi_west->Get(Form("h_dphi_rphi_posz"));
  hfit_P_west_neg[i] = (TH1*) infile_fit_dphi_west->Get(Form("h_dphi_rphi_negz"));
 
  TFile* infile_fit_dz_west = new TFile(Form("west/Rootfiles/hist_dz_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dz_west->cd();
  hfit_Z_west_pos[i] = (TH1*) infile_fit_dz_west->Get(Form("h_dz_zr_posz"));
  hfit_Z_west_neg[i] = (TH1*) infile_fit_dz_west->Get(Form("h_dz_zr_negz"));

  TFile* infile_fit_dr_east = new TFile(Form("east/Rootfiles/hist_dr_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dr_east->cd();
  hfit_R_rphi_east_pos[i] = (TH1*) infile_fit_dr_east->Get(Form("h_dr_rphi_posz"));
  hfit_R_rphi_east_neg[i] = (TH1*) infile_fit_dr_east->Get(Form("h_dr_rphi_negz"));
  hfit_R_zr_east_pos[i] = (TH1*) infile_fit_dr_east->Get(Form("h_dr_zr_posz"));
  hfit_R_zr_east_neg[i] = (TH1*) infile_fit_dr_east->Get(Form("h_dr_zr_negz"));

  TFile* infile_fit_dphi_east = new TFile(Form("east/Rootfiles/hist_dphi_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dphi_east->cd();
  hfit_P_east_pos[i] = (TH1*) infile_fit_dphi_east->Get(Form("h_dphi_rphi_posz"));
  hfit_P_east_neg[i] = (TH1*) infile_fit_dphi_east->Get(Form("h_dphi_rphi_negz"));
 
  TFile* infile_fit_dz_east = new TFile(Form("east/Rootfiles/hist_dz_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dz_east->cd();
  hfit_Z_east_pos[i] = (TH1*) infile_fit_dz_east->Get(Form("h_dz_zr_posz"));
  hfit_Z_east_neg[i] = (TH1*) infile_fit_dz_east->Get(Form("h_dz_zr_negz"));

  h_R_central_pos[i]->SetLineColor(kRed); h_R_central_pos[i]->SetLineWidth(1); h_R_central_pos[i]->SetFillColor(0); h_R_central_pos[i]->SetMarkerColor(kRed);
  h_P_central_pos[i]->SetLineColor(kRed); h_P_central_pos[i]->SetLineWidth(1); h_P_central_pos[i]->SetFillColor(0); h_P_central_pos[i]->SetMarkerColor(kRed);
  h_Z_central_pos[i]->SetLineColor(kRed); h_Z_central_pos[i]->SetLineWidth(1); h_Z_central_pos[i]->SetFillColor(0); h_Z_central_pos[i]->SetMarkerColor(kRed);
  h_R_central_neg[i]->SetLineColor(kRed); h_R_central_neg[i]->SetLineWidth(1); h_R_central_neg[i]->SetFillColor(0); h_R_central_neg[i]->SetMarkerColor(kRed);
  h_P_central_neg[i]->SetLineColor(kRed); h_P_central_neg[i]->SetLineWidth(1); h_P_central_neg[i]->SetFillColor(0); h_P_central_neg[i]->SetMarkerColor(kRed);
  h_Z_central_neg[i]->SetLineColor(kRed); h_Z_central_neg[i]->SetLineWidth(1); h_Z_central_neg[i]->SetFillColor(0); h_Z_central_neg[i]->SetMarkerColor(kRed);
  hfit_R_rphi_central_pos[i]->SetLineColor(kOrange); hfit_R_rphi_central_pos[i]->SetLineWidth(1); hfit_R_rphi_central_pos[i]->SetFillColor(0); hfit_R_rphi_central_pos[i]->SetMarkerColor(kOrange);
  hfit_R_zr_central_pos[i]->SetLineColor(kMagenta); hfit_R_zr_central_pos[i]->SetLineWidth(1); hfit_R_zr_central_pos[i]->SetFillColor(0); hfit_R_zr_central_pos[i]->SetMarkerColor(kMagenta);
  hfit_P_central_pos[i]->SetLineColor(kOrange); hfit_P_central_pos[i]->SetLineWidth(1); hfit_P_central_pos[i]->SetFillColor(0); hfit_P_central_pos[i]->SetMarkerColor(kOrange);
  hfit_Z_central_pos[i]->SetLineColor(kMagenta); hfit_Z_central_pos[i]->SetLineWidth(1); hfit_Z_central_pos[i]->SetFillColor(0); hfit_Z_central_pos[i]->SetMarkerColor(kMagenta);
  hfit_R_rphi_central_neg[i]->SetLineColor(kOrange); hfit_R_rphi_central_neg[i]->SetLineWidth(1); hfit_R_rphi_central_neg[i]->SetFillColor(0); hfit_R_rphi_central_neg[i]->SetMarkerColor(kOrange);
  hfit_R_zr_central_neg[i]->SetLineColor(kMagenta); hfit_R_zr_central_neg[i]->SetLineWidth(1); hfit_R_zr_central_neg[i]->SetFillColor(0); hfit_R_zr_central_neg[i]->SetMarkerColor(kMagenta);
  hfit_P_central_neg[i]->SetLineColor(kOrange); hfit_P_central_neg[i]->SetLineWidth(1); hfit_P_central_neg[i]->SetFillColor(0); hfit_P_central_neg[i]->SetMarkerColor(kOrange);
  hfit_Z_central_neg[i]->SetLineColor(kMagenta); hfit_Z_central_neg[i]->SetLineWidth(1); hfit_Z_central_neg[i]->SetFillColor(0); hfit_Z_central_neg[i]->SetMarkerColor(kMagenta);

  h_R_west_pos[i]->SetLineColor(kGreen); h_R_west_pos[i]->SetLineWidth(1); h_R_west_pos[i]->SetFillColor(0); h_R_west_pos[i]->SetMarkerColor(kGreen);
  h_P_west_pos[i]->SetLineColor(kGreen); h_P_west_pos[i]->SetLineWidth(1); h_P_west_pos[i]->SetFillColor(0); h_P_west_pos[i]->SetMarkerColor(kGreen);
  h_Z_west_pos[i]->SetLineColor(kGreen); h_Z_west_pos[i]->SetLineWidth(1); h_Z_west_pos[i]->SetFillColor(0); h_Z_west_pos[i]->SetMarkerColor(kGreen);
  h_R_west_neg[i]->SetLineColor(kGreen); h_R_west_neg[i]->SetLineWidth(1); h_R_west_neg[i]->SetFillColor(0); h_R_west_neg[i]->SetMarkerColor(kGreen);
  h_P_west_neg[i]->SetLineColor(kGreen); h_P_west_neg[i]->SetLineWidth(1); h_P_west_neg[i]->SetFillColor(0); h_P_west_neg[i]->SetMarkerColor(kGreen);
  h_Z_west_neg[i]->SetLineColor(kGreen); h_Z_west_neg[i]->SetLineWidth(1); h_Z_west_neg[i]->SetFillColor(0); h_Z_west_neg[i]->SetMarkerColor(kGreen);
  hfit_R_rphi_west_pos[i]->SetLineColor(kGreen+4); hfit_R_rphi_west_pos[i]->SetLineWidth(1); hfit_R_rphi_west_pos[i]->SetFillColor(0); hfit_R_rphi_west_pos[i]->SetMarkerColor(kGreen+4);
  hfit_R_zr_west_pos[i]->SetLineColor(kYellow); hfit_R_zr_west_pos[i]->SetLineWidth(1); hfit_R_zr_west_pos[i]->SetFillColor(0); hfit_R_zr_west_pos[i]->SetMarkerColor(kYellow);
  hfit_P_west_pos[i]->SetLineColor(kGreen+4); hfit_P_west_pos[i]->SetLineWidth(1); hfit_P_west_pos[i]->SetFillColor(0); hfit_P_west_pos[i]->SetMarkerColor(kGreen+4);
  hfit_Z_west_pos[i]->SetLineColor(kYellow); hfit_Z_west_pos[i]->SetLineWidth(1); hfit_Z_west_pos[i]->SetFillColor(0); hfit_Z_west_pos[i]->SetMarkerColor(kYellow);
  hfit_R_rphi_west_neg[i]->SetLineColor(kGreen+4); hfit_R_rphi_west_neg[i]->SetLineWidth(1); hfit_R_rphi_west_neg[i]->SetFillColor(0); hfit_R_rphi_west_neg[i]->SetMarkerColor(kGreen+4);
  hfit_R_zr_west_neg[i]->SetLineColor(kYellow); hfit_R_zr_west_neg[i]->SetLineWidth(1); hfit_R_zr_west_neg[i]->SetFillColor(0); hfit_R_zr_west_neg[i]->SetMarkerColor(kYellow);
  hfit_P_west_neg[i]->SetLineColor(kGreen+4); hfit_P_west_neg[i]->SetLineWidth(1); hfit_P_west_neg[i]->SetFillColor(0); hfit_P_west_neg[i]->SetMarkerColor(kGreen+4);
  hfit_Z_west_neg[i]->SetLineColor(kYellow); hfit_Z_west_neg[i]->SetLineWidth(1); hfit_Z_west_neg[i]->SetFillColor(0); hfit_Z_west_neg[i]->SetMarkerColor(kYellow);

  h_R_east_pos[i]->SetLineColor(kBlue); h_R_east_pos[i]->SetLineWidth(1); h_R_east_pos[i]->SetFillColor(0); h_R_east_pos[i]->SetMarkerColor(kBlue);
  h_P_east_pos[i]->SetLineColor(kBlue); h_P_east_pos[i]->SetLineWidth(1); h_P_east_pos[i]->SetFillColor(0); h_P_east_pos[i]->SetMarkerColor(kBlue);
  h_Z_east_pos[i]->SetLineColor(kBlue); h_Z_east_pos[i]->SetLineWidth(1); h_Z_east_pos[i]->SetFillColor(0); h_Z_east_pos[i]->SetMarkerColor(kBlue);
  h_R_east_neg[i]->SetLineColor(kBlue); h_R_east_neg[i]->SetLineWidth(1); h_R_east_neg[i]->SetFillColor(0); h_R_east_neg[i]->SetMarkerColor(kBlue);
  h_P_east_neg[i]->SetLineColor(kBlue); h_P_east_neg[i]->SetLineWidth(1); h_P_east_neg[i]->SetFillColor(0); h_P_east_neg[i]->SetMarkerColor(kBlue);
  h_Z_east_neg[i]->SetLineColor(kBlue); h_Z_east_neg[i]->SetLineWidth(1); h_Z_east_neg[i]->SetFillColor(0); h_Z_east_neg[i]->SetMarkerColor(kBlue);
  hfit_R_rphi_east_pos[i]->SetLineColor(kViolet+2); hfit_R_rphi_east_pos[i]->SetLineWidth(1); hfit_R_rphi_east_pos[i]->SetFillColor(0); hfit_R_rphi_east_pos[i]->SetMarkerColor(kViolet+2);
  hfit_R_zr_east_pos[i]->SetLineColor(kCyan); hfit_R_zr_east_pos[i]->SetLineWidth(1); hfit_R_zr_east_pos[i]->SetFillColor(0); hfit_R_zr_east_pos[i]->SetMarkerColor(kCyan);
  hfit_P_east_pos[i]->SetLineColor(kViolet+2); hfit_P_east_pos[i]->SetLineWidth(1); hfit_P_east_pos[i]->SetFillColor(0); hfit_P_east_pos[i]->SetMarkerColor(kViolet+2);
  hfit_Z_east_pos[i]->SetLineColor(kCyan); hfit_Z_east_pos[i]->SetLineWidth(1); hfit_Z_east_pos[i]->SetFillColor(0); hfit_Z_east_pos[i]->SetMarkerColor(kCyan);
  hfit_R_rphi_east_neg[i]->SetLineColor(kViolet+2); hfit_R_rphi_east_neg[i]->SetLineWidth(1); hfit_R_rphi_east_neg[i]->SetFillColor(0); hfit_R_rphi_east_neg[i]->SetMarkerColor(kViolet+2);
  hfit_R_zr_east_neg[i]->SetLineColor(kCyan); hfit_R_zr_east_neg[i]->SetLineWidth(1); hfit_R_zr_east_neg[i]->SetFillColor(0); hfit_R_zr_east_neg[i]->SetMarkerColor(kCyan);
  hfit_P_east_neg[i]->SetLineColor(kViolet+2); hfit_P_east_neg[i]->SetLineWidth(1); hfit_P_east_neg[i]->SetFillColor(0); hfit_P_east_neg[i]->SetMarkerColor(kViolet+2);
  hfit_Z_east_neg[i]->SetLineColor(kCyan); hfit_Z_east_neg[i]->SetLineWidth(1); hfit_Z_east_neg[i]->SetFillColor(0); hfit_Z_east_neg[i]->SetMarkerColor(kCyan);

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_central_pos[i],h_R_west_pos[i],h_R_east_pos[i],h_R_central_neg[i],h_R_west_neg[i],h_R_east_neg[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_central_pos[i],h_P_west_pos[i],h_P_east_pos[i],h_P_central_neg[i],h_P_west_neg[i],h_P_east_neg[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_central_pos[i],h_Z_west_pos[i],h_Z_east_pos[i],h_Z_central_neg[i],h_Z_west_neg[i],h_Z_east_neg[i]});

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

TCanvas* can = new TCanvas("can","",6400,2400);
can->Divide(4,2);
can->cd(1);
for (int i=0; i<nrun; i++)
{
  h_P_central_pos[i]->Draw("hist,e,same");
  h_P_west_pos[i]->Draw("hist,e,same");
  h_P_east_pos[i]->Draw("hist,e,same");
  hfit_P_central_pos[i]->Draw("hist,e,same");
  hfit_P_west_pos[i]->Draw("hist,e,same");
  hfit_P_east_pos[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_nco_i_P[i]->Draw();
  l_nco_o_P[i]->Draw();
  l_nci_i_P[i]->Draw();
  l_nci_o_P[i]->Draw();
}
can->cd(2);
for (int i=0; i<nrun; i++)
{
  h_R_central_pos[i]->Draw("hist,e,same");
  h_R_west_pos[i]->Draw("hist,e,same");
  h_R_east_pos[i]->Draw("hist,e,same");
  hfit_R_rphi_central_pos[i]->Draw("hist,e,same");
  hfit_R_rphi_west_pos[i]->Draw("hist,e,same");
  hfit_R_rphi_east_pos[i]->Draw("hist,e,same");
  hfit_R_zr_central_pos[i]->Draw("hist,e,same");
  hfit_R_zr_west_pos[i]->Draw("hist,e,same");
  hfit_R_zr_east_pos[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_nco_i_R[i]->Draw();
  l_nco_o_R[i]->Draw();
  l_nci_i_R[i]->Draw();
  l_nci_o_R[i]->Draw();
}
can->cd(3);
for (int i=0; i<nrun; i++)
{
  h_Z_central_pos[i]->Draw("hist,e,same");
  h_Z_west_pos[i]->Draw("hist,e,same");
  h_Z_east_pos[i]->Draw("hist,e,same");
  hfit_Z_central_pos[i]->Draw("hist,e,same");
  hfit_Z_west_pos[i]->Draw("hist,e,same");
  hfit_Z_east_pos[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_nco_i_Z[i]->Draw();
  l_nco_o_Z[i]->Draw();
  l_nci_i_Z[i]->Draw();
  l_nci_o_Z[i]->Draw();
}
can->cd(4);
TLegend *legend = new TLegend(0.0, 0.0, 0.9, 0.9);
for (int i=0; i<nrun; i++)
{
  legend->SetHeader(Form("Run %d, R = %d cm",runs[i],(int)selectR));
  legend->AddEntry(h_R_central_pos[i], Form("Central tiles, matrix"), "l");
  legend->AddEntry(h_R_west_pos[i], Form("West tiles, matrix"), "l");
  legend->AddEntry(h_R_east_pos[i], Form("East tiles, matrix"), "l");
  legend->AddEntry(hfit_R_rphi_central_pos[i], Form("Central tiles, fit R-#phi"), "l");
  legend->AddEntry(hfit_R_rphi_west_pos[i], Form("West tiles, fit R-#phi"), "l");
  legend->AddEntry(hfit_R_rphi_east_pos[i], Form("East tiles, fit R-#phi"), "l");
  legend->AddEntry(hfit_R_zr_central_pos[i], Form("Central tiles, fit Z-R"), "l");
  legend->AddEntry(hfit_R_zr_west_pos[i], Form("West tiles, fit Z-R"), "l");
  legend->AddEntry(hfit_R_zr_east_pos[i], Form("East tiles, fit Z-R"), "l");
}
legend->SetTextSize(0.1);
legend->Draw();
can->cd(5);
for (int i=0; i<nrun; i++) 
{
  h_P_central_neg[i]->Draw("hist,e,same");
  h_P_west_neg[i]->Draw("hist,e,same");
  h_P_east_neg[i]->Draw("hist,e,same");
  hfit_P_central_neg[i]->Draw("hist,e,same");
  hfit_P_west_neg[i]->Draw("hist,e,same");
  hfit_P_east_neg[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_sco_i_P[i]->Draw();
  l_sco_o_P[i]->Draw();
  l_sci_i_P[i]->Draw();
  l_sci_o_P[i]->Draw();
}
can->cd(6);
for (int i=0; i<nrun; i++)
{
  h_R_central_neg[i]->Draw("hist,e,same");
  h_R_west_neg[i]->Draw("hist,e,same");
  h_R_east_neg[i]->Draw("hist,e,same");
  hfit_R_rphi_central_neg[i]->Draw("hist,e,same");
  hfit_R_rphi_west_neg[i]->Draw("hist,e,same");
  hfit_R_rphi_east_neg[i]->Draw("hist,e,same");
  hfit_R_zr_central_neg[i]->Draw("hist,e,same");
  hfit_R_zr_west_neg[i]->Draw("hist,e,same");
  hfit_R_zr_east_neg[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_sco_i_R[i]->Draw();
  l_sco_o_R[i]->Draw();
  l_sci_i_R[i]->Draw();
  l_sci_o_R[i]->Draw();
}
can->cd(7);
for (int i=0; i<nrun; i++)
{
  h_Z_central_neg[i]->Draw("hist,e,same");
  h_Z_west_neg[i]->Draw("hist,e,same");
  h_Z_east_neg[i]->Draw("hist,e,same");
  hfit_Z_central_neg[i]->Draw("hist,e,same");
  hfit_Z_west_neg[i]->Draw("hist,e,same");
  hfit_Z_east_neg[i]->Draw("hist,e,same");
}
for (int i=0; i<nrun; i++)
{
  l_sco_i_Z[i]->Draw();
  l_sco_o_Z[i]->Draw();
  l_sci_i_Z[i]->Draw();
  l_sci_o_Z[i]->Draw();
}
gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d.pdf",(int)selectR));

delete can;
}

}
