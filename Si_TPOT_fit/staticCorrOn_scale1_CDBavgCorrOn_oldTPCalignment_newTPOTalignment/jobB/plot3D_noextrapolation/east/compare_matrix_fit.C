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

void compare_matrix_fit()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

std::pair<double, double> thetarange_NCO = {0.609928,0.917104};
std::pair<double, double> thetarange_NCI = {0.0276558,0.566484};
std::pair<double, double> thetarange_SCI = {-0.568656,-0.0325273};
//std::pair<double, double> thetarange_SCO = {-0.788212,-0.613404};
std::pair<double, double> thetarange_SCO = {-0.92,-0.613404};

std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
for (const auto& selectR : selectRs)
{

const int nrun = 1;
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

TH1 *h_R_rphi_pos[nrun], *h_P_rphi_pos[nrun];
TH1 *h_R_zr_pos[nrun], *h_Z_zr_pos[nrun];
TH1 *hfit_R_rphi_pos[nrun], *hfit_R_zr_pos[nrun], *hfit_P_pos[nrun], *hfit_Z_pos[nrun];

TH1 *h_R_rphi_neg[nrun], *h_P_rphi_neg[nrun];
TH1 *h_R_zr_neg[nrun], *h_Z_zr_neg[nrun];
TH1 *hfit_R_rphi_neg[nrun], *hfit_R_zr_neg[nrun], *hfit_P_neg[nrun], *hfit_Z_neg[nrun];

TLine *l_nco_o_R[nrun], *l_nco_i_R[nrun], *l_nci_o_R[nrun], *l_nci_i_R[nrun];
TLine *l_nco_o_P[nrun], *l_nco_i_P[nrun], *l_nci_o_P[nrun], *l_nci_i_P[nrun];
TLine *l_nco_o_Z[nrun], *l_nco_i_Z[nrun], *l_nci_o_Z[nrun], *l_nci_i_Z[nrun];

TLine *l_sci_i_R[nrun], *l_sci_o_R[nrun], *l_sco_i_R[nrun], *l_sco_o_R[nrun];
TLine *l_sci_i_P[nrun], *l_sci_o_P[nrun], *l_sco_i_P[nrun], *l_sco_o_P[nrun];
TLine *l_sci_i_Z[nrun], *l_sci_o_Z[nrun], *l_sco_i_Z[nrun], *l_sco_o_Z[nrun];

double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;

for (int i = 0; i < nrun; i++)
{

  TFile* infile_dr_rphi = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_phi.root",(int)selectR),"");
  infile_dr_rphi->cd();
  h_R_rphi_pos[i] = (TH1*) infile_dr_rphi->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h_R_rphi_neg[i] = (TH1*) infile_dr_rphi->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  TFile* infile_dphi_rphi = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d_phi.root",(int)selectR),"");
  infile_dphi_rphi->cd();
  h_P_rphi_pos[i] = (TH1*) infile_dphi_rphi->Get(Form("hIntDistortionP_posz_%d",runs[i]));
  h_P_rphi_neg[i] = (TH1*) infile_dphi_rphi->Get(Form("hIntDistortionP_negz_%d",runs[i]));
 
  TFile* infile_dr_zr = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_z.root",(int)selectR),"");
  infile_dr_zr->cd();
  h_R_zr_pos[i] = (TH1*) infile_dr_zr->Get(Form("hIntDistortionR_posz_%d",runs[i]));
  h_R_zr_neg[i] = (TH1*) infile_dr_zr->Get(Form("hIntDistortionR_negz_%d",runs[i]));

  TFile* infile_dz_zr = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d_z.root",(int)selectR),"");
  infile_dz_zr->cd();
  h_Z_zr_pos[i] = (TH1*) infile_dz_zr->Get(Form("hIntDistortionZ_posz_%d",runs[i]));
  h_Z_zr_neg[i] = (TH1*) infile_dz_zr->Get(Form("hIntDistortionZ_negz_%d",runs[i]));

  TFile* infile_fit_dr = new TFile(Form("Rootfiles/hist_dr_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dr->cd();
  hfit_R_rphi_pos[i] = (TH1*) infile_fit_dr->Get(Form("h_dr_rphi_posz"));
  hfit_R_rphi_neg[i] = (TH1*) infile_fit_dr->Get(Form("h_dr_rphi_negz"));
  hfit_R_zr_pos[i] = (TH1*) infile_fit_dr->Get(Form("h_dr_zr_posz"));
  hfit_R_zr_neg[i] = (TH1*) infile_fit_dr->Get(Form("h_dr_zr_negz"));

  TFile* infile_fit_dphi = new TFile(Form("Rootfiles/hist_dphi_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dphi->cd();
  hfit_P_pos[i] = (TH1*) infile_fit_dphi->Get(Form("h_dphi_rphi_posz"));
  hfit_P_neg[i] = (TH1*) infile_fit_dphi->Get(Form("h_dphi_rphi_negz"));
 
  TFile* infile_fit_dz = new TFile(Form("Rootfiles/hist_dz_1Dz_R%d_fit.root",(int)selectR),"");
  infile_fit_dz->cd();
  hfit_Z_pos[i] = (TH1*) infile_fit_dz->Get(Form("h_dz_zr_posz"));
  hfit_Z_neg[i] = (TH1*) infile_fit_dz->Get(Form("h_dz_zr_negz"));

  h_R_rphi_pos[i]->SetLineColor(kOrange+2); h_R_rphi_pos[i]->SetLineWidth(1); h_R_rphi_pos[i]->SetFillColor(0); h_R_rphi_pos[i]->SetMarkerColor(kOrange+2);
  h_P_rphi_pos[i]->SetLineColor(kOrange+2); h_P_rphi_pos[i]->SetLineWidth(1); h_P_rphi_pos[i]->SetFillColor(0); h_P_rphi_pos[i]->SetMarkerColor(kOrange+2);
  h_R_zr_pos[i]->SetLineColor(kMagenta+2); h_R_zr_pos[i]->SetLineWidth(1); h_R_zr_pos[i]->SetFillColor(0); h_R_zr_pos[i]->SetMarkerColor(kMagenta+2);
  h_Z_zr_pos[i]->SetLineColor(kMagenta+2); h_Z_zr_pos[i]->SetLineWidth(1); h_Z_zr_pos[i]->SetFillColor(0); h_Z_zr_pos[i]->SetMarkerColor(kMagenta+2);

  h_R_rphi_neg[i]->SetLineColor(kOrange+2); h_R_rphi_neg[i]->SetLineWidth(1); h_R_rphi_neg[i]->SetFillColor(0); h_R_rphi_neg[i]->SetMarkerColor(kOrange+2);
  h_P_rphi_neg[i]->SetLineColor(kOrange+2); h_P_rphi_neg[i]->SetLineWidth(1); h_P_rphi_neg[i]->SetFillColor(0); h_P_rphi_neg[i]->SetMarkerColor(kOrange+2);
  h_R_zr_neg[i]->SetLineColor(kMagenta+2); h_R_zr_neg[i]->SetLineWidth(1); h_R_zr_neg[i]->SetFillColor(0); h_R_zr_neg[i]->SetMarkerColor(kMagenta+2);
  h_Z_zr_neg[i]->SetLineColor(kMagenta+2); h_Z_zr_neg[i]->SetLineWidth(1); h_Z_zr_neg[i]->SetFillColor(0); h_Z_zr_neg[i]->SetMarkerColor(kMagenta+2);

  hfit_R_rphi_pos[i]->SetLineColor(kOrange); hfit_R_rphi_pos[i]->SetLineWidth(1); hfit_R_rphi_pos[i]->SetFillColor(0); hfit_R_rphi_pos[i]->SetMarkerColor(kOrange);
  hfit_R_zr_pos[i]->SetLineColor(kMagenta); hfit_R_zr_pos[i]->SetLineWidth(1); hfit_R_zr_pos[i]->SetFillColor(0); hfit_R_zr_pos[i]->SetMarkerColor(kMagenta);
  hfit_P_pos[i]->SetLineColor(kOrange); hfit_P_pos[i]->SetLineWidth(1); hfit_P_pos[i]->SetFillColor(0); hfit_P_pos[i]->SetMarkerColor(kOrange);
  hfit_Z_pos[i]->SetLineColor(kMagenta); hfit_Z_pos[i]->SetLineWidth(1); hfit_Z_pos[i]->SetFillColor(0); hfit_Z_pos[i]->SetMarkerColor(kMagenta);

  hfit_R_rphi_neg[i]->SetLineColor(kOrange); hfit_R_rphi_neg[i]->SetLineWidth(1); hfit_R_rphi_neg[i]->SetFillColor(0); hfit_R_rphi_neg[i]->SetMarkerColor(kOrange);
  hfit_R_zr_neg[i]->SetLineColor(kMagenta); hfit_R_zr_neg[i]->SetLineWidth(1); hfit_R_zr_neg[i]->SetFillColor(0); hfit_R_zr_neg[i]->SetMarkerColor(kMagenta);
  hfit_P_neg[i]->SetLineColor(kOrange); hfit_P_neg[i]->SetLineWidth(1); hfit_P_neg[i]->SetFillColor(0); hfit_P_neg[i]->SetMarkerColor(kOrange);
  hfit_Z_neg[i]->SetLineColor(kMagenta); hfit_Z_neg[i]->SetLineWidth(1); hfit_Z_neg[i]->SetFillColor(0); hfit_Z_neg[i]->SetMarkerColor(kMagenta);

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_rphi_pos[i],h_R_zr_pos[i],h_R_rphi_neg[i],h_R_zr_neg[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_rphi_pos[i],h_P_rphi_neg[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_zr_pos[i],h_Z_zr_neg[i]});

  l_nco_o_R[i] = new TLine(selectR*tan(thetarange_NCO.second),yrange_R.first,selectR*tan(thetarange_NCO.second),yrange_R.second);
  l_nco_o_P[i] = new TLine(selectR*tan(thetarange_NCO.second),yrange_P.first,selectR*tan(thetarange_NCO.second),yrange_P.second);
  l_nco_o_Z[i] = new TLine(selectR*tan(thetarange_NCO.second),yrange_Z.first,selectR*tan(thetarange_NCO.second),yrange_Z.second);

  l_nco_i_R[i] = new TLine(selectR*tan(thetarange_NCO.first),yrange_R.first,selectR*tan(thetarange_NCO.first),yrange_R.second);
  l_nco_i_P[i] = new TLine(selectR*tan(thetarange_NCO.first),yrange_P.first,selectR*tan(thetarange_NCO.first),yrange_P.second);
  l_nco_i_Z[i] = new TLine(selectR*tan(thetarange_NCO.first),yrange_Z.first,selectR*tan(thetarange_NCO.first),yrange_Z.second);

  l_nci_o_R[i] = new TLine(selectR*tan(thetarange_NCI.second),yrange_R.first,selectR*tan(thetarange_NCI.second),yrange_R.second);
  l_nci_o_P[i] = new TLine(selectR*tan(thetarange_NCI.second),yrange_P.first,selectR*tan(thetarange_NCI.second),yrange_P.second);
  l_nci_o_Z[i] = new TLine(selectR*tan(thetarange_NCI.second),yrange_Z.first,selectR*tan(thetarange_NCI.second),yrange_Z.second);

  l_nci_i_R[i] = new TLine(selectR*tan(thetarange_NCI.first),yrange_R.first,selectR*tan(thetarange_NCI.first),yrange_R.second);
  l_nci_i_P[i] = new TLine(selectR*tan(thetarange_NCI.first),yrange_P.first,selectR*tan(thetarange_NCI.first),yrange_P.second);
  l_nci_i_Z[i] = new TLine(selectR*tan(thetarange_NCI.first),yrange_Z.first,selectR*tan(thetarange_NCI.first),yrange_Z.second);

  l_sci_i_R[i] = new TLine(selectR*tan(thetarange_SCI.second),yrange_R.first,selectR*tan(thetarange_SCI.second),yrange_R.second);
  l_sci_i_P[i] = new TLine(selectR*tan(thetarange_SCI.second),yrange_P.first,selectR*tan(thetarange_SCI.second),yrange_P.second);
  l_sci_i_Z[i] = new TLine(selectR*tan(thetarange_SCI.second),yrange_Z.first,selectR*tan(thetarange_SCI.second),yrange_Z.second);

  l_sci_o_R[i] = new TLine(selectR*tan(thetarange_SCI.first),yrange_R.first,selectR*tan(thetarange_SCI.first),yrange_R.second);
  l_sci_o_P[i] = new TLine(selectR*tan(thetarange_SCI.first),yrange_P.first,selectR*tan(thetarange_SCI.first),yrange_P.second);
  l_sci_o_Z[i] = new TLine(selectR*tan(thetarange_SCI.first),yrange_Z.first,selectR*tan(thetarange_SCI.first),yrange_Z.second);

  l_sco_i_R[i] = new TLine(selectR*tan(thetarange_SCO.second),yrange_R.first,selectR*tan(thetarange_SCO.second),yrange_R.second);
  l_sco_i_P[i] = new TLine(selectR*tan(thetarange_SCO.second),yrange_P.first,selectR*tan(thetarange_SCO.second),yrange_P.second);
  l_sco_i_Z[i] = new TLine(selectR*tan(thetarange_SCO.second),yrange_Z.first,selectR*tan(thetarange_SCO.second),yrange_Z.second);

  l_sco_o_R[i] = new TLine(selectR*tan(thetarange_SCO.first),yrange_R.first,selectR*tan(thetarange_SCO.first),yrange_R.second);
  l_sco_o_P[i] = new TLine(selectR*tan(thetarange_SCO.first),yrange_P.first,selectR*tan(thetarange_SCO.first),yrange_P.second);
  l_sco_o_Z[i] = new TLine(selectR*tan(thetarange_SCO.first),yrange_Z.first,selectR*tan(thetarange_SCO.first),yrange_Z.second);
}

TCanvas* can = new TCanvas("can","",6400,2400);
can->Divide(4,2);
can->cd(1);
for (int i=0; i<nrun; i++)
{
  h_P_rphi_pos[i]->Draw("hist,e,same");
  hfit_P_pos[i]->Draw("hist,e,same");
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
  h_R_rphi_pos[i]->Draw("hist,e,same");
  h_R_zr_pos[i]->Draw("hist,e,same");
  hfit_R_rphi_pos[i]->Draw("hist,e,same");
  hfit_R_zr_pos[i]->Draw("hist,e,same");
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
  //h_Z_pos[i]->SetMinimum(-0.5);
  //h_Z_pos[i]->SetMaximum(2.5);
  h_Z_zr_pos[i]->Draw("hist,e,same");
  hfit_Z_pos[i]->Draw("hist,e,same");
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
  legend->AddEntry(h_R_rphi_pos[i], Form("Matrix Inversion (#phi)"), "l");
  legend->AddEntry(h_R_zr_pos[i], Form("Matrix Inversion (z)"), "l");
  legend->AddEntry(hfit_R_rphi_pos[i], Form("Fit R-#phi"), "l");
  legend->AddEntry(hfit_R_zr_pos[i], Form("Fit Z-R"), "l");
}
legend->SetTextSize(0.1);
legend->Draw();
can->cd(5);
for (int i=0; i<nrun; i++) 
{
  h_P_rphi_neg[i]->Draw("hist,e,same");
  hfit_P_neg[i]->Draw("hist,e,same");
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
  h_R_rphi_neg[i]->Draw("hist,e,same");
  h_R_zr_neg[i]->Draw("hist,e,same");
  hfit_R_rphi_neg[i]->Draw("hist,e,same");
  hfit_R_zr_neg[i]->Draw("hist,e,same");
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
  //h_Z_neg[i]->SetMinimum(-2.5);
  //h_Z_neg[i]->SetMaximum(0.5);
  h_Z_zr_neg[i]->Draw("hist,e,same");
  hfit_Z_neg[i]->Draw("hist,e,same");
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
can->SaveAs(Form("figure/compare_fit_vs_matrix_vsZ_from3D_atR%d.pdf",(int)selectR));

delete can;
}

}
