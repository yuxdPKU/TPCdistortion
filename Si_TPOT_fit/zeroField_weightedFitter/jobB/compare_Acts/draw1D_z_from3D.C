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

int run = 52077;

TLine *l_nco_o_R, *l_nco_i_R, *l_nci_o_R, *l_nci_i_R;
TLine *l_nco_o_P, *l_nco_i_P, *l_nci_o_P, *l_nci_i_P;
TLine *l_nco_o_Z, *l_nco_i_Z, *l_nci_o_Z, *l_nci_i_Z;
TLine *l_sci_i_R, *l_sci_o_R, *l_sco_i_R, *l_sco_o_R;
TLine *l_sci_i_P, *l_sci_o_P, *l_sco_i_P, *l_sco_o_P;
TLine *l_sci_i_Z, *l_sci_o_Z, *l_sco_i_Z, *l_sco_o_Z;

TH1 *h_R_pos_Acts, *h_P_pos_Acts, *h_Z_pos_Acts;
TH1 *h_R_neg_Acts, *h_P_neg_Acts, *h_Z_neg_Acts;
TH1 *h_R_pos_WF, *h_P_pos_WF, *h_Z_pos_WF;
TH1 *h_R_neg_WF, *h_P_neg_WF, *h_Z_neg_WF;

TFile* ifile_WF_dr = new TFile(Form("../Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"");
ifile_WF_dr->cd();
h_R_pos_WF = (TH1*) ifile_WF_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_WF = (TH1*) ifile_WF_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_WF->SetName(Form("hIntDistortionR_posz_%d_WF",run));
h_R_neg_WF->SetName(Form("hIntDistortionR_negz_%d_WF",run));

TFile* ifile_WF_rdphi = new TFile(Form("../Rootfiles/hist_rdphi_1Dz_R%d.root",(int)selectR),"");
ifile_WF_rdphi->cd();
h_P_pos_WF = (TH1*) ifile_WF_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_WF = (TH1*) ifile_WF_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_WF->SetName(Form("hIntDistortionP_posz_%d_WF",run));
h_P_neg_WF->SetName(Form("hIntDistortionP_negz_%d_WF",run));

TFile* ifile_WF_dz = new TFile(Form("../Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"");
ifile_WF_dz->cd();
h_Z_pos_WF = (TH1*) ifile_WF_dz->Get(Form("hIntDistortionZ_posz_%d",run));
h_Z_neg_WF = (TH1*) ifile_WF_dz->Get(Form("hIntDistortionZ_negz_%d",run));
h_Z_pos_WF->SetName(Form("hIntDistortionZ_posz_%d_WF",run));
h_Z_neg_WF->SetName(Form("hIntDistortionZ_negz_%d_WF",run));

TFile* ifile_Acts_dr = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_dr_1Dz_R%d.root",(int)selectR),"");
ifile_Acts_dr->cd();
h_R_pos_Acts = (TH1*) ifile_Acts_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_Acts = (TH1*) ifile_Acts_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_Acts->SetName(Form("hIntDistortionR_posz_%d_Acts",run));
h_R_neg_Acts->SetName(Form("hIntDistortionR_negz_%d_Acts",run));

TFile* ifile_Acts_rdphi = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_rdphi_1Dz_R%d.root",(int)selectR),"");
ifile_Acts_rdphi->cd();
h_P_pos_Acts = (TH1*) ifile_Acts_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_Acts = (TH1*) ifile_Acts_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_Acts->SetName(Form("hIntDistortionP_posz_%d_Acts",run));
h_P_neg_Acts->SetName(Form("hIntDistortionP_negz_%d_Acts",run));

TFile* ifile_Acts_dz = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_dz_1Dz_R%d.root",(int)selectR),"");
ifile_Acts_dz->cd();
h_Z_pos_Acts = (TH1*) ifile_Acts_dz->Get(Form("hIntDistortionZ_posz_%d",run));
h_Z_neg_Acts = (TH1*) ifile_Acts_dz->Get(Form("hIntDistortionZ_negz_%d",run));
h_Z_pos_Acts->SetName(Form("hIntDistortionZ_posz_%d_Acts",run));
h_Z_neg_Acts->SetName(Form("hIntDistortionZ_negz_%d_Acts",run));

h_R_pos_WF->SetLineColor(2); h_R_pos_WF->SetLineWidth(1); h_R_pos_WF->SetFillColor(0); h_R_pos_WF->SetMarkerColor(2);
h_P_pos_WF->SetLineColor(2); h_P_pos_WF->SetLineWidth(1); h_P_pos_WF->SetFillColor(0); h_P_pos_WF->SetMarkerColor(2);
h_Z_pos_WF->SetLineColor(2); h_Z_pos_WF->SetLineWidth(1); h_Z_pos_WF->SetFillColor(0); h_Z_pos_WF->SetMarkerColor(2);
h_R_neg_WF->SetLineColor(2); h_R_neg_WF->SetLineWidth(1); h_R_neg_WF->SetFillColor(0); h_R_neg_WF->SetMarkerColor(2);
h_P_neg_WF->SetLineColor(2); h_P_neg_WF->SetLineWidth(1); h_P_neg_WF->SetFillColor(0); h_P_neg_WF->SetMarkerColor(2);
h_Z_neg_WF->SetLineColor(2); h_Z_neg_WF->SetLineWidth(1); h_Z_neg_WF->SetFillColor(0); h_Z_neg_WF->SetMarkerColor(2);

h_R_pos_Acts->SetLineColor(4); h_R_pos_Acts->SetLineWidth(1); h_R_pos_Acts->SetFillColor(0); h_R_pos_Acts->SetMarkerColor(4);
h_P_pos_Acts->SetLineColor(4); h_P_pos_Acts->SetLineWidth(1); h_P_pos_Acts->SetFillColor(0); h_P_pos_Acts->SetMarkerColor(4);
h_Z_pos_Acts->SetLineColor(4); h_Z_pos_Acts->SetLineWidth(1); h_Z_pos_Acts->SetFillColor(0); h_Z_pos_Acts->SetMarkerColor(4);
h_R_neg_Acts->SetLineColor(4); h_R_neg_Acts->SetLineWidth(1); h_R_neg_Acts->SetFillColor(0); h_R_neg_Acts->SetMarkerColor(4);
h_P_neg_Acts->SetLineColor(4); h_P_neg_Acts->SetLineWidth(1); h_P_neg_Acts->SetFillColor(0); h_P_neg_Acts->SetMarkerColor(4);
h_Z_neg_Acts->SetLineColor(4); h_Z_neg_Acts->SetLineWidth(1); h_Z_neg_Acts->SetFillColor(0); h_Z_neg_Acts->SetMarkerColor(4);


std::pair<double,double> yrange_R_neg = SetCommonYRange({h_R_neg_WF, h_R_neg_Acts});
std::pair<double,double> yrange_P_neg = SetCommonYRange({h_P_neg_WF, h_P_neg_Acts});
std::pair<double,double> yrange_Z_neg = SetCommonYRange({h_Z_neg_WF, h_Z_neg_Acts});
std::pair<double,double> yrange_R_pos = SetCommonYRange({h_R_pos_WF, h_R_pos_Acts});
std::pair<double,double> yrange_P_pos = SetCommonYRange({h_P_pos_WF, h_P_pos_Acts});
std::pair<double,double> yrange_Z_pos = SetCommonYRange({h_Z_pos_WF, h_Z_pos_Acts});


l_nco_o_R = new TLine(selectR*tan(phirange_NCO.second),yrange_R_pos.first,selectR*tan(phirange_NCO.second),yrange_R_pos.second);
l_nco_o_P = new TLine(selectR*tan(phirange_NCO.second),yrange_P_pos.first,selectR*tan(phirange_NCO.second),yrange_P_pos.second);
l_nco_o_Z = new TLine(selectR*tan(phirange_NCO.second),yrange_Z_pos.first,selectR*tan(phirange_NCO.second),yrange_Z_pos.second);

l_nco_i_R = new TLine(selectR*tan(phirange_NCO.first),yrange_R_pos.first,selectR*tan(phirange_NCO.first),yrange_R_pos.second);
l_nco_i_P = new TLine(selectR*tan(phirange_NCO.first),yrange_P_pos.first,selectR*tan(phirange_NCO.first),yrange_P_pos.second);
l_nco_i_Z = new TLine(selectR*tan(phirange_NCO.first),yrange_Z_pos.first,selectR*tan(phirange_NCO.first),yrange_Z_pos.second);

l_nci_o_R = new TLine(selectR*tan(phirange_NCI.second),yrange_R_pos.first,selectR*tan(phirange_NCI.second),yrange_R_pos.second);
l_nci_o_P = new TLine(selectR*tan(phirange_NCI.second),yrange_P_pos.first,selectR*tan(phirange_NCI.second),yrange_P_pos.second);
l_nci_o_Z = new TLine(selectR*tan(phirange_NCI.second),yrange_Z_pos.first,selectR*tan(phirange_NCI.second),yrange_Z_pos.second);

l_nci_i_R = new TLine(selectR*tan(phirange_NCI.first),yrange_R_pos.first,selectR*tan(phirange_NCI.first),yrange_R_pos.second);
l_nci_i_P = new TLine(selectR*tan(phirange_NCI.first),yrange_P_pos.first,selectR*tan(phirange_NCI.first),yrange_P_pos.second);
l_nci_i_Z = new TLine(selectR*tan(phirange_NCI.first),yrange_Z_pos.first,selectR*tan(phirange_NCI.first),yrange_Z_pos.second);

l_sci_i_R = new TLine(selectR*tan(phirange_SCI.second),yrange_R_neg.first,selectR*tan(phirange_SCI.second),yrange_R_neg.second);
l_sci_i_P = new TLine(selectR*tan(phirange_SCI.second),yrange_P_neg.first,selectR*tan(phirange_SCI.second),yrange_P_neg.second);
l_sci_i_Z = new TLine(selectR*tan(phirange_SCI.second),yrange_Z_neg.first,selectR*tan(phirange_SCI.second),yrange_Z_neg.second);

l_sci_o_R = new TLine(selectR*tan(phirange_SCI.first),yrange_R_neg.first,selectR*tan(phirange_SCI.first),yrange_R_neg.second);
l_sci_o_P = new TLine(selectR*tan(phirange_SCI.first),yrange_P_neg.first,selectR*tan(phirange_SCI.first),yrange_P_neg.second);
l_sci_o_Z = new TLine(selectR*tan(phirange_SCI.first),yrange_Z_neg.first,selectR*tan(phirange_SCI.first),yrange_Z_neg.second);

l_sco_i_R = new TLine(selectR*tan(phirange_SCO.second),yrange_R_neg.first,selectR*tan(phirange_SCO.second),yrange_R_neg.second);
l_sco_i_P = new TLine(selectR*tan(phirange_SCO.second),yrange_P_neg.first,selectR*tan(phirange_SCO.second),yrange_P_neg.second);
l_sco_i_Z = new TLine(selectR*tan(phirange_SCO.second),yrange_Z_neg.first,selectR*tan(phirange_SCO.second),yrange_Z_neg.second);

l_sco_o_R = new TLine(selectR*tan(phirange_SCO.first),yrange_R_neg.first,selectR*tan(phirange_SCO.first),yrange_R_neg.second);
l_sco_o_P = new TLine(selectR*tan(phirange_SCO.first),yrange_P_neg.first,selectR*tan(phirange_SCO.first),yrange_P_neg.second);
l_sco_o_Z = new TLine(selectR*tan(phirange_SCO.first),yrange_Z_neg.first,selectR*tan(phirange_SCO.first),yrange_Z_neg.second);

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos_Acts->Draw("hist,e,same");
h_P_pos_WF->Draw("hist,e,same");
l_nco_i_P->Draw();
l_nco_o_P->Draw();
l_nci_i_P->Draw();
l_nci_o_P->Draw();
can->cd(2);
gPad->SetLogy(0);
h_R_pos_Acts->Draw("hist,e,same");
h_R_pos_WF->Draw("hist,e,same");
l_nco_i_R->Draw();
l_nco_o_R->Draw();
l_nci_i_R->Draw();
l_nci_o_R->Draw();
can->cd(3);
gPad->SetLogy(0);
h_Z_pos_Acts->Draw("hist,e,same");
h_Z_pos_WF->Draw("hist,e,same");
l_nco_i_Z->Draw();
l_nco_o_Z->Draw();
l_nci_i_Z->Draw();
l_nci_o_Z->Draw();
can->cd(4);
TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d, R=%d cm",run,(int)selectR));
legend->AddEntry(h_P_pos_Acts, Form("ActsTrkFitter"), "l");
legend->AddEntry(h_P_pos_WF, Form("WeightedFitter"), "l");
legend->SetTextSize(0.1);
legend->Draw();
can->cd(5);
gPad->SetLogy(0);
h_P_neg_Acts->Draw("hist,e,same");
h_P_neg_WF->Draw("hist,e,same");
l_sco_i_P->Draw();
l_sco_o_P->Draw();
l_sci_i_P->Draw();
l_sci_o_P->Draw();
can->cd(6);
gPad->SetLogy(0);
h_R_neg_Acts->Draw("hist,e,same");
h_R_neg_WF->Draw("hist,e,same");
l_sco_i_R->Draw();
l_sco_o_R->Draw();
l_sci_i_R->Draw();
l_sci_o_R->Draw();
can->cd(7);
gPad->SetLogy(0);
h_Z_neg_Acts->Draw("hist,e,same");
h_Z_neg_WF->Draw("hist,e,same");
l_sco_i_Z->Draw();
l_sco_o_Z->Draw();
l_sci_i_Z->Draw();
l_sci_o_Z->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsZ_from3D_atR%d.pdf",(int)selectR));

delete can;
}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
//        yMin = TMath::Min(yMin, h->GetMinimum());
//        yMax = TMath::Max(yMax, h->GetMaximum());
        yMin = TMath::Min(yMin, h->GetBinContent(h->GetMinimumBin()));
        yMax = TMath::Max(yMax, h->GetBinContent(h->GetMaximumBin()));
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin > 0) yMin *= 0.9;
    if (yMin < 0) yMin *= 1.1;
    if (yMax > 0) yMax *= 1.1;
    if (yMax < 0) yMax *= 0.9;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
