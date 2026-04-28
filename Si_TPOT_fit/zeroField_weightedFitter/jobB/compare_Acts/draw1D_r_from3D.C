#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, bool convert_P_2_RP)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      double value = h3->GetBinContent(xbin, i, zbin);
      double error = h3->GetBinError(xbin, i, zbin);
      if (convert_P_2_RP)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
	//cout<<"rcenter = "<<rcenter<<" , i = "<<i<<endl;
        value *= rcenter;
        error *= rcenter;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void plot1D(TH2* h2, TH1* h1)
{
  int xbin = h2->GetXaxis()->FindBin(4.55073);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      int xminbin = h2->GetXaxis()->FindBin(4.55073);
      int xmaxbin = h2->GetXaxis()->FindBin(4.84711);
      float bincontent=0, binerror=0;
      int nbin=0;
      for (int j = xminbin; j <= xmaxbin; j++)
      {
        bincontent+=h2->GetBinContent(j, i);
        binerror+=h2->GetBinError(j, i);
	nbin++;
      }
      bincontent/=nbin;
      binerror/=nbin;
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void draw1D_r_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

// runNo - CAD-MBD-NS CAD-ZDC-NS
// 53877 - 400khz 4758.536231884057
// 53876 - 430khz 5082.326086956524
// 53756 - 380khz 4471.421428571427
// 53744 - 300khz 3581.862318840581
// 53630 - 550khz 6849.317241379308
// 53534 - 250khz 3013.5338345864657
// 53285 - 70khz 787.2765957446811

int run = 52077;

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

TH1 *h_R_pos_Acts, *h_P_pos_Acts, *h_Z_pos_Acts;
TH1 *h_R_neg_Acts, *h_P_neg_Acts, *h_Z_neg_Acts;
TH1 *h_R_pos_WF, *h_P_pos_WF, *h_Z_pos_WF;
TH1 *h_R_neg_WF, *h_P_neg_WF, *h_Z_neg_WF;
TH1 *h_R_pos_Lamination, *h_P_pos_Lamination, *h_Z_pos_Lamination;
TH1 *h_R_neg_Lamination, *h_P_neg_Lamination, *h_Z_neg_Lamination;

TFile* ifile_WF_dr = new TFile(Form("../Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"");
ifile_WF_dr->cd();
h_R_pos_WF = (TH1*) ifile_WF_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_WF = (TH1*) ifile_WF_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_WF->SetName(Form("hIntDistortionR_posz_%d_WF",run));
h_R_neg_WF->SetName(Form("hIntDistortionR_negz_%d_WF",run));

TFile* ifile_WF_rdphi = new TFile(Form("../Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"");
ifile_WF_rdphi->cd();
h_P_pos_WF = (TH1*) ifile_WF_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_WF = (TH1*) ifile_WF_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_WF->SetName(Form("hIntDistortionP_posz_%d_WF",run));
h_P_neg_WF->SetName(Form("hIntDistortionP_negz_%d_WF",run));

TFile* ifile_WF_dz = new TFile(Form("../Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"");
ifile_WF_dz->cd();
h_Z_pos_WF = (TH1*) ifile_WF_dz->Get(Form("hIntDistortionZ_posz_%d",run));
h_Z_neg_WF = (TH1*) ifile_WF_dz->Get(Form("hIntDistortionZ_negz_%d",run));
h_Z_pos_WF->SetName(Form("hIntDistortionZ_posz_%d_WF",run));
h_Z_neg_WF->SetName(Form("hIntDistortionZ_negz_%d_WF",run));

TFile* ifile_Lamination_dr = new TFile(Form("../Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"");
ifile_Lamination_dr->cd();
h_R_pos_Lamination = (TH1*) ifile_Lamination_dr->Get(Form("hIntDistortionR_posz_%d_cdb",run));
h_R_neg_Lamination = (TH1*) ifile_Lamination_dr->Get(Form("hIntDistortionR_negz_%d_cdb",run));
h_R_pos_Lamination->SetName(Form("hIntDistortionR_posz_%d_Lamination",run));
h_R_neg_Lamination->SetName(Form("hIntDistortionR_negz_%d_Lamination",run));

TFile* ifile_Lamination_rdphi = new TFile(Form("../Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"");
ifile_Lamination_rdphi->cd();
h_P_pos_Lamination = (TH1*) ifile_Lamination_rdphi->Get(Form("hIntDistortionP_posz_%d_cdb",run));
h_P_neg_Lamination = (TH1*) ifile_Lamination_rdphi->Get(Form("hIntDistortionP_negz_%d_cdb",run));
h_P_pos_Lamination->SetName(Form("hIntDistortionP_posz_%d_Lamination",run));
h_P_neg_Lamination->SetName(Form("hIntDistortionP_negz_%d_Lamination",run));

TFile* ifile_Lamination_dz = new TFile(Form("../Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"");
ifile_Lamination_dz->cd();
h_Z_pos_Lamination = (TH1*) ifile_Lamination_dz->Get(Form("hIntDistortionZ_posz_%d_cdb",run));
h_Z_neg_Lamination = (TH1*) ifile_Lamination_dz->Get(Form("hIntDistortionZ_negz_%d_cdb",run));
h_Z_pos_Lamination->SetName(Form("hIntDistortionZ_posz_%d_Lamination",run));
h_Z_neg_Lamination->SetName(Form("hIntDistortionZ_negz_%d_Lamination",run));

TFile* ifile_Acts_dr = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"");
ifile_Acts_dr->cd();
h_R_pos_Acts = (TH1*) ifile_Acts_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_Acts = (TH1*) ifile_Acts_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_Acts->SetName(Form("hIntDistortionR_posz_%d_Acts",run));
h_R_neg_Acts->SetName(Form("hIntDistortionR_negz_%d_Acts",run));

TFile* ifile_Acts_rdphi = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"");
ifile_Acts_rdphi->cd();
h_P_pos_Acts = (TH1*) ifile_Acts_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_Acts = (TH1*) ifile_Acts_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_Acts->SetName(Form("hIntDistortionP_posz_%d_Acts",run));
h_P_neg_Acts->SetName(Form("hIntDistortionP_negz_%d_Acts",run));

TFile* ifile_Acts_dz = new TFile(Form("../../../zeroField/jobB/Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"");
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

h_R_pos_Lamination->SetLineColor(1); h_R_pos_Lamination->SetLineWidth(1); h_R_pos_Lamination->SetFillColor(0); h_R_pos_Lamination->SetMarkerColor(1);
h_P_pos_Lamination->SetLineColor(1); h_P_pos_Lamination->SetLineWidth(1); h_P_pos_Lamination->SetFillColor(0); h_P_pos_Lamination->SetMarkerColor(1);
h_Z_pos_Lamination->SetLineColor(1); h_Z_pos_Lamination->SetLineWidth(1); h_Z_pos_Lamination->SetFillColor(0); h_Z_pos_Lamination->SetMarkerColor(1);
h_R_neg_Lamination->SetLineColor(1); h_R_neg_Lamination->SetLineWidth(1); h_R_neg_Lamination->SetFillColor(0); h_R_neg_Lamination->SetMarkerColor(1);
h_P_neg_Lamination->SetLineColor(1); h_P_neg_Lamination->SetLineWidth(1); h_P_neg_Lamination->SetFillColor(0); h_P_neg_Lamination->SetMarkerColor(1);
h_Z_neg_Lamination->SetLineColor(1); h_Z_neg_Lamination->SetLineWidth(1); h_Z_neg_Lamination->SetFillColor(0); h_Z_neg_Lamination->SetMarkerColor(1);

h_R_pos_Acts->SetLineColor(4); h_R_pos_Acts->SetLineWidth(1); h_R_pos_Acts->SetFillColor(0); h_R_pos_Acts->SetMarkerColor(4);
h_P_pos_Acts->SetLineColor(4); h_P_pos_Acts->SetLineWidth(1); h_P_pos_Acts->SetFillColor(0); h_P_pos_Acts->SetMarkerColor(4);
h_Z_pos_Acts->SetLineColor(4); h_Z_pos_Acts->SetLineWidth(1); h_Z_pos_Acts->SetFillColor(0); h_Z_pos_Acts->SetMarkerColor(4);
h_R_neg_Acts->SetLineColor(4); h_R_neg_Acts->SetLineWidth(1); h_R_neg_Acts->SetFillColor(0); h_R_neg_Acts->SetMarkerColor(4);
h_P_neg_Acts->SetLineColor(4); h_P_neg_Acts->SetLineWidth(1); h_P_neg_Acts->SetFillColor(0); h_P_neg_Acts->SetMarkerColor(4);
h_Z_neg_Acts->SetLineColor(4); h_Z_neg_Acts->SetLineWidth(1); h_Z_neg_Acts->SetFillColor(0); h_Z_neg_Acts->SetMarkerColor(4);


std::pair<double,double> yrange_R_neg = SetCommonYRange({h_R_neg_WF, h_R_neg_Acts,h_R_neg_Lamination});
std::pair<double,double> yrange_P_neg = SetCommonYRange({h_P_neg_WF, h_P_neg_Acts, h_P_neg_Lamination});
std::pair<double,double> yrange_Z_neg = SetCommonYRange({h_Z_neg_WF, h_Z_neg_Acts, h_Z_neg_Lamination});
std::pair<double,double> yrange_R_pos = SetCommonYRange({h_R_pos_WF, h_R_pos_Acts, h_R_pos_Lamination});
std::pair<double,double> yrange_P_pos = SetCommonYRange({h_P_pos_WF, h_P_pos_Acts, h_P_pos_Lamination});
std::pair<double,double> yrange_Z_pos = SetCommonYRange({h_Z_pos_WF, h_Z_pos_Acts, h_Z_pos_Lamination});

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos_Acts->Draw("hist,e,same");
h_P_pos_WF->Draw("hist,e,same");
h_P_pos_Lamination->Draw("hist,same");
can->cd(2);
gPad->SetLogy(0);
h_R_pos_Acts->Draw("hist,e,same");
h_R_pos_WF->Draw("hist,e,same");
h_R_pos_Lamination->Draw("hist,same");
can->cd(3);
gPad->SetLogy(0);
h_Z_pos_Acts->Draw("hist,e,same");
h_Z_pos_WF->Draw("hist,e,same");
h_Z_pos_Lamination->Draw("hist,same");
can->cd(4);
TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d, |Z|=%d cm",run,(int)selectZ));
legend->AddEntry(h_P_pos_Acts, Form("ActsTrkFitter"), "l");
legend->AddEntry(h_P_pos_WF, Form("WeightedFitter"), "l");
legend->AddEntry(h_P_pos_Lamination, Form("Lamination"), "l");
legend->SetTextSize(0.1);
legend->Draw();
can->cd(5);
gPad->SetLogy(0);
h_P_neg_Acts->Draw("hist,e,same");
h_P_neg_WF->Draw("hist,e,same");
h_P_neg_Lamination->Draw("hist,same");
can->cd(6);
gPad->SetLogy(0);
h_R_neg_Acts->Draw("hist,e,same");
h_R_neg_WF->Draw("hist,e,same");
h_R_neg_Lamination->Draw("hist,same");
can->cd(7);
gPad->SetLogy(0);
h_Z_neg_Acts->Draw("hist,e,same");
h_Z_neg_WF->Draw("hist,e,same");
h_Z_neg_Lamination->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d.pdf",(int)selectZ));

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
