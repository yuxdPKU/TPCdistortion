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

TH1 *h_R_pos_noAvgCorr, *h_P_pos_noAvgCorr, *h_Z_pos_noAvgCorr;
TH1 *h_R_neg_noAvgCorr, *h_P_neg_noAvgCorr, *h_Z_neg_noAvgCorr;
TH1 *h_R_pos_withAvgCorr, *h_P_pos_withAvgCorr, *h_Z_pos_withAvgCorr;
TH1 *h_R_neg_withAvgCorr, *h_P_neg_withAvgCorr, *h_Z_neg_withAvgCorr;

TFile* ifile_withAvgCorr_dr = new TFile(Form("../Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"");
ifile_withAvgCorr_dr->cd();
h_R_pos_withAvgCorr = (TH1*) ifile_withAvgCorr_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_withAvgCorr = (TH1*) ifile_withAvgCorr_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_withAvgCorr->SetName(Form("hIntDistortionR_posz_%d_withAvgCorr",run));
h_R_neg_withAvgCorr->SetName(Form("hIntDistortionR_negz_%d_withAvgCorr",run));

TFile* ifile_withAvgCorr_rdphi = new TFile(Form("../Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"");
ifile_withAvgCorr_rdphi->cd();
h_P_pos_withAvgCorr = (TH1*) ifile_withAvgCorr_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_withAvgCorr = (TH1*) ifile_withAvgCorr_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_withAvgCorr->SetName(Form("hIntDistortionP_posz_%d_withAvgCorr",run));
h_P_neg_withAvgCorr->SetName(Form("hIntDistortionP_negz_%d_withAvgCorr",run));

TFile* ifile_withAvgCorr_dz = new TFile(Form("../Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"");
ifile_withAvgCorr_dz->cd();
h_Z_pos_withAvgCorr = (TH1*) ifile_withAvgCorr_dz->Get(Form("hIntDistortionZ_posz_%d",run));
h_Z_neg_withAvgCorr = (TH1*) ifile_withAvgCorr_dz->Get(Form("hIntDistortionZ_negz_%d",run));
h_Z_pos_withAvgCorr->SetName(Form("hIntDistortionZ_posz_%d_withAvgCorr",run));
h_Z_neg_withAvgCorr->SetName(Form("hIntDistortionZ_negz_%d_withAvgCorr",run));

TFile* ifile_noAvgCorr_dr = new TFile(Form("../../../zeroField_weightedFitter/jobB/Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"");
ifile_noAvgCorr_dr->cd();
h_R_pos_noAvgCorr = (TH1*) ifile_noAvgCorr_dr->Get(Form("hIntDistortionR_posz_%d",run));
h_R_neg_noAvgCorr = (TH1*) ifile_noAvgCorr_dr->Get(Form("hIntDistortionR_negz_%d",run));
h_R_pos_noAvgCorr->SetName(Form("hIntDistortionR_posz_%d_noAvgCorr",run));
h_R_neg_noAvgCorr->SetName(Form("hIntDistortionR_negz_%d_noAvgCorr",run));

TFile* ifile_noAvgCorr_rdphi = new TFile(Form("../../../zeroField_weightedFitter/jobB/Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"");
ifile_noAvgCorr_rdphi->cd();
h_P_pos_noAvgCorr = (TH1*) ifile_noAvgCorr_rdphi->Get(Form("hIntDistortionP_posz_%d",run));
h_P_neg_noAvgCorr = (TH1*) ifile_noAvgCorr_rdphi->Get(Form("hIntDistortionP_negz_%d",run));
h_P_pos_noAvgCorr->SetName(Form("hIntDistortionP_posz_%d_noAvgCorr",run));
h_P_neg_noAvgCorr->SetName(Form("hIntDistortionP_negz_%d_noAvgCorr",run));

TFile* ifile_noAvgCorr_dz = new TFile(Form("../../../zeroField_weightedFitter/jobB/Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"");
ifile_noAvgCorr_dz->cd();
h_Z_pos_noAvgCorr = (TH1*) ifile_noAvgCorr_dz->Get(Form("hIntDistortionZ_posz_%d",run));
h_Z_neg_noAvgCorr = (TH1*) ifile_noAvgCorr_dz->Get(Form("hIntDistortionZ_negz_%d",run));
h_Z_pos_noAvgCorr->SetName(Form("hIntDistortionZ_posz_%d_noAvgCorr",run));
h_Z_neg_noAvgCorr->SetName(Form("hIntDistortionZ_negz_%d_noAvgCorr",run));

h_R_pos_withAvgCorr->SetLineColor(2); h_R_pos_withAvgCorr->SetLineWidth(1); h_R_pos_withAvgCorr->SetFillColor(0); h_R_pos_withAvgCorr->SetMarkerColor(2);
h_P_pos_withAvgCorr->SetLineColor(2); h_P_pos_withAvgCorr->SetLineWidth(1); h_P_pos_withAvgCorr->SetFillColor(0); h_P_pos_withAvgCorr->SetMarkerColor(2);
h_Z_pos_withAvgCorr->SetLineColor(2); h_Z_pos_withAvgCorr->SetLineWidth(1); h_Z_pos_withAvgCorr->SetFillColor(0); h_Z_pos_withAvgCorr->SetMarkerColor(2);
h_R_neg_withAvgCorr->SetLineColor(2); h_R_neg_withAvgCorr->SetLineWidth(1); h_R_neg_withAvgCorr->SetFillColor(0); h_R_neg_withAvgCorr->SetMarkerColor(2);
h_P_neg_withAvgCorr->SetLineColor(2); h_P_neg_withAvgCorr->SetLineWidth(1); h_P_neg_withAvgCorr->SetFillColor(0); h_P_neg_withAvgCorr->SetMarkerColor(2);
h_Z_neg_withAvgCorr->SetLineColor(2); h_Z_neg_withAvgCorr->SetLineWidth(1); h_Z_neg_withAvgCorr->SetFillColor(0); h_Z_neg_withAvgCorr->SetMarkerColor(2);

h_R_pos_noAvgCorr->SetLineColor(4); h_R_pos_noAvgCorr->SetLineWidth(1); h_R_pos_noAvgCorr->SetFillColor(0); h_R_pos_noAvgCorr->SetMarkerColor(4);
h_P_pos_noAvgCorr->SetLineColor(4); h_P_pos_noAvgCorr->SetLineWidth(1); h_P_pos_noAvgCorr->SetFillColor(0); h_P_pos_noAvgCorr->SetMarkerColor(4);
h_Z_pos_noAvgCorr->SetLineColor(4); h_Z_pos_noAvgCorr->SetLineWidth(1); h_Z_pos_noAvgCorr->SetFillColor(0); h_Z_pos_noAvgCorr->SetMarkerColor(4);
h_R_neg_noAvgCorr->SetLineColor(4); h_R_neg_noAvgCorr->SetLineWidth(1); h_R_neg_noAvgCorr->SetFillColor(0); h_R_neg_noAvgCorr->SetMarkerColor(4);
h_P_neg_noAvgCorr->SetLineColor(4); h_P_neg_noAvgCorr->SetLineWidth(1); h_P_neg_noAvgCorr->SetFillColor(0); h_P_neg_noAvgCorr->SetMarkerColor(4);
h_Z_neg_noAvgCorr->SetLineColor(4); h_Z_neg_noAvgCorr->SetLineWidth(1); h_Z_neg_noAvgCorr->SetFillColor(0); h_Z_neg_noAvgCorr->SetMarkerColor(4);


std::pair<double,double> yrange_R_neg = SetCommonYRange({h_R_neg_withAvgCorr, h_R_neg_noAvgCorr});
std::pair<double,double> yrange_P_neg = SetCommonYRange({h_P_neg_withAvgCorr, h_P_neg_noAvgCorr});
std::pair<double,double> yrange_Z_neg = SetCommonYRange({h_Z_neg_withAvgCorr, h_Z_neg_noAvgCorr});
std::pair<double,double> yrange_R_pos = SetCommonYRange({h_R_pos_withAvgCorr, h_R_pos_noAvgCorr});
std::pair<double,double> yrange_P_pos = SetCommonYRange({h_P_pos_withAvgCorr, h_P_pos_noAvgCorr});
std::pair<double,double> yrange_Z_pos = SetCommonYRange({h_Z_pos_withAvgCorr, h_Z_pos_noAvgCorr});

TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos_noAvgCorr->Draw("hist,e,same");
h_P_pos_withAvgCorr->Draw("hist,e,same");
can->cd(2);
gPad->SetLogy(0);
h_R_pos_noAvgCorr->Draw("hist,e,same");
h_R_pos_withAvgCorr->Draw("hist,e,same");
can->cd(3);
gPad->SetLogy(0);
h_Z_pos_noAvgCorr->Draw("hist,e,same");
h_Z_pos_withAvgCorr->Draw("hist,e,same");
can->cd(4);
TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d, |Z|=%d cm",run,(int)selectZ));
legend->AddEntry(h_P_pos_noAvgCorr, Form("No Avg Corr"), "l");
legend->AddEntry(h_P_pos_withAvgCorr, Form("With Avg Corr"), "l");
legend->SetTextSize(0.1);
legend->Draw();
can->cd(5);
gPad->SetLogy(0);
h_P_neg_noAvgCorr->Draw("hist,e,same");
h_P_neg_withAvgCorr->Draw("hist,e,same");
can->cd(6);
gPad->SetLogy(0);
h_R_neg_noAvgCorr->Draw("hist,e,same");
h_R_neg_withAvgCorr->Draw("hist,e,same");
can->cd(7);
gPad->SetLogy(0);
h_Z_neg_noAvgCorr->Draw("hist,e,same");
h_Z_neg_withAvgCorr->Draw("hist,e,same");

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
