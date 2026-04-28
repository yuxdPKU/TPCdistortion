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

const int nrun = 2;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nrun] = {79507,79516};

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

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

  h_N_pos_method1[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_N_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_N_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_R_pos_method1[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_P_pos_method1[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_Z_pos_method1[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_Z_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_Z_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_N_neg_method1[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_N_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_N_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_R_neg_method1[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_P_neg_method1[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_Z_neg_method1[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_Z_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_Z_prz_neg_method1[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method1[i],h_N_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos_method1[i],h_R_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos_method1[i],h_P_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_Z_prz_pos_method1[i],h_Z_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg_method1[i],h_N_neg_method1[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg_method1[i],h_R_neg_method1[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg_method1[i],h_P_neg_method1[i],-selectZ, false);
  plot1D_Zbin(h_Z_prz_neg_method1[i],h_Z_neg_method1[i],-selectZ, false);
  h_N_pos_method1[i]->SetLineColor(2); h_N_pos_method1[i]->SetLineWidth(1); h_N_pos_method1[i]->SetFillColor(0); h_N_pos_method1[i]->SetMarkerColor(2);
  h_R_pos_method1[i]->SetLineColor(2); h_R_pos_method1[i]->SetLineWidth(1); h_R_pos_method1[i]->SetFillColor(0); h_R_pos_method1[i]->SetMarkerColor(2);
  h_P_pos_method1[i]->SetLineColor(2); h_P_pos_method1[i]->SetLineWidth(1); h_P_pos_method1[i]->SetFillColor(0); h_P_pos_method1[i]->SetMarkerColor(2);
  h_Z_pos_method1[i]->SetLineColor(2); h_Z_pos_method1[i]->SetLineWidth(1); h_Z_pos_method1[i]->SetFillColor(0); h_Z_pos_method1[i]->SetMarkerColor(2);
  h_N_neg_method1[i]->SetLineColor(2); h_N_neg_method1[i]->SetLineWidth(1); h_N_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_R_neg_method1[i]->SetLineColor(2); h_R_neg_method1[i]->SetLineWidth(1); h_R_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_P_neg_method1[i]->SetLineColor(2); h_P_neg_method1[i]->SetLineWidth(1); h_P_neg_method1[i]->SetFillColor(0); h_P_neg_method1[i]->SetMarkerColor(2);
  h_Z_neg_method1[i]->SetLineColor(2); h_Z_neg_method1[i]->SetLineWidth(1); h_Z_neg_method1[i]->SetFillColor(0); h_Z_neg_method1[i]->SetMarkerColor(2);

  file_3D_map_method2[i] = new TFile(Form("../run3pp_newAlignment_closuretest/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_posz");
  h_R_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_negz");

  h_N_pos_method2[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_N_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_N_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_R_pos_method2[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_P_pos_method2[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_Z_pos_method2[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_Z_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_Z_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_N_neg_method2[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_N_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_N_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_R_neg_method2[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_P_neg_method2[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_Z_neg_method2[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_Z_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_Z_prz_neg_method2[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method2[i],h_N_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos_method2[i],h_R_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos_method2[i],h_P_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_Z_prz_pos_method2[i],h_Z_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg_method2[i],h_N_neg_method2[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg_method2[i],h_R_neg_method2[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg_method2[i],h_P_neg_method2[i],-selectZ, false);
  plot1D_Zbin(h_Z_prz_neg_method2[i],h_Z_neg_method2[i],-selectZ, false);
  h_N_pos_method2[i]->SetLineColor((4)); h_N_pos_method2[i]->SetLineWidth(1); h_N_pos_method2[i]->SetFillColor(0); h_N_pos_method2[i]->SetMarkerColor((4));
  h_R_pos_method2[i]->SetLineColor((4)); h_R_pos_method2[i]->SetLineWidth(1); h_R_pos_method2[i]->SetFillColor(0); h_R_pos_method2[i]->SetMarkerColor((4));
  h_P_pos_method2[i]->SetLineColor((4)); h_P_pos_method2[i]->SetLineWidth(1); h_P_pos_method2[i]->SetFillColor(0); h_P_pos_method2[i]->SetMarkerColor((4));
  h_Z_pos_method2[i]->SetLineColor((4)); h_Z_pos_method2[i]->SetLineWidth(1); h_Z_pos_method2[i]->SetFillColor(0); h_Z_pos_method2[i]->SetMarkerColor((4));
  h_N_neg_method2[i]->SetLineColor((4)); h_N_neg_method2[i]->SetLineWidth(1); h_N_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor((4));
  h_R_neg_method2[i]->SetLineColor((4)); h_R_neg_method2[i]->SetLineWidth(1); h_R_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor((4));
  h_P_neg_method2[i]->SetLineColor((4)); h_P_neg_method2[i]->SetLineWidth(1); h_P_neg_method2[i]->SetFillColor(0); h_P_neg_method2[i]->SetMarkerColor((4));
  h_Z_neg_method2[i]->SetLineColor((4)); h_Z_neg_method2[i]->SetLineWidth(1); h_Z_neg_method2[i]->SetFillColor(0); h_Z_neg_method2[i]->SetMarkerColor((4));

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_method1[i], h_R_pos_method1[i], h_R_neg_method2[i], h_R_pos_method2[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_method1[i], h_P_pos_method1[i], h_P_neg_method2[i], h_P_pos_method2[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method2[i], h_Z_pos_method2[i]});
}

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
}
ofile3->Write();

for (int i=0; i<nrun; i++)
{
TCanvas* can = new TCanvas("can","",2400,1200);
can->Divide(3,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->Draw("hist,e,same");
h_P_pos_method2[i]->Draw("hist,e,same");
can->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->Draw("hist,e,same");
h_R_pos_method2[i]->Draw("hist,e,same");
can->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
h_Z_pos_method2[i]->Draw("hist,e,same");
can->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->Draw("hist,e,same");
h_P_neg_method2[i]->Draw("hist,e,same");
can->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->Draw("hist,e,same");
h_R_neg_method2[i]->Draw("hist,e,same");
can->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
h_Z_neg_method2[i]->Draw("hist,e,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d.pdf",(int)selectZ,runs[i]));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d, closure test",runs[i]));
legend->AddEntry(h_P_pos_method1[i], Form("w/o TPOT-based correction"), "l");
legend->AddEntry(h_P_pos_method2[i], Form("w/ TPOT-based correction"), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d.pdf",runs[i]));

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
