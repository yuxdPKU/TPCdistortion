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

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  cdbfilename[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  cout<<"CDB file path: "<<cdbfilename[i].c_str()<<endl;
  cdbfile[i] = new TFile(cdbfilename[i].c_str(),"");
  hcdb_R_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_posz");
  hcdb_P_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_posz");
  hcdb_Z_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_posz");
  hcdb_N_pr_pos[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");
  hcdb_R_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_negz");
  hcdb_P_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_negz");
  hcdb_Z_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_negz");
  hcdb_N_pr_neg[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");

  hcdb_N_pos[i] = new TH1F(Form("hentries_posz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_N_neg[i] = new TH1F(Form("hentries_negz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmax());
  plot1D(hcdb_N_pr_pos[i],hcdb_N_pos[i]);
  plot1D(hcdb_R_pr_pos[i],hcdb_R_pos[i]);
  plot1D(hcdb_P_pr_pos[i],hcdb_P_pos[i]);
  plot1D(hcdb_Z_pr_pos[i],hcdb_Z_pos[i]);
  plot1D(hcdb_N_pr_neg[i],hcdb_N_neg[i]);
  plot1D(hcdb_R_pr_neg[i],hcdb_R_neg[i]);
  plot1D(hcdb_P_pr_neg[i],hcdb_P_neg[i]);
  plot1D(hcdb_Z_pr_neg[i],hcdb_Z_neg[i]);
  hcdb_N_pos[i]->SetLineColor(1); hcdb_N_pos[i]->SetLineWidth(1); hcdb_N_pos[i]->SetFillColor(0); hcdb_N_pos[i]->SetMarkerColor(1);
  hcdb_R_pos[i]->SetLineColor(1); hcdb_R_pos[i]->SetLineWidth(1); hcdb_R_pos[i]->SetFillColor(0); hcdb_R_pos[i]->SetMarkerColor(1);
  hcdb_P_pos[i]->SetLineColor(1); hcdb_P_pos[i]->SetLineWidth(1); hcdb_P_pos[i]->SetFillColor(0); hcdb_P_pos[i]->SetMarkerColor(1);
  hcdb_Z_pos[i]->SetLineColor(1); hcdb_Z_pos[i]->SetLineWidth(1); hcdb_Z_pos[i]->SetFillColor(0); hcdb_Z_pos[i]->SetMarkerColor(1);
  hcdb_N_neg[i]->SetLineColor(1); hcdb_N_neg[i]->SetLineWidth(1); hcdb_N_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_R_neg[i]->SetLineColor(1); hcdb_R_neg[i]->SetLineWidth(1); hcdb_R_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_P_neg[i]->SetLineColor(1); hcdb_P_neg[i]->SetLineWidth(1); hcdb_P_neg[i]->SetFillColor(0); hcdb_P_neg[i]->SetMarkerColor(1);
  hcdb_Z_neg[i]->SetLineColor(1); hcdb_Z_neg[i]->SetLineWidth(1); hcdb_Z_neg[i]->SetFillColor(0); hcdb_Z_neg[i]->SetMarkerColor(1);
}

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun];
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

  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos[i],h_N_pos[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos[i],h_R_pos[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos[i],h_P_pos[i],selectZ, true);
  plot1D_Zbin(h_Z_prz_pos[i],h_Z_pos[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg[i],h_N_neg[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg[i],h_R_neg[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg[i],h_P_neg[i],-selectZ, true);
  plot1D_Zbin(h_Z_prz_neg[i],h_Z_neg[i],-selectZ, false);
  h_N_pos[i]->SetLineColor(2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(2);
  h_R_pos[i]->SetLineColor(2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(2);
  h_P_pos[i]->SetLineColor(2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(2);
  h_Z_pos[i]->SetLineColor(2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(2);
  h_N_neg[i]->SetLineColor(2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(2);
  h_R_neg[i]->SetLineColor(2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(2);
  h_P_neg[i]->SetLineColor(2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(2);
  h_Z_neg[i]->SetLineColor(2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(2);

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg[i], h_R_pos[i], hcdb_R_neg[i], hcdb_R_pos[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg[i], h_P_pos[i], hcdb_P_neg[i], hcdb_P_pos[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg[i], h_Z_pos[i], hcdb_Z_neg[i], hcdb_Z_pos[i]});
}

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos[i]->Write();
  hcdb_R_pos[i]->Write();
  h_R_neg[i]->Write();
  hcdb_R_neg[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos[i]->Write();
  hcdb_P_pos[i]->Write();
  h_P_neg[i]->Write();
  hcdb_P_neg[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos[i]->Write();
  hcdb_Z_pos[i]->Write();
  h_Z_neg[i]->Write();
  hcdb_Z_neg[i]->Write();
}
ofile3->Write();

for (int i=0; i<nrun; i++)
{
TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
h_P_pos[i]->Draw("hist,e,same");
hcdb_P_pos[i]->Draw("hist,same");
can->cd(2);
gPad->SetLogy(0);
h_R_pos[i]->Draw("hist,e,same");
hcdb_R_pos[i]->Draw("hist,same");
can->cd(3);
gPad->SetLogy(0);
h_Z_pos[i]->Draw("hist,e,same");
hcdb_Z_pos[i]->Draw("hist,same");
can->cd(4);
gPad->SetLogy(1);
h_N_pos[i]->Draw("hist,e,same");
//hcdb_N_pos[i]->Draw("hist,same");
can->cd(5);
gPad->SetLogy(0);
h_P_neg[i]->Draw("hist,e,same");
hcdb_P_neg[i]->Draw("hist,same");
can->cd(6);
gPad->SetLogy(0);
h_R_neg[i]->Draw("hist,e,same");
hcdb_R_neg[i]->Draw("hist,same");
can->cd(7);
gPad->SetLogy(0);
h_Z_neg[i]->Draw("hist,e,same");
hcdb_Z_neg[i]->Draw("hist,same");
can->cd(8);
gPad->SetLogy(1);
h_N_neg[i]->Draw("hist,e,same");
//hcdb_N_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d.pdf",(int)selectZ,runs[i]));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->AddEntry(h_P_pos[i], Form("Run %d",runs[i]), "l");
legend->AddEntry(hcdb_P_pos[i], Form("Run %d, Lamination",runs[i]), "l");
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
