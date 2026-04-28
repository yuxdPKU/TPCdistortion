#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);
TH1* SubtractHistograms(const TH1* h1, const TH1* h2, const char* name = "diff");
void draw1Dmap(TString filename, TString tag, double selectZ, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color=1, bool convert_RP_2_P=false);
void draw1Dmap_2D(TString filename, TString tag, TH1*& h_R_pos, TH1*& h_R_neg, TH1*& h_P_pos, TH1*& h_P_neg, TH1*& h_Z_pos, TH1*& h_Z_neg, int color);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, bool convert_RP_2_P=false)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      double value = h3->GetBinContent(xbin, i, zbin);
      double error = h3->GetBinError(xbin, i, zbin);
      if (convert_RP_2_P)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
        value /= rcenter;
        error /= rcenter;
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

const int nrun = 1;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nrun] = {79516};

//get map from cdb lamination study
std::string cdbfilename_lamination[nrun];
TH1 *hcdb_lamination_N_pos[nrun], *hcdb_lamination_R_pos[nrun], *hcdb_lamination_P_pos[nrun], *hcdb_lamination_Z_pos[nrun];
TH1 *hcdb_lamination_N_neg[nrun], *hcdb_lamination_R_neg[nrun], *hcdb_lamination_P_neg[nrun], *hcdb_lamination_Z_neg[nrun];

//get map from cdb moduleedge study
std::string cdbfilename_moduleedge[nrun];
TH1 *hcdb_moduleedge_N_pos[nrun], *hcdb_moduleedge_R_pos[nrun], *hcdb_moduleedge_P_pos[nrun], *hcdb_moduleedge_Z_pos[nrun];
TH1 *hcdb_moduleedge_N_neg[nrun], *hcdb_moduleedge_R_neg[nrun], *hcdb_moduleedge_P_neg[nrun], *hcdb_moduleedge_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  cdbfilename_lamination[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  cout<<"Lamination CDB file path: "<<cdbfilename_lamination[i].c_str()<<endl;

  cdbfilename_moduleedge[i] = CDBInterface::instance()->getUrl("TPC_Module_Edge");
  cout<<"Module-edge CDB file path: "<<cdbfilename_moduleedge[i].c_str()<<endl;

  draw1Dmap_2D(cdbfilename_lamination[i].c_str(), Form("%d_lamination_cdb",runs[i]), hcdb_lamination_R_pos[i], hcdb_lamination_R_neg[i], hcdb_lamination_P_pos[i], hcdb_lamination_P_neg[i], hcdb_lamination_Z_pos[i], hcdb_lamination_Z_neg[i], 1);

  draw1Dmap_2D(cdbfilename_moduleedge[i].c_str(), Form("%d_moduleedge_cdb",runs[i]), hcdb_moduleedge_R_pos[i], hcdb_moduleedge_R_neg[i], hcdb_moduleedge_P_pos[i], hcdb_moduleedge_P_neg[i], hcdb_moduleedge_Z_pos[i], hcdb_moduleedge_Z_neg[i], 1);
}

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
for (const auto& selectZ : selectZs)
{

TH1 *h_N_pos_method1[nrun], *h_R_pos_method1[nrun], *h_P_pos_method1[nrun], *h_Z_pos_method1[nrun];
TH1 *h_N_neg_method1[nrun], *h_R_neg_method1[nrun], *h_P_neg_method1[nrun], *h_Z_neg_method1[nrun];

TH1 *h_N_pos_method2[nrun], *h_R_pos_method2[nrun], *h_P_pos_method2[nrun], *h_Z_pos_method2[nrun];
TH1 *h_N_neg_method2[nrun], *h_R_neg_method2[nrun], *h_P_neg_method2[nrun], *h_Z_neg_method2[nrun];

TH1 *h_N_pos_method3[nrun], *h_R_pos_method3[nrun], *h_P_pos_method3[nrun], *h_Z_pos_method3[nrun];
TH1 *h_N_neg_method3[nrun], *h_R_neg_method3[nrun], *h_P_neg_method3[nrun], *h_Z_neg_method3[nrun];

TH1 *h_N_pos_method4[nrun], *h_R_pos_method4[nrun], *h_P_pos_method4[nrun], *h_Z_pos_method4[nrun];
TH1 *h_N_neg_method4[nrun], *h_R_neg_method4[nrun], *h_P_neg_method4[nrun], *h_Z_neg_method4[nrun];

TH1 *h_N_pos_static_moduleedge[nrun], *h_R_pos_static_moduleedge[nrun], *h_P_pos_static_moduleedge[nrun], *h_Z_pos_static_moduleedge[nrun];
TH1 *h_N_neg_static_moduleedge[nrun], *h_R_neg_static_moduleedge[nrun], *h_P_neg_static_moduleedge[nrun], *h_Z_neg_static_moduleedge[nrun];

TH1 *h_N_pos_static[nrun], *h_R_pos_static[nrun], *h_P_pos_static[nrun], *h_Z_pos_static[nrun];
TH1 *h_N_neg_static[nrun], *h_R_neg_static[nrun], *h_P_neg_static[nrun], *h_Z_neg_static[nrun];

TH1 *h_N_pos_moduleedge[nrun], *h_R_pos_moduleedge[nrun], *h_P_pos_moduleedge[nrun], *h_Z_pos_moduleedge[nrun];
TH1 *h_N_neg_moduleedge[nrun], *h_R_neg_moduleedge[nrun], *h_P_neg_moduleedge[nrun], *h_Z_neg_moduleedge[nrun];

//get map from cdb static study
std::string cdbfilename_static[nrun];
TH1 *hcdb_static_N_pos[nrun], *hcdb_static_R_pos[nrun], *hcdb_static_P_pos[nrun], *hcdb_static_Z_pos[nrun];
TH1 *hcdb_static_N_neg[nrun], *hcdb_static_R_neg[nrun], *hcdb_static_P_neg[nrun], *hcdb_static_Z_neg[nrun];

TH1 *h_DC_N_pos[nrun], *h_DC_R_pos[nrun], *h_DC_P_pos[nrun], *h_DC_Z_pos[nrun];
TH1 *h_DC_N_neg[nrun], *h_DC_R_neg[nrun], *h_DC_P_neg[nrun], *h_DC_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{

  draw1Dmap(Form("../run3pp_newSiFieldonAlignment_newTPOTzfAlignment/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_allcorrectionapplied",runs[i]), selectZ, h_R_pos_method1[i], h_R_neg_method1[i], h_P_pos_method1[i], h_P_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method1[i], 2);

  draw1Dmap(Form("../run3pp_newSiFieldonAlignment_newTPOTzfAlignment_disableOtherCorr/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_nocorrectionapplied",runs[i]), selectZ, h_R_pos_method2[i], h_R_neg_method2[i], h_P_pos_method2[i], h_P_neg_method2[i], h_Z_pos_method2[i], h_Z_neg_method2[i], 4);

  draw1Dmap(Form("../run3pp_newSiFieldonAlignment_newTPOTzfAlignment_disableModuleEdgeCorr/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_nocorrectionapplied",runs[i]), selectZ, h_R_pos_method3[i], h_R_neg_method3[i], h_P_pos_method3[i], h_P_neg_method3[i], h_Z_pos_method3[i], h_Z_neg_method3[i], 6);

  draw1Dmap(Form("../run3pp_newSiFieldonAlignment_newTPOTzfAlignment_disableStaticCorr/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]), Form("%d_nocorrectionapplied",runs[i]), selectZ, h_R_pos_method4[i], h_R_neg_method4[i], h_P_pos_method4[i], h_P_neg_method4[i], h_Z_pos_method4[i], h_Z_neg_method4[i], 7);

  cdbfilename_static[i] = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/static_only_inverted_10-new.root";
  cout<<"Static CDB file path: "<<cdbfilename_static[i].c_str()<<endl;

  draw1Dmap(cdbfilename_static[i].c_str(), Form("%d_static_cdb",runs[i]), selectZ, hcdb_static_R_pos[i], hcdb_static_R_neg[i], hcdb_static_P_pos[i], hcdb_static_P_neg[i], hcdb_static_Z_pos[i], hcdb_static_Z_neg[i], 1, true);

  draw1Dmap(Form("/sphenix/u/hangal/Public/forXudong/dc_distortions_run%d/_0p9.distortion_map.hist.root",runs[i]), Form("%d_nocorrectionapplied",runs[i]), selectZ, h_DC_R_pos[i], h_DC_R_neg[i], h_DC_P_pos[i], h_DC_P_neg[i], h_DC_Z_pos[i], h_DC_Z_neg[i], 1);

  h_R_pos_static_moduleedge[i] = SubtractHistograms(h_R_pos_method2[i], h_R_pos_method1[i], "R_pos_static_moduleedge");
  h_R_neg_static_moduleedge[i] = SubtractHistograms(h_R_neg_method2[i], h_R_neg_method1[i], "R_neg_static_moduleedge");
  h_P_pos_static_moduleedge[i] = SubtractHistograms(h_P_pos_method2[i], h_P_pos_method1[i], "P_pos_static_moduleedge");
  h_P_neg_static_moduleedge[i] = SubtractHistograms(h_P_neg_method2[i], h_P_neg_method1[i], "P_neg_static_moduleedge");
  h_Z_pos_static_moduleedge[i] = SubtractHistograms(h_Z_pos_method2[i], h_Z_pos_method1[i], "Z_pos_static_moduleedge");
  h_Z_neg_static_moduleedge[i] = SubtractHistograms(h_Z_neg_method2[i], h_Z_neg_method1[i], "Z_neg_static_moduleedge");

  h_R_pos_static_moduleedge[i]->SetLineColor(6); h_R_pos_static_moduleedge[i]->SetLineWidth(1); h_R_pos_static_moduleedge[i]->SetFillColor(0); h_R_pos_static_moduleedge[i]->SetMarkerColor(6);
  h_P_pos_static_moduleedge[i]->SetLineColor(6); h_P_pos_static_moduleedge[i]->SetLineWidth(1); h_P_pos_static_moduleedge[i]->SetFillColor(0); h_P_pos_static_moduleedge[i]->SetMarkerColor(6);
  h_Z_pos_static_moduleedge[i]->SetLineColor(6); h_Z_pos_static_moduleedge[i]->SetLineWidth(1); h_Z_pos_static_moduleedge[i]->SetFillColor(0); h_Z_pos_static_moduleedge[i]->SetMarkerColor(6);
  h_R_neg_static_moduleedge[i]->SetLineColor(6); h_R_neg_static_moduleedge[i]->SetLineWidth(1); h_R_neg_static_moduleedge[i]->SetFillColor(0); h_R_neg_static_moduleedge[i]->SetMarkerColor(6);
  h_P_neg_static_moduleedge[i]->SetLineColor(6); h_P_neg_static_moduleedge[i]->SetLineWidth(1); h_P_neg_static_moduleedge[i]->SetFillColor(0); h_P_neg_static_moduleedge[i]->SetMarkerColor(6);
  h_Z_neg_static_moduleedge[i]->SetLineColor(6); h_Z_neg_static_moduleedge[i]->SetLineWidth(1); h_Z_neg_static_moduleedge[i]->SetFillColor(0); h_Z_neg_static_moduleedge[i]->SetMarkerColor(6);

  h_R_pos_static[i] = SubtractHistograms(h_R_pos_method4[i], h_R_pos_method1[i], "R_pos_static");
  h_R_neg_static[i] = SubtractHistograms(h_R_neg_method4[i], h_R_neg_method1[i], "R_neg_static");
  h_P_pos_static[i] = SubtractHistograms(h_P_pos_method4[i], h_P_pos_method1[i], "P_pos_static");
  h_P_neg_static[i] = SubtractHistograms(h_P_neg_method4[i], h_P_neg_method1[i], "P_neg_static");
  h_Z_pos_static[i] = SubtractHistograms(h_Z_pos_method4[i], h_Z_pos_method1[i], "Z_pos_static");
  h_Z_neg_static[i] = SubtractHistograms(h_Z_neg_method4[i], h_Z_neg_method1[i], "Z_neg_static");

  h_R_pos_static[i]->SetLineColor(6); h_R_pos_static[i]->SetLineWidth(1); h_R_pos_static[i]->SetFillColor(0); h_R_pos_static[i]->SetMarkerColor(6);
  h_P_pos_static[i]->SetLineColor(6); h_P_pos_static[i]->SetLineWidth(1); h_P_pos_static[i]->SetFillColor(0); h_P_pos_static[i]->SetMarkerColor(6);
  h_Z_pos_static[i]->SetLineColor(6); h_Z_pos_static[i]->SetLineWidth(1); h_Z_pos_static[i]->SetFillColor(0); h_Z_pos_static[i]->SetMarkerColor(6);
  h_R_neg_static[i]->SetLineColor(6); h_R_neg_static[i]->SetLineWidth(1); h_R_neg_static[i]->SetFillColor(0); h_R_neg_static[i]->SetMarkerColor(6);
  h_P_neg_static[i]->SetLineColor(6); h_P_neg_static[i]->SetLineWidth(1); h_P_neg_static[i]->SetFillColor(0); h_P_neg_static[i]->SetMarkerColor(6);
  h_Z_neg_static[i]->SetLineColor(6); h_Z_neg_static[i]->SetLineWidth(1); h_Z_neg_static[i]->SetFillColor(0); h_Z_neg_static[i]->SetMarkerColor(6);

  h_R_pos_moduleedge[i] = SubtractHistograms(h_R_pos_method3[i], h_R_pos_method1[i], "R_pos_moduleedge");
  h_R_neg_moduleedge[i] = SubtractHistograms(h_R_neg_method3[i], h_R_neg_method1[i], "R_neg_moduleedge");
  h_P_pos_moduleedge[i] = SubtractHistograms(h_P_pos_method3[i], h_P_pos_method1[i], "P_pos_moduleedge");
  h_P_neg_moduleedge[i] = SubtractHistograms(h_P_neg_method3[i], h_P_neg_method1[i], "P_neg_moduleedge");
  h_Z_pos_moduleedge[i] = SubtractHistograms(h_Z_pos_method3[i], h_Z_pos_method1[i], "Z_pos_moduleedge");
  h_Z_neg_moduleedge[i] = SubtractHistograms(h_Z_neg_method3[i], h_Z_neg_method1[i], "Z_neg_moduleedge");

  h_R_pos_moduleedge[i]->SetLineColor(6); h_R_pos_moduleedge[i]->SetLineWidth(1); h_R_pos_moduleedge[i]->SetFillColor(0); h_R_pos_moduleedge[i]->SetMarkerColor(6);
  h_P_pos_moduleedge[i]->SetLineColor(6); h_P_pos_moduleedge[i]->SetLineWidth(1); h_P_pos_moduleedge[i]->SetFillColor(0); h_P_pos_moduleedge[i]->SetMarkerColor(6);
  h_Z_pos_moduleedge[i]->SetLineColor(6); h_Z_pos_moduleedge[i]->SetLineWidth(1); h_Z_pos_moduleedge[i]->SetFillColor(0); h_Z_pos_moduleedge[i]->SetMarkerColor(6);
  h_R_neg_moduleedge[i]->SetLineColor(6); h_R_neg_moduleedge[i]->SetLineWidth(1); h_R_neg_moduleedge[i]->SetFillColor(0); h_R_neg_moduleedge[i]->SetMarkerColor(6);
  h_P_neg_moduleedge[i]->SetLineColor(6); h_P_neg_moduleedge[i]->SetLineWidth(1); h_P_neg_moduleedge[i]->SetFillColor(0); h_P_neg_moduleedge[i]->SetMarkerColor(6);
  h_Z_neg_moduleedge[i]->SetLineColor(6); h_Z_neg_moduleedge[i]->SetLineWidth(1); h_Z_neg_moduleedge[i]->SetFillColor(0); h_Z_neg_moduleedge[i]->SetMarkerColor(6);
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
	hists_P.push_back(h_P_neg_method4[i]);
	hists_P.push_back(h_P_pos_method4[i]);
	hists_P.push_back(h_P_neg_static_moduleedge[i]);
	hists_P.push_back(h_P_pos_static_moduleedge[i]);
	hists_P.push_back(h_P_neg_static[i]);
	hists_P.push_back(h_P_pos_static[i]);
	hists_P.push_back(h_P_neg_moduleedge[i]);
	hists_P.push_back(h_P_pos_moduleedge[i]);
	//hists_P.push_back(hcdb_lamination_P_pos[i]);
	//hists_P.push_back(hcdb_lamination_P_neg[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_R.push_back(h_R_neg_method1[i]);
	hists_R.push_back(h_R_pos_method1[i]);
	hists_R.push_back(h_R_neg_method2[i]);
	hists_R.push_back(h_R_pos_method2[i]);
	hists_R.push_back(h_R_neg_method3[i]);
	hists_R.push_back(h_R_pos_method3[i]);
	hists_R.push_back(h_R_neg_method4[i]);
	hists_R.push_back(h_R_pos_method4[i]);
	hists_R.push_back(h_R_neg_static_moduleedge[i]);
	hists_R.push_back(h_R_pos_static_moduleedge[i]);
	hists_R.push_back(h_R_neg_static[i]);
	hists_R.push_back(h_R_pos_static[i]);
	hists_R.push_back(h_R_neg_moduleedge[i]);
	hists_R.push_back(h_R_pos_moduleedge[i]);
	//hists_R.push_back(hcdb_lamination_R_pos[i]);
	//hists_R.push_back(hcdb_lamination_R_neg[i]);
}
for (int i=0; i<nrun; i++)
{
	hists_Z.push_back(h_Z_neg_method1[i]);
	hists_Z.push_back(h_Z_pos_method1[i]);
	hists_Z.push_back(h_Z_neg_method2[i]);
	hists_Z.push_back(h_Z_pos_method2[i]);
	hists_Z.push_back(h_Z_neg_method3[i]);
	hists_Z.push_back(h_Z_pos_method3[i]);
	hists_Z.push_back(h_Z_neg_method4[i]);
	hists_Z.push_back(h_Z_pos_method4[i]);
	hists_Z.push_back(h_Z_neg_static_moduleedge[i]);
	hists_Z.push_back(h_Z_pos_static_moduleedge[i]);
	hists_Z.push_back(h_Z_neg_static[i]);
	hists_Z.push_back(h_Z_pos_static[i]);
	hists_Z.push_back(h_Z_neg_moduleedge[i]);
	hists_Z.push_back(h_Z_pos_moduleedge[i]);
	//hists_Z.push_back(hcdb_lamination_Z_pos[i]);
	//hists_Z.push_back(hcdb_lamination_Z_neg[i]);
}
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

TFile* ofile = new TFile(Form("Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_pos_method3[i]->Write();
  h_R_pos_method4[i]->Write();
  h_R_pos_static_moduleedge[i]->Write();
  h_R_pos_static[i]->Write();
  h_R_pos_moduleedge[i]->Write();
  hcdb_lamination_R_pos[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
  h_R_neg_method3[i]->Write();
  h_R_neg_method4[i]->Write();
  h_R_neg_static_moduleedge[i]->Write();
  h_R_neg_static[i]->Write();
  h_R_neg_moduleedge[i]->Write();
  hcdb_lamination_R_neg[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("Rootfiles/hist_dphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_pos_method3[i]->Write();
  h_P_pos_method4[i]->Write();
  h_P_pos_static_moduleedge[i]->Write();
  h_P_pos_static[i]->Write();
  h_P_pos_moduleedge[i]->Write();
  hcdb_lamination_P_pos[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
  h_P_neg_method3[i]->Write();
  h_P_neg_method4[i]->Write();
  h_P_neg_static_moduleedge[i]->Write();
  h_P_neg_static[i]->Write();
  h_P_neg_moduleedge[i]->Write();
  hcdb_lamination_P_neg[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_pos_method3[i]->Write();
  h_Z_pos_method4[i]->Write();
  h_Z_pos_static_moduleedge[i]->Write();
  h_Z_pos_static[i]->Write();
  h_Z_pos_moduleedge[i]->Write();
  hcdb_lamination_Z_pos[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
  h_Z_neg_method3[i]->Write();
  h_Z_neg_method4[i]->Write();
  h_Z_neg_static_moduleedge[i]->Write();
  h_Z_neg_static[i]->Write();
  h_Z_neg_moduleedge[i]->Write();
  hcdb_lamination_Z_neg[i]->Write();
}
ofile3->Write();

for (int i=0; i<nrun; i++)
{
TCanvas* can_TPOT = new TCanvas("can_TPOT","",2400,1200);
can_TPOT->Divide(3,2);
can_TPOT->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->Draw("hist,e,same");
h_P_pos_method2[i]->Draw("hist,e,same");
h_P_pos_method3[i]->Draw("hist,e,same");
h_P_pos_method4[i]->Draw("hist,e,same");
can_TPOT->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->Draw("hist,e,same");
h_R_pos_method2[i]->Draw("hist,e,same");
h_R_pos_method3[i]->Draw("hist,e,same");
h_R_pos_method4[i]->Draw("hist,e,same");
can_TPOT->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
h_Z_pos_method2[i]->Draw("hist,e,same");
h_Z_pos_method3[i]->Draw("hist,e,same");
h_Z_pos_method4[i]->Draw("hist,e,same");
can_TPOT->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->Draw("hist,e,same");
h_P_neg_method2[i]->Draw("hist,e,same");
h_P_neg_method3[i]->Draw("hist,e,same");
h_P_neg_method4[i]->Draw("hist,e,same");
can_TPOT->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->Draw("hist,e,same");
h_R_neg_method2[i]->Draw("hist,e,same");
h_R_neg_method3[i]->Draw("hist,e,same");
h_R_neg_method4[i]->Draw("hist,e,same");
can_TPOT->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
h_Z_neg_method2[i]->Draw("hist,e,same");
h_Z_neg_method3[i]->Draw("hist,e,same");
h_Z_neg_method4[i]->Draw("hist,e,same");

gPad->RedrawAxis();

can_TPOT->Update();
can_TPOT->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_TPOT.pdf",(int)selectZ,runs[i]));

TLegend *legend_TPOT = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_TPOT->SetHeader(Form("Run %d",runs[i]));
legend_TPOT->AddEntry(h_P_pos_method1[i], Form("Module-edge & Static on -- SC"), "l");
legend_TPOT->AddEntry(h_P_pos_method2[i], Form("Module-edge & Static off -- ME + STATIC + SC"), "l");
legend_TPOT->AddEntry(h_P_pos_method3[i], Form("Module-edge off, Static on -- ME + SC"), "l");
legend_TPOT->AddEntry(h_P_pos_method4[i], Form("Module-edge on, Static off -- STATIC + SC"), "l");
legend_TPOT->Draw();
TCanvas* can_TPOT_leg = new TCanvas("can_TPOT_leg","",2500,1000);
legend_TPOT->Draw();
can_TPOT_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_TPOT.pdf",runs[i]));

delete can_TPOT;
delete can_TPOT_leg;

TCanvas* can_TPOT_vs_Lamination = new TCanvas("can_TPOT_vs_Lamination","",2400,1200);
can_TPOT_vs_Lamination->Divide(3,2);
can_TPOT_vs_Lamination->cd(1);
gPad->SetLogy(0);
h_P_pos_method1[i]->SetMinimum(-7e-3);
h_P_pos_method1[i]->Draw("hist,e,same");
hcdb_lamination_P_pos[i]->Draw("hist,same");
can_TPOT_vs_Lamination->cd(2);
gPad->SetLogy(0);
h_R_pos_method1[i]->SetMinimum(-1);
h_R_pos_method1[i]->SetMaximum(0.5);
h_R_pos_method1[i]->Draw("hist,e,same");
hcdb_lamination_R_pos[i]->Draw("hist,same");
can_TPOT_vs_Lamination->cd(3);
gPad->SetLogy(0);
h_Z_pos_method1[i]->Draw("hist,e,same");
hcdb_lamination_Z_pos[i]->Draw("hist,same");
can_TPOT_vs_Lamination->cd(4);
gPad->SetLogy(0);
h_P_neg_method1[i]->SetMinimum(-7e-3);
h_P_neg_method1[i]->Draw("hist,e,same");
hcdb_lamination_P_neg[i]->Draw("hist,same");
can_TPOT_vs_Lamination->cd(5);
gPad->SetLogy(0);
h_R_neg_method1[i]->SetMinimum(-1);
h_R_neg_method1[i]->SetMaximum(0.5);
h_R_neg_method1[i]->Draw("hist,e,same");
hcdb_lamination_R_neg[i]->Draw("hist,same");
can_TPOT_vs_Lamination->cd(6);
gPad->SetLogy(0);
h_Z_neg_method1[i]->Draw("hist,e,same");
hcdb_lamination_Z_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can_TPOT_vs_Lamination->Update();
can_TPOT_vs_Lamination->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_TPOT_vs_Lamination.pdf",(int)selectZ,runs[i]));

TLegend *legend_TPOT_vs_Lamination = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_TPOT_vs_Lamination->SetHeader(Form("Run %d",runs[i]));
legend_TPOT_vs_Lamination->AddEntry(h_P_pos_method1[i], Form("Module-edge & Static on -- SC"), "l");
legend_TPOT_vs_Lamination->AddEntry(hcdb_lamination_P_pos[i], Form("Lamination"), "l");
legend_TPOT_vs_Lamination->Draw();
TCanvas* can_TPOT_vs_Lamination_leg = new TCanvas("can_TPOT_vs_Lamination_leg","",2500,1000);
legend_TPOT_vs_Lamination->Draw();
can_TPOT_vs_Lamination_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_TPOT_vs_Lamination.pdf",runs[i]));

delete can_TPOT_vs_Lamination;
delete can_TPOT_vs_Lamination_leg;

TCanvas* can_TPOT_moduleedge = new TCanvas("can_TPOT_moduleedge","",2400,1200);
can_TPOT_moduleedge->Divide(3,2);
can_TPOT_moduleedge->cd(1);
gPad->SetLogy(0);
h_P_pos_moduleedge[i]->SetMinimum(-5e-3);
h_P_pos_moduleedge[i]->SetMaximum(5e-3);
h_P_pos_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_P_pos[i]->Draw("hist,same");
can_TPOT_moduleedge->cd(2);
gPad->SetLogy(0);
h_R_pos_moduleedge[i]->SetMinimum(-0.01);
h_R_pos_moduleedge[i]->SetMaximum(0.01);
h_R_pos_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_R_pos[i]->Draw("hist,same");
can_TPOT_moduleedge->cd(3);
gPad->SetLogy(0);
h_Z_pos_moduleedge[i]->SetMinimum(-0.01);
h_Z_pos_moduleedge[i]->SetMaximum(0.01);
h_Z_pos_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_Z_pos[i]->Draw("hist,same");
can_TPOT_moduleedge->cd(4);
gPad->SetLogy(0);
h_P_neg_moduleedge[i]->SetMinimum(-5e-3);
h_P_neg_moduleedge[i]->SetMaximum(5e-3);
h_P_neg_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_P_neg[i]->Draw("hist,same");
can_TPOT_moduleedge->cd(5);
gPad->SetLogy(0);
h_R_neg_moduleedge[i]->SetMinimum(-0.01);
h_R_neg_moduleedge[i]->SetMaximum(0.01);
h_R_neg_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_R_neg[i]->Draw("hist,same");
can_TPOT_moduleedge->cd(6);
gPad->SetLogy(0);
h_Z_neg_moduleedge[i]->SetMinimum(-0.01);
h_Z_neg_moduleedge[i]->SetMaximum(0.01);
h_Z_neg_moduleedge[i]->Draw("hist,e,same");
hcdb_moduleedge_Z_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can_TPOT_moduleedge->Update();
can_TPOT_moduleedge->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_TPOT_vs_ME.pdf",(int)selectZ,runs[i]));

TLegend *legend_TPOT_moduleedge = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_TPOT_moduleedge->SetHeader(Form("Run %d",runs[i]));
legend_TPOT_moduleedge->AddEntry(h_P_pos_moduleedge[i], Form("TPOT -- ME"), "l");
legend_TPOT_moduleedge->AddEntry(hcdb_moduleedge_P_pos[i], Form("CDB -- ME"), "l");
legend_TPOT_moduleedge->Draw();
TCanvas* can_TPOT_moduleedge_leg = new TCanvas("can_TPOT_moduleedge_leg","",2500,1000);
legend_TPOT_moduleedge->Draw();
can_TPOT_moduleedge_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_TPOT_vs_ME.pdf",runs[i]));

delete can_TPOT_moduleedge;
delete can_TPOT_moduleedge_leg;

TCanvas* can_TPOT_static = new TCanvas("can_TPOT_static","",2400,1200);
can_TPOT_static->Divide(3,2);
can_TPOT_static->cd(1);
gPad->SetLogy(0);
h_P_pos_static[i]->SetMinimum(-70e-3);
h_P_pos_static[i]->SetMaximum(0e-3);
h_P_pos_static[i]->Draw("hist,e,same");
hcdb_static_P_pos[i]->Draw("hist,same");
can_TPOT_static->cd(2);
gPad->SetLogy(0);
h_R_pos_static[i]->SetMinimum(0);
h_R_pos_static[i]->SetMaximum(3);
h_R_pos_static[i]->Draw("hist,e,same");
hcdb_static_R_pos[i]->Draw("hist,same");
can_TPOT_static->cd(3);
gPad->SetLogy(0);
h_Z_pos_static[i]->SetMinimum(-0.4);
h_Z_pos_static[i]->SetMaximum(0.2);
h_Z_pos_static[i]->Draw("hist,e,same");
hcdb_static_Z_pos[i]->Draw("hist,same");
can_TPOT_static->cd(4);
gPad->SetLogy(0);
h_P_neg_static[i]->SetMinimum(-70e-3);
h_P_neg_static[i]->SetMaximum(0e-3);
h_P_neg_static[i]->Draw("hist,e,same");
hcdb_static_P_neg[i]->Draw("hist,same");
can_TPOT_static->cd(5);
gPad->SetLogy(0);
h_R_neg_static[i]->SetMinimum(0);
h_R_neg_static[i]->SetMaximum(3);
h_R_neg_static[i]->Draw("hist,e,same");
hcdb_static_R_neg[i]->Draw("hist,same");
can_TPOT_static->cd(6);
gPad->SetLogy(0);
h_Z_neg_static[i]->SetMinimum(-0.4);
h_Z_neg_static[i]->SetMaximum(0.2);
h_Z_neg_static[i]->Draw("hist,e,same");
hcdb_static_Z_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can_TPOT_static->Update();
can_TPOT_static->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_TPOT_vs_static.pdf",(int)selectZ,runs[i]));

TLegend *legend_TPOT_static = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_TPOT_static->SetHeader(Form("Run %d",runs[i]));
legend_TPOT_static->AddEntry(h_P_pos_static[i], Form("TPOT -- static"), "l");
legend_TPOT_static->AddEntry(hcdb_static_P_pos[i], Form("CDB -- static"), "l");
legend_TPOT_static->Draw();
TCanvas* can_TPOT_static_leg = new TCanvas("can_TPOT_static_leg","",2500,1000);
legend_TPOT_static->Draw();
can_TPOT_static_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_TPOT_vs_static.pdf",runs[i]));

delete can_TPOT_static;
delete can_TPOT_static_leg;

TCanvas* can_TPOT_DC = new TCanvas("can_TPOT_DC","",2400,1200);
can_TPOT_DC->Divide(3,2);
can_TPOT_DC->cd(1);
gPad->SetLogy(0);
h_P_pos_method2[i]->SetMinimum(-20e-3);
h_P_pos_method2[i]->SetMaximum(20e-3);
h_P_pos_method2[i]->Draw("hist,e,same");
h_DC_P_pos[i]->Draw("hist,same");
can_TPOT_DC->cd(2);
gPad->SetLogy(0);
h_R_pos_method2[i]->SetMinimum(0);
h_R_pos_method2[i]->SetMaximum(3);
h_R_pos_method2[i]->Draw("hist,e,same");
h_DC_R_pos[i]->Draw("hist,same");
can_TPOT_DC->cd(3);
gPad->SetLogy(0);
h_Z_pos_method2[i]->SetMinimum(-0.4);
h_Z_pos_method2[i]->SetMaximum(0.3);
h_Z_pos_method2[i]->Draw("hist,e,same");
h_DC_Z_pos[i]->Draw("hist,same");
can_TPOT_DC->cd(4);
gPad->SetLogy(0);
h_P_neg_method2[i]->SetMinimum(-20e-3);
h_P_neg_method2[i]->SetMaximum(20e-3);
h_P_neg_method2[i]->Draw("hist,e,same");
h_DC_P_neg[i]->Draw("hist,same");
can_TPOT_DC->cd(5);
gPad->SetLogy(0);
h_R_neg_method2[i]->SetMinimum(0);
h_R_neg_method2[i]->SetMaximum(3);
h_R_neg_method2[i]->Draw("hist,e,same");
h_DC_R_neg[i]->Draw("hist,same");
can_TPOT_DC->cd(6);
gPad->SetLogy(0);
h_Z_neg_method2[i]->SetMinimum(-0.4);
h_Z_neg_method2[i]->SetMaximum(0.3);
h_Z_neg_method2[i]->Draw("hist,e,same");
h_DC_Z_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can_TPOT_DC->Update();
can_TPOT_DC->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%d_TPOT_vs_DC.pdf",(int)selectZ,runs[i]));

TLegend *legend_TPOT_DC = new TLegend(0.1, 0.1, 0.9, 0.9);
legend_TPOT_DC->SetHeader(Form("Run %d",runs[i]));
legend_TPOT_DC->AddEntry(h_P_pos_method2[i], Form("TPOT -- ME + static + SC"), "l");
legend_TPOT_DC->AddEntry(h_DC_P_pos[i], Form("DC -- static + SC"), "l");
legend_TPOT_DC->Draw();
TCanvas* can_TPOT_DC_leg = new TCanvas("can_TPOT_DC_leg","",2500,1000);
legend_TPOT_DC->Draw();
can_TPOT_DC_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%d_TPOT_vs_DC.pdf",runs[i]));

delete can_TPOT_DC;
delete can_TPOT_DC_leg;


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

TH1* SubtractHistograms(const TH1* h1, const TH1* h2, const char* name = "diff")
{
    if (!h1 || !h2) return nullptr;
    if (h1->GetNbinsX() != h2->GetNbinsX()) {
        std::cerr << "Error: histograms have different number of bins." << std::endl;
        return nullptr;
    }
    TH1* hdiff = (TH1*)h1->Clone(name);
    hdiff->Reset();
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double val1 = h1->GetBinContent(i);
        double err1 = h1->GetBinError(i);
        double val2 = h2->GetBinContent(i);
        double err2 = h2->GetBinError(i);
        double diff = val1 - val2;
        double err = sqrt(err1*err1 + err2*err2);
        hdiff->SetBinContent(i, diff);
        hdiff->SetBinError(i, err);
    }
    return hdiff;
}

void draw1Dmap(TString filename, TString tag, double selectZ,
  TH1*& h_R_pos, TH1*& h_R_neg,
  TH1*& h_P_pos, TH1*& h_P_neg,
  TH1*& h_Z_pos, TH1*& h_Z_neg,
  int color=1,
  bool convert_RP_2_P=false)
{
  TFile* file_3D_map = new TFile(filename,"");
  if (!file_3D_map || file_3D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  TH3* h_R_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionR_posz");
  TH3* h_P_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionP_posz");
  TH3* h_Z_prz_pos = (TH3*) file_3D_map->Get("hIntDistortionZ_posz");
  //TH3* h_N_prz_pos = (TH3*) file_3D_map->Get("hentries_posz");
  TH3* h_R_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionR_negz");
  TH3* h_P_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionP_negz");
  TH3* h_Z_prz_neg = (TH3*) file_3D_map->Get("hIntDistortionZ_negz");
  //TH3* h_N_prz_neg = (TH3*) file_3D_map->Get("hentries_negz");

  if (!h_R_prz_pos || !h_P_prz_pos || !h_Z_prz_pos ||
      !h_R_prz_neg || !h_P_prz_neg || !h_Z_prz_neg) {
      std::cerr << "Error: missing 3D histograms in file" << std::endl;
      delete file_3D_map;
      return;
  }

  //h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos->GetYaxis()->GetNbins(),h_N_prz_pos->GetYaxis()->GetXmin(),h_N_prz_pos->GetYaxis()->GetXmax());
  h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos->GetYaxis()->GetNbins(),h_R_prz_pos->GetYaxis()->GetXmin(),h_R_prz_pos->GetYaxis()->GetXmax());
  h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos->GetYaxis()->GetNbins(),h_P_prz_pos->GetYaxis()->GetXmin(),h_P_prz_pos->GetYaxis()->GetXmax());
  h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos->GetYaxis()->GetNbins(),h_Z_prz_pos->GetYaxis()->GetXmin(),h_Z_prz_pos->GetYaxis()->GetXmax());
  //h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg->GetYaxis()->GetNbins(),h_N_prz_neg->GetYaxis()->GetXmin(),h_N_prz_neg->GetYaxis()->GetXmax());
  h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg->GetYaxis()->GetNbins(),h_R_prz_neg->GetYaxis()->GetXmin(),h_R_prz_neg->GetYaxis()->GetXmax());
  h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg->GetYaxis()->GetNbins(),h_P_prz_neg->GetYaxis()->GetXmin(),h_P_prz_neg->GetYaxis()->GetXmax());
  h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg->GetYaxis()->GetNbins(),h_Z_prz_neg->GetYaxis()->GetXmin(),h_Z_prz_neg->GetYaxis()->GetXmax());

  // do not save in the file_3D_map
  // directly saved in the stack memory
  //h_N_pos->SetDirectory(0);
  h_R_pos->SetDirectory(0);
  h_P_pos->SetDirectory(0);
  h_Z_pos->SetDirectory(0);
  //h_N_neg->SetDirectory(0);
  h_R_neg->SetDirectory(0);
  h_P_neg->SetDirectory(0);
  h_Z_neg->SetDirectory(0);

  //plot1D_Zbin(h_N_prz_pos,h_N_pos,selectZ, false);
  plot1D_Zbin(h_R_prz_pos,h_R_pos,selectZ, false);
  plot1D_Zbin(h_P_prz_pos,h_P_pos,selectZ, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_pos,h_Z_pos,selectZ, false);
  //plot1D_Zbin(h_N_prz_neg,h_N_neg,-selectZ, false);
  plot1D_Zbin(h_R_prz_neg,h_R_neg,-selectZ, false);
  plot1D_Zbin(h_P_prz_neg,h_P_neg,-selectZ, convert_RP_2_P);
  plot1D_Zbin(h_Z_prz_neg,h_Z_neg,-selectZ, false);

  //h_N_pos->SetLineColor(color); h_N_pos->SetLineWidth(1); h_N_pos->SetFillColor(0); h_N_pos->SetMarkerColor(color);
  h_R_pos->SetLineColor(color); h_R_pos->SetLineWidth(1); h_R_pos->SetFillColor(0); h_R_pos->SetMarkerColor(color);
  h_P_pos->SetLineColor(color); h_P_pos->SetLineWidth(1); h_P_pos->SetFillColor(0); h_P_pos->SetMarkerColor(color);
  h_Z_pos->SetLineColor(color); h_Z_pos->SetLineWidth(1); h_Z_pos->SetFillColor(0); h_Z_pos->SetMarkerColor(color);
  //h_N_neg->SetLineColor(color); h_N_neg->SetLineWidth(1); h_N_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_R_neg->SetLineColor(color); h_R_neg->SetLineWidth(1); h_R_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_P_neg->SetLineColor(color); h_P_neg->SetLineWidth(1); h_P_neg->SetFillColor(0); h_P_neg->SetMarkerColor(color);
  h_Z_neg->SetLineColor(color); h_Z_neg->SetLineWidth(1); h_Z_neg->SetFillColor(0); h_Z_neg->SetMarkerColor(color);

  delete file_3D_map;

  return;
}

void draw1Dmap_2D(TString filename, TString tag,
  TH1*& h_R_pos, TH1*& h_R_neg,
  TH1*& h_P_pos, TH1*& h_P_neg,
  TH1*& h_Z_pos, TH1*& h_Z_neg,
  int color=1)
{
  TFile* file_2D_map = new TFile(filename,"");
  if (!file_2D_map || file_2D_map->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
  }

  TH2* h_R_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionR_posz");
  TH2* h_P_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionP_posz");
  TH2* h_Z_pr_pos = (TH2*) file_2D_map->Get("hIntDistortionZ_posz");
  //TH2* h_N_pr_pos = (TH2*) file_2D_map->Get("hentries_posz");
  TH2* h_R_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionR_negz");
  TH2* h_P_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionP_negz");
  TH2* h_Z_pr_neg = (TH2*) file_2D_map->Get("hIntDistortionZ_negz");
  //TH2* h_N_pr_neg = (TH2*) file_2D_map->Get("hentries_negz");

  if (!h_R_pr_pos || !h_P_pr_pos || !h_Z_pr_pos ||
      !h_R_pr_neg || !h_P_pr_neg || !h_Z_pr_neg) {
      std::cerr << "Error: missing 2D histograms in file" << std::endl;
      delete file_2D_map;
      return;
  }

  //h_N_pos = new TH1F(Form("hentries_posz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_pos->GetYaxis()->GetNbins(),h_N_pr_pos->GetYaxis()->GetXmin(),h_N_pr_pos->GetYaxis()->GetXmax());
  h_R_pos = new TH1F(Form("hIntDistortionR_posz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_pos->GetYaxis()->GetNbins(),h_R_pr_pos->GetYaxis()->GetXmin(),h_R_pr_pos->GetYaxis()->GetXmax());
  h_P_pos = new TH1F(Form("hIntDistortionP_posz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_pos->GetYaxis()->GetNbins(),h_P_pr_pos->GetYaxis()->GetXmin(),h_P_pr_pos->GetYaxis()->GetXmax());
  h_Z_pos = new TH1F(Form("hIntDistortionZ_posz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_pos->GetYaxis()->GetNbins(),h_Z_pr_pos->GetYaxis()->GetXmin(),h_Z_pr_pos->GetYaxis()->GetXmax());
  //h_N_neg = new TH1F(Form("hentries_negz_%s",tag.Data()),Form("N;R (cm);N"),h_N_pr_neg->GetYaxis()->GetNbins(),h_N_pr_neg->GetYaxis()->GetXmin(),h_N_pr_neg->GetYaxis()->GetXmax());
  h_R_neg = new TH1F(Form("hIntDistortionR_negz_%s",tag.Data()),Form("dR;R (cm);dR (cm)"),h_R_pr_neg->GetYaxis()->GetNbins(),h_R_pr_neg->GetYaxis()->GetXmin(),h_R_pr_neg->GetYaxis()->GetXmax());
  h_P_neg = new TH1F(Form("hIntDistortionP_negz_%s",tag.Data()),Form("dphi;R (cm);dphi (rad)"),h_P_pr_neg->GetYaxis()->GetNbins(),h_P_pr_neg->GetYaxis()->GetXmin(),h_P_pr_neg->GetYaxis()->GetXmax());
  h_Z_neg = new TH1F(Form("hIntDistortionZ_negz_%s",tag.Data()),Form("dz;R (cm);dz (cm)"),h_Z_pr_neg->GetYaxis()->GetNbins(),h_Z_pr_neg->GetYaxis()->GetXmin(),h_Z_pr_neg->GetYaxis()->GetXmax());

  // do not save in the file_2D_map
  // directly saved in the stack memory
  //h_N_pos->SetDirectory(0);
  h_R_pos->SetDirectory(0);
  h_P_pos->SetDirectory(0);
  h_Z_pos->SetDirectory(0);
  //h_N_neg->SetDirectory(0);
  h_R_neg->SetDirectory(0);
  h_P_neg->SetDirectory(0);
  h_Z_neg->SetDirectory(0);

  //plot1D(h_N_pr_pos,h_N_pos);
  plot1D(h_R_pr_pos,h_R_pos);
  plot1D(h_P_pr_pos,h_P_pos);
  plot1D(h_Z_pr_pos,h_Z_pos);
  //plot1D(h_N_pr_neg,h_N_neg);
  plot1D(h_R_pr_neg,h_R_neg);
  plot1D(h_P_pr_neg,h_P_neg);
  plot1D(h_Z_pr_neg,h_Z_neg);

  //h_N_pos->SetLineColor(color); h_N_pos->SetLineWidth(1); h_N_pos->SetFillColor(0); h_N_pos->SetMarkerColor(color);
  h_R_pos->SetLineColor(color); h_R_pos->SetLineWidth(1); h_R_pos->SetFillColor(0); h_R_pos->SetMarkerColor(color);
  h_P_pos->SetLineColor(color); h_P_pos->SetLineWidth(1); h_P_pos->SetFillColor(0); h_P_pos->SetMarkerColor(color);
  h_Z_pos->SetLineColor(color); h_Z_pos->SetLineWidth(1); h_Z_pos->SetFillColor(0); h_Z_pos->SetMarkerColor(color);
  //h_N_neg->SetLineColor(color); h_N_neg->SetLineWidth(1); h_N_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_R_neg->SetLineColor(color); h_R_neg->SetLineWidth(1); h_R_neg->SetFillColor(0); h_R_neg->SetMarkerColor(color);
  h_P_neg->SetLineColor(color); h_P_neg->SetLineWidth(1); h_P_neg->SetFillColor(0); h_P_neg->SetMarkerColor(color);
  h_Z_neg->SetLineColor(color); h_Z_neg->SetLineWidth(1); h_Z_neg->SetFillColor(0); h_Z_neg->SetMarkerColor(color);

  delete file_2D_map;

  return;
}