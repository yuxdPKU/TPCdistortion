#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

void plot1D_Zbin(TH3* h3, TH1* h1, float z)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
	  //cout<<"pbin = "<<xbin<<", rbin = "<<i<<" , zbin = "<<zbin<<" , content = "<<h3->GetBinContent(xbin, i, zbin)<<endl;
      h1->SetBinContent(i, h3->GetBinContent(xbin, i, zbin));
      h1->SetBinError(i, h3->GetBinError(xbin, i, zbin));
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

const int nrun = 1;
//int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
//int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
//int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
//int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
int mbdrates[nrun] = {400};
int runs[nrun] = {53877};

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun], *hcdb_RP_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun], *hcdb_RP_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun], *hcdb_RP_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun], *hcdb_RP_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  auto rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", runs[i]);
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  cdbfilename[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
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
  hcdb_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cdb",runs[i]),Form("dphi;R (cm);dphi (rad)"),hcdb_P_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_N_neg[i] = new TH1F(Form("hentries_negz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cdb",runs[i]),Form("dphi;R (cm);dphi (rad)"),hcdb_P_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmax());
  plot1D(hcdb_N_pr_pos[i],hcdb_N_pos[i]);
  plot1D(hcdb_R_pr_pos[i],hcdb_R_pos[i]);
  plot1D(hcdb_P_pr_pos[i],hcdb_P_pos[i]);
  plot1D(hcdb_Z_pr_pos[i],hcdb_Z_pos[i]);
  plot1D(hcdb_N_pr_neg[i],hcdb_N_neg[i]);
  plot1D(hcdb_R_pr_neg[i],hcdb_R_neg[i]);
  plot1D(hcdb_P_pr_neg[i],hcdb_P_neg[i]);
  plot1D(hcdb_Z_pr_neg[i],hcdb_Z_neg[i]);
  hcdb_N_pos[i]->SetLineColor(1); hcdb_N_pos[i]->SetLineWidth(1); hcdb_N_pos[i]->SetFillColor(0); hcdb_N_pos[i]->SetMarkerColor(i+2);
  hcdb_R_pos[i]->SetLineColor(1); hcdb_R_pos[i]->SetLineWidth(1); hcdb_R_pos[i]->SetFillColor(0); hcdb_R_pos[i]->SetMarkerColor(i+2);
  hcdb_P_pos[i]->SetLineColor(1); hcdb_P_pos[i]->SetLineWidth(1); hcdb_P_pos[i]->SetFillColor(0); hcdb_P_pos[i]->SetMarkerColor(i+2);
  hcdb_Z_pos[i]->SetLineColor(1); hcdb_Z_pos[i]->SetLineWidth(1); hcdb_Z_pos[i]->SetFillColor(0); hcdb_Z_pos[i]->SetMarkerColor(i+2);
  hcdb_N_neg[i]->SetLineColor(1); hcdb_N_neg[i]->SetLineWidth(1); hcdb_N_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(i+2);
  hcdb_R_neg[i]->SetLineColor(1); hcdb_R_neg[i]->SetLineWidth(1); hcdb_R_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(i+2);
  hcdb_P_neg[i]->SetLineColor(1); hcdb_P_neg[i]->SetLineWidth(1); hcdb_P_neg[i]->SetFillColor(0); hcdb_P_neg[i]->SetMarkerColor(i+2);
  hcdb_Z_neg[i]->SetLineColor(1); hcdb_Z_neg[i]->SetLineWidth(1); hcdb_Z_neg[i]->SetFillColor(0); hcdb_Z_neg[i]->SetMarkerColor(i+2);
}

//std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
std::vector<float> selectZs={5};
for (const auto& selectZ : selectZs)
{

TFile* file_3D_map[nrun];
TH3 *h_N_prz_pos[nrun], *h_R_prz_pos[nrun], *h_P_prz_pos[nrun], *h_Z_prz_pos[nrun], *h_RP_prz_pos[nrun];
TH3 *h_N_prz_neg[nrun], *h_R_prz_neg[nrun], *h_P_prz_neg[nrun], *h_Z_prz_neg[nrun], *h_RP_prz_neg[nrun];
TH1 *h_N_pos[nrun], *h_R_pos[nrun], *h_P_pos[nrun], *h_Z_pos[nrun], *h_RP_pos[nrun];
TH1 *h_N_neg[nrun], *h_R_neg[nrun], *h_P_neg[nrun], *h_Z_neg[nrun], *h_RP_neg[nrun];
double ymax_N=-100, ymin_N=100;
double ymax_R=-100, ymin_R=100;
double ymax_P=-100, ymin_P=100;
double ymax_Z=-100, ymin_Z=100;
double ymax_RP=-100, ymin_RP=100;
for (int i = 0; i < nrun; i++)
{
  file_3D_map[i] = new TFile(Form("../Rootfiles/Distortions_2D_mm_%d_rz.root",runs[i]),"");
  h_R_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_posz");
  h_RP_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_posz");
  h_N_prz_pos[i] = (TH3*) file_3D_map[i]->Get("hentries_posz");
  h_R_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionZ_negz");
  h_RP_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hIntDistortionRP_negz");
  h_N_prz_neg[i] = (TH3*) file_3D_map[i]->Get("hentries_negz");

  h_N_pos[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos[i]->GetYaxis()->GetNbins(),h_N_prz_pos[i]->GetYaxis()->GetXmin(),h_N_prz_pos[i]->GetYaxis()->GetXmax());
  h_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos[i]->GetYaxis()->GetNbins(),h_R_prz_pos[i]->GetYaxis()->GetXmin(),h_R_prz_pos[i]->GetYaxis()->GetXmax());
  h_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("dphi @ Z=%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_pos[i]->GetYaxis()->GetNbins(),h_P_prz_pos[i]->GetYaxis()->GetXmin(),h_P_prz_pos[i]->GetYaxis()->GetXmax());
  h_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos[i]->GetYaxis()->GetNbins(),h_Z_prz_pos[i]->GetYaxis()->GetXmin(),h_Z_prz_pos[i]->GetYaxis()->GetXmax());
  h_RP_pos[i] = new TH1F(Form("hIntDistortionRP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_RP_prz_pos[i]->GetYaxis()->GetNbins(),h_RP_prz_pos[i]->GetYaxis()->GetXmin(),h_RP_prz_pos[i]->GetYaxis()->GetXmax());
  h_N_neg[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg[i]->GetYaxis()->GetNbins(),h_N_prz_neg[i]->GetYaxis()->GetXmin(),h_N_prz_neg[i]->GetYaxis()->GetXmax());
  h_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg[i]->GetYaxis()->GetNbins(),h_R_prz_neg[i]->GetYaxis()->GetXmin(),h_R_prz_neg[i]->GetYaxis()->GetXmax());
  h_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("dphi @ Z=-%d cm;R (cm);dphi (rad)",(int)selectZ),h_P_prz_neg[i]->GetYaxis()->GetNbins(),h_P_prz_neg[i]->GetYaxis()->GetXmin(),h_P_prz_neg[i]->GetYaxis()->GetXmax());
  h_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg[i]->GetYaxis()->GetNbins(),h_Z_prz_neg[i]->GetYaxis()->GetXmin(),h_Z_prz_neg[i]->GetYaxis()->GetXmax());
  h_RP_neg[i] = new TH1F(Form("hIntDistortionRP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_RP_prz_neg[i]->GetYaxis()->GetNbins(),h_RP_prz_neg[i]->GetYaxis()->GetXmin(),h_RP_prz_neg[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos[i],h_N_pos[i],selectZ);
  plot1D_Zbin(h_R_prz_pos[i],h_R_pos[i],selectZ);
  plot1D_Zbin(h_P_prz_pos[i],h_P_pos[i],selectZ);
  plot1D_Zbin(h_Z_prz_pos[i],h_Z_pos[i],selectZ);
  plot1D_Zbin(h_RP_prz_pos[i],h_RP_pos[i],selectZ);
  plot1D_Zbin(h_N_prz_neg[i],h_N_neg[i],-selectZ);
  plot1D_Zbin(h_R_prz_neg[i],h_R_neg[i],-selectZ);
  plot1D_Zbin(h_P_prz_neg[i],h_P_neg[i],-selectZ);
  plot1D_Zbin(h_Z_prz_neg[i],h_Z_neg[i],-selectZ);
  plot1D_Zbin(h_RP_prz_neg[i],h_RP_neg[i],-selectZ);
  h_N_pos[i]->SetLineColor(i+2); h_N_pos[i]->SetLineWidth(1); h_N_pos[i]->SetFillColor(0); h_N_pos[i]->SetMarkerColor(i+2);
  h_R_pos[i]->SetLineColor(i+2); h_R_pos[i]->SetLineWidth(1); h_R_pos[i]->SetFillColor(0); h_R_pos[i]->SetMarkerColor(i+2);
  h_P_pos[i]->SetLineColor(i+2); h_P_pos[i]->SetLineWidth(1); h_P_pos[i]->SetFillColor(0); h_P_pos[i]->SetMarkerColor(i+2);
  h_Z_pos[i]->SetLineColor(i+2); h_Z_pos[i]->SetLineWidth(1); h_Z_pos[i]->SetFillColor(0); h_Z_pos[i]->SetMarkerColor(i+2);
  h_RP_pos[i]->SetLineColor(i+2); h_RP_pos[i]->SetLineWidth(1); h_RP_pos[i]->SetFillColor(0); h_RP_pos[i]->SetMarkerColor(i+2);
  h_N_neg[i]->SetLineColor(i+2); h_N_neg[i]->SetLineWidth(1); h_N_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_R_neg[i]->SetLineColor(i+2); h_R_neg[i]->SetLineWidth(1); h_R_neg[i]->SetFillColor(0); h_R_neg[i]->SetMarkerColor(i+2);
  h_P_neg[i]->SetLineColor(i+2); h_P_neg[i]->SetLineWidth(1); h_P_neg[i]->SetFillColor(0); h_P_neg[i]->SetMarkerColor(i+2);
  h_Z_neg[i]->SetLineColor(i+2); h_Z_neg[i]->SetLineWidth(1); h_Z_neg[i]->SetFillColor(0); h_Z_neg[i]->SetMarkerColor(i+2);
  h_RP_neg[i]->SetLineColor(i+2); h_RP_neg[i]->SetLineWidth(1); h_RP_neg[i]->SetFillColor(0); h_RP_neg[i]->SetMarkerColor(i+2);

  if (h_N_neg[i]->GetMaximum()>ymax_N) ymax_N = h_N_neg[i]->GetMaximum();
  if (h_N_pos[i]->GetMaximum()>ymax_N) ymax_N = h_N_pos[i]->GetMaximum();
  if (h_N_neg[i]->GetMinimum()<ymin_N) ymin_N = h_N_neg[i]->GetMinimum();
  if (h_N_pos[i]->GetMinimum()<ymin_N) ymin_N = h_N_pos[i]->GetMinimum();
  if (hcdb_N_neg[i]->GetMaximum()>ymax_N) ymax_N = hcdb_N_neg[i]->GetMaximum();
  if (hcdb_N_pos[i]->GetMaximum()>ymax_N) ymax_N = hcdb_N_pos[i]->GetMaximum();
  if (hcdb_N_neg[i]->GetMinimum()<ymin_N) ymin_N = hcdb_N_neg[i]->GetMinimum();
  if (hcdb_N_pos[i]->GetMinimum()<ymin_N) ymin_N = hcdb_N_pos[i]->GetMinimum();

  if (h_R_neg[i]->GetMaximum()>ymax_R) ymax_R = h_R_neg[i]->GetMaximum();
  if (h_R_pos[i]->GetMaximum()>ymax_R) ymax_R = h_R_pos[i]->GetMaximum();
  if (h_R_neg[i]->GetMinimum()<ymin_R) ymin_R = h_R_neg[i]->GetMinimum();
  if (h_R_pos[i]->GetMinimum()<ymin_R) ymin_R = h_R_pos[i]->GetMinimum();
  if (hcdb_R_neg[i]->GetMaximum()>ymax_R) ymax_R = hcdb_R_neg[i]->GetMaximum();
  if (hcdb_R_pos[i]->GetMaximum()>ymax_R) ymax_R = hcdb_R_pos[i]->GetMaximum();
  if (hcdb_R_neg[i]->GetMinimum()<ymin_R) ymin_R = hcdb_R_neg[i]->GetMinimum();
  if (hcdb_R_pos[i]->GetMinimum()<ymin_R) ymin_R = hcdb_R_pos[i]->GetMinimum();

  if (h_P_neg[i]->GetMaximum()>ymax_P) ymax_P = h_P_neg[i]->GetMaximum();
  if (h_P_pos[i]->GetMaximum()>ymax_P) ymax_P = h_P_pos[i]->GetMaximum();
  if (h_P_neg[i]->GetMinimum()<ymin_P) ymin_P = h_P_neg[i]->GetMinimum();
  if (h_P_pos[i]->GetMinimum()<ymin_P) ymin_P = h_P_pos[i]->GetMinimum();

  if (h_Z_neg[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_neg[i]->GetMaximum();
  if (h_Z_pos[i]->GetMaximum()>ymax_Z) ymax_Z = h_Z_pos[i]->GetMaximum();
  if (h_Z_neg[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_neg[i]->GetMinimum();
  if (h_Z_pos[i]->GetMinimum()<ymin_Z) ymin_Z = h_Z_pos[i]->GetMinimum();
  if (hcdb_Z_neg[i]->GetMaximum()>ymax_Z) ymax_Z = hcdb_Z_neg[i]->GetMaximum();
  if (hcdb_Z_pos[i]->GetMaximum()>ymax_Z) ymax_Z = hcdb_Z_pos[i]->GetMaximum();
  if (hcdb_Z_neg[i]->GetMinimum()<ymin_Z) ymin_Z = hcdb_Z_neg[i]->GetMinimum();
  if (hcdb_Z_pos[i]->GetMinimum()<ymin_Z) ymin_Z = hcdb_Z_pos[i]->GetMinimum();

  if (h_RP_neg[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_neg[i]->GetMaximum();
  if (h_RP_pos[i]->GetMaximum()>ymax_RP) ymax_RP = h_RP_pos[i]->GetMaximum();
  if (h_RP_neg[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_neg[i]->GetMinimum();
  if (h_RP_pos[i]->GetMinimum()<ymin_RP) ymin_RP = h_RP_pos[i]->GetMinimum();
  if (hcdb_P_neg[i]->GetMaximum()>ymax_RP) ymax_RP = hcdb_P_neg[i]->GetMaximum();
  if (hcdb_P_pos[i]->GetMaximum()>ymax_RP) ymax_RP = hcdb_P_pos[i]->GetMaximum();
  if (hcdb_P_neg[i]->GetMinimum()<ymin_RP) ymin_RP = hcdb_P_neg[i]->GetMinimum();
  if (hcdb_P_pos[i]->GetMinimum()<ymin_RP) ymin_RP = hcdb_P_pos[i]->GetMinimum();

}

if (ymax_N>0) ymax_N *= 1.1; else ymax_N = 1;
if (ymax_R>0) ymax_R *= 1.1; else ymax_R = 0;
if (ymax_P>0) ymax_P *= 1.1; else ymax_P = 0;
if (ymax_Z>0) ymax_Z *= 1.1; else ymax_Z = 0;
if (ymax_RP>0) ymax_RP *= 1.1; else ymax_RP = 0;

if (ymin_N<0) ymin_N *= 1.1; else ymin_N = 0.9;
if (ymin_R<0) ymin_R *= 1.1; else ymin_R = 0.0;
if (ymin_P<0) ymin_P *= 1.1; else ymin_P = 0.0;
if (ymin_Z<0) ymin_Z *= 1.1; else ymin_Z = 0.0;
if (ymin_RP<0) ymin_RP *= 1.1; else ymin_RP = 0.0;

for (int i = 0; i < nrun; i++)
{
  h_N_pos[i]->SetMinimum(ymin_N); h_N_pos[i]->SetMaximum(ymax_N);
  h_R_pos[i]->SetMinimum(ymin_R); h_R_pos[i]->SetMaximum(ymax_R);
  h_P_pos[i]->SetMinimum(ymin_P); h_P_pos[i]->SetMaximum(ymax_P);
  h_Z_pos[i]->SetMinimum(ymin_Z); h_Z_pos[i]->SetMaximum(ymax_Z);
  h_RP_pos[i]->SetMinimum(ymin_RP); h_RP_pos[i]->SetMaximum(ymax_RP);

  h_N_neg[i]->SetMinimum(ymin_N); h_N_neg[i]->SetMaximum(ymax_N);
  h_R_neg[i]->SetMinimum(ymin_R); h_R_neg[i]->SetMaximum(ymax_R);
  h_P_neg[i]->SetMinimum(ymin_P); h_P_neg[i]->SetMaximum(ymax_P);
  h_Z_neg[i]->SetMinimum(ymin_Z); h_Z_neg[i]->SetMaximum(ymax_Z);
  h_RP_neg[i]->SetMinimum(ymin_RP); h_RP_neg[i]->SetMaximum(ymax_RP);
}

TFile* infile_pos = new TFile("../QA/hist_rdphi_vs_r_pos.root","");
TGraphErrors *graph_full_pos = (TGraphErrors*) infile_pos->Get("graph_full");
TGraphErrors *graph_sub_pos = (TGraphErrors*) infile_pos->Get("graph_sub");
graph_full_pos->SetLineColor(kViolet);
graph_sub_pos->SetLineColor(kBlue);
TFile* infile_neg = new TFile("../QA/hist_rdphi_vs_r_neg.root","");
TGraphErrors *graph_full_neg = (TGraphErrors*) infile_neg->Get("graph_full");
TGraphErrors *graph_sub_neg = (TGraphErrors*) infile_neg->Get("graph_sub");
graph_full_neg->SetLineColor(kViolet);
graph_sub_neg->SetLineColor(kBlue);

TLegend *legend_pos = new TLegend(0.45, 0.6, 0.9, 0.9);
legend_pos->SetHeader("Run 53877");
for (int i=0; i<nrun; i++) legend_pos->AddEntry(h_P_pos[i], Form("Si-TPOT fit map Z=5cm"), "l");
for (int i=0; i<nrun; i++) legend_pos->AddEntry(hcdb_P_pos[i], Form("CDB lamination fit map north side"), "l");
legend_pos->AddEntry(graph_full_pos,"Raw residual 0<Z<10 cm, fit full region");
legend_pos->AddEntry(graph_sub_pos,"Raw residual 0<Z<10 cm, fit main peak");
legend_pos->SetTextSize(0.03);

TLegend *legend_neg = new TLegend(0.45, 0.6, 0.9, 0.9);
legend_neg->SetHeader("Run 53877");
for (int i=0; i<nrun; i++) legend_neg->AddEntry(h_P_neg[i], Form("Si-TPOT fit map Z=-5cm"), "l");
for (int i=0; i<nrun; i++) legend_neg->AddEntry(hcdb_P_neg[i], Form("CDB lamination fit map south side"), "l");
legend_neg->AddEntry(graph_full_neg,"Raw residual -10<Z<0 cm, fit full region");
legend_neg->AddEntry(graph_sub_neg,"Raw residual -10<Z<0 cm, fit main peak");
legend_neg->SetTextSize(0.03);

TCanvas* can = new TCanvas("can","",800,1200);
can->Divide(1,2);
can->cd(1);
gPad->SetRightMargin(0.01);
for (int i=0; i<nrun; i++) h_RP_pos[i]->SetMinimum(-0.65);
for (int i=0; i<nrun; i++) h_RP_pos[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) hcdb_P_pos[i]->Draw("e,same");
graph_full_pos->Draw("P,same");
graph_sub_pos->Draw("P,same");
legend_pos->Draw();
can->cd(2);
gPad->SetRightMargin(0.01);
for (int i=0; i<nrun; i++) h_RP_neg[i]->SetMinimum(-0.65);
for (int i=0; i<nrun; i++) h_RP_neg[i]->Draw("hist,e,same");
for (int i=0; i<nrun; i++) hcdb_P_neg[i]->Draw("e,same");
graph_full_neg->Draw("P,same");
graph_sub_neg->Draw("P,same");
legend_neg->Draw();

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d.pdf",(int)selectZ));

delete can;
}

}
