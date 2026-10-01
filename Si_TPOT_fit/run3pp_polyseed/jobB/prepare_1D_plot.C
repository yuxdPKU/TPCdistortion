#include "utilities.h"

void prepare_1D_plot(int run = 79516)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString phiType[3] = { "central", "west", "east" };

  const int nfiles = 2;
  oneDHist hists[nfiles];
  twoDHist hists_rz[nfiles], hists_pr[nfiles], hists_pz[nfiles];
  TString names[nfiles] = 
  {
    "MI_data_reco_FITacts_USEMMSfalse",
    "MI_data_reco_FITacts_USEMMStrue"
  };

  for (int itype = 0; itype < 3; itype++)
  {
    double selectPhi = phi_type2val(phiType[itype]);
    if (std::isnan(selectPhi)) return;

    std::cout << "Use phiType = " << phiType[itype]
          << ", phi = " << selectPhi
          << std::endl;

    // plot distortion vs R, at fixed Z and phi
    std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atZ(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectZ, selectPhi, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("Rootfiles/hist_1Dr_Z%d_%s.root",(int)selectZ,phiType[itype].Data()),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
    }
    
    // plot distortion vs Z, at fixed R and phi
    std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    for (const auto& selectR : selectRs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atR(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectR, selectPhi, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("Rootfiles/hist_1Dz_R%d_%s.root",(int)selectR,phiType[itype].Data()),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
    }

    // plot distortion at R-Z plane, at fixed phi
    for (int i=0; i<nfiles; i++)
    {
      draw2Dmap_atP(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectPhi, hists_rz[i]);
    }
    TFile* ofile = new TFile(Form("Rootfiles/hist_2Drz_%s.root",phiType[itype].Data()),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_rz[i].h_N_pos->Write(); hists_rz[i].h_N_neg->Write();
      hists_rz[i].h_R_pos->Write(); hists_rz[i].h_R_neg->Write();
      hists_rz[i].h_P_pos->Write(); hists_rz[i].h_P_neg->Write();
      hists_rz[i].h_Z_pos->Write(); hists_rz[i].h_Z_neg->Write();
    }
    ofile->Write();
  }

  // plot distortion vs phi, at fixed R and Z
  std::vector<float> selectRs={70, 60, 50, 40, 30};
  std::vector<float> selectZs={5, 10, 15, 30, 60};
  for (const auto& selectR : selectRs)
  {
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        draw1Dmap_atRZ(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectR, selectZ, hists[i], kBlack);
      }

      TFile* ofile = new TFile(Form("Rootfiles/hist_1Dp_R%d_Z%d.root",(int)selectR,(int)selectZ),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write();
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
    }
  }

  // plot distortion at phi-R plane, at fixed Z
  for (const auto& selectZ : selectZs)
  {
    for (int i=0; i<nfiles; i++)
    {
      draw2Dmap_atZ(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectZ, hists_pr[i]);
    }

    TFile* ofile = new TFile(Form("Rootfiles/hist_2Dpr_Z%d.root",(int)selectZ),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_pr[i].h_N_pos->Write(); hists_pr[i].h_N_neg->Write();
      hists_pr[i].h_R_pos->Write(); hists_pr[i].h_R_neg->Write();
      hists_pr[i].h_P_pos->Write(); hists_pr[i].h_P_neg->Write();
      hists_pr[i].h_Z_pos->Write(); hists_pr[i].h_Z_neg->Write();
    }
    ofile->Write();
  }

  // plot distortion at phi-z plane, at fixed R
  for (const auto& selectR : selectRs)
  {
    for (int i=0; i<nfiles; i++)
    {
      draw2Dmap_atR(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectR, hists_pz[i]);
    }

    TFile* ofile = new TFile(Form("Rootfiles/hist_2Dpz_R%d.root",(int)selectR),"recreate");
    ofile->cd();
    for (int i=0; i<nfiles; i++)
    {
      hists_pz[i].h_N_pos->Write(); hists_pz[i].h_N_neg->Write();
      hists_pz[i].h_R_pos->Write(); hists_pz[i].h_R_neg->Write();
      hists_pz[i].h_P_pos->Write(); hists_pz[i].h_P_neg->Write();
      hists_pz[i].h_Z_pos->Write(); hists_pz[i].h_Z_neg->Write();
    }
    ofile->Write();
  }

}
