#include "utilities.h"

void prepare_1D_plot(int run = 29)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString phiType[3] = { "central", "west", "east" };

  const int nfiles = 13;
  oneDHist hists[nfiles];
  TString names[nfiles] = 
  {
    "Input",
    "MI_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue",
    "MI_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue",
    "MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue",
    "MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue",
    "MI_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue",
    "MI_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse",
    "MI_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue"
  };

  for (int itype = 0; itype < 3; itype++)
  {
    double selectPhi = phi_type2val(phiType[itype]);
    if (std::isnan(selectPhi)) return;

    std::cout << "Use phiType = " << phiType[itype]
          << ", phi = " << selectPhi
          << std::endl;

    std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
    for (const auto& selectZ : selectZs)
    {
      for (int i=0; i<nfiles; i++)
      {
        if (i>0)
        {
          draw1Dmap_atZ(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectZ, selectPhi, hists[i], kBlack);
        }
        else if (i==0)
        {
          draw1Dmap_atZ(Form("/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_converted.root"), names[i], selectZ, selectPhi, hists[i], kBlack, true, 1, false);
        }
      }

      TFile* ofile = new TFile(Form("Rootfiles/hist_1Dr_Z%d_%s.root",(int)selectZ,phiType[itype].Data()),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        if (i>0) { hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write(); }
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();

    }
    
    std::vector<float> selectRs={75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    for (const auto& selectR : selectRs)
    {
      for (int i=0; i<nfiles; i++)
      {
        if (i>0)
        {
          draw1Dmap_atR(Form("./Rootfiles/Distortions_full_mm_%d_%s.root", run, names[i].Data()), names[i], selectR, selectPhi, hists[i], kBlack);
        }
        else if (i==0)
        {
          draw1Dmap_atR(Form("/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps/average_minus_static_distortion_converted.root"), names[i], selectR, selectPhi, hists[i], kBlack, true, 1, false);
        }
      }

      TFile* ofile = new TFile(Form("Rootfiles/hist_1Dz_R%d_%s.root",(int)selectR,phiType[itype].Data()),"recreate");
      ofile->cd();
      for (int i=0; i<nfiles; i++)
      {
        if (i>0) { hists[i].h_N_pos->Write(); hists[i].h_N_neg->Write(); }
        hists[i].h_R_pos->Write(); hists[i].h_R_neg->Write();
        hists[i].h_P_pos->Write(); hists[i].h_P_neg->Write();
        hists[i].h_Z_pos->Write(); hists[i].h_Z_neg->Write();
      }
      ofile->Write();
    
    }

  }

}
