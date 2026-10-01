#include "utilities.h"

void draw1D_r_fitter()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  std::vector<TString> phiTypes = {"central", "west", "east"};
  for (const auto& phiType : phiTypes)
  {
  oneDHist hists_mmsfalse{};
  oneDHist hists_mmstrue{};

  std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
  for (const auto& selectZ : selectZs)
  {
    get1Dhist(Form("Rootfiles/hist_1Dr_Z%d_%s.root", (int)selectZ, phiType.Data()), "MI_data_reco_FITacts_USEMMSfalse", hists_mmsfalse, kBlack);
    get1Dhist(Form("Rootfiles/hist_1Dr_Z%d_%s.root", (int)selectZ, phiType.Data()), "MI_data_reco_FITacts_USEMMStrue", hists_mmstrue, kRed - 7);

    if (!hists_mmsfalse.h_R_pos || !hists_mmstrue.h_R_pos) {
      std::cerr << "Error: some histograms are null for Z = " << selectZ << std::endl;
      continue;
    }

    std::vector<oneDHist> allHists = {
      hists_mmsfalse,
      hists_mmstrue
    };

    std::vector<TH1*> hists_P; hists_P.clear();
    std::vector<TH1*> hists_R; hists_R.clear();
    std::vector<TH1*> hists_Z; hists_Z.clear();

    AddOneDHistVectorToVectors(allHists, hists_P, hists_R, hists_Z);

    std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
    std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
    std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);

    TCanvas* can = new TCanvas("can","",2400,1200);
    can->Divide(4,2);
    can->cd(1);
    gPad->SetLogy(1);
    hists_mmsfalse.h_N_pos->Draw("hist,e,same");
    hists_mmstrue.h_N_pos->Draw("hist,e,same");
    can->cd(2);
    //hists_acts.h_P_pos->SetMaximum(15e-3);
    hists_mmsfalse.h_P_pos->Draw("hist,e,same");
    hists_mmstrue.h_P_pos->Draw("hist,e,same");
    can->cd(3);
    gPad->SetLogy(0);
    hists_mmsfalse.h_R_pos->Draw("hist,e,same");
    hists_mmstrue.h_R_pos->Draw("hist,e,same");
    can->cd(4);
    gPad->SetLogy(0);
    hists_mmsfalse.h_Z_pos->Draw("hist,e,same");
    hists_mmstrue.h_Z_pos->Draw("hist,e,same");
    can->cd(5);
    gPad->SetLogy(1);
    hists_mmsfalse.h_N_neg->Draw("hist,e,same");
    hists_mmstrue.h_N_neg->Draw("hist,e,same");
    can->cd(6);
    //hists_mmsfalse.h_P_neg->SetMaximum(15e-3);
    hists_mmsfalse.h_P_neg->Draw("hist,e,same");
    hists_mmstrue.h_P_neg->Draw("hist,e,same");
    can->cd(7);
    gPad->SetLogy(0);
    hists_mmsfalse.h_R_neg->Draw("hist,e,same");
    hists_mmstrue.h_R_neg->Draw("hist,e,same");
    can->cd(8);
    gPad->SetLogy(0);
    hists_mmsfalse.h_Z_neg->Draw("hist,e,same");
    hists_mmstrue.h_Z_neg->Draw("hist,e,same");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("figure/resid_vsR_from3D_atZ%d_%s.pdf",(int)selectZ,phiType.Data()));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("ana573_2026p003_v001, Run 79516, poly seed");
    legend->AddEntry(hists_mmsfalse.h_P_pos, Form("No TPOT"), "F");
    legend->AddEntry(hists_mmstrue.h_P_pos, Form("With TPOT"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_%s.pdf",phiType.Data()));

    delete can;
    delete can_leg;

  }
  }

}
