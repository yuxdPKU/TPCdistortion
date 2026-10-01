#include "utilities.h"

void draw1D_phi_fitter()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  std::vector<float> selectRs={70, 60, 50, 40, 30};
  for (const auto& selectR : selectRs)
  {

  std::vector<float> selectZs={5, 10, 15, 30, 60};
  for (const auto& selectZ : selectZs)
  {
    oneDHist hists_mmsfalse{};
    get1Dhist(Form("Rootfiles/hist_1Dp_R%d_Z%d.root", (int) selectR, (int)selectZ), "MI_data_reco_FITacts_USEMMSfalse", hists_mmsfalse, kBlack);

    if (!hists_mmsfalse.h_R_pos) {
      std::cerr << "Error: some histograms are null for R = " << selectR << ", Z = " << selectZ << std::endl;
      continue;
    }

    std::vector<oneDHist> allHists = {
      hists_mmsfalse
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
    can->cd(2);
    hists_mmsfalse.h_P_pos->Draw("hist,e,same");
    can->cd(3);
    gPad->SetLogy(0);
    hists_mmsfalse.h_R_pos->Draw("hist,e,same");
    can->cd(4);
    gPad->SetLogy(0);
    hists_mmsfalse.h_Z_pos->Draw("hist,e,same");
    can->cd(5);
    gPad->SetLogy(1);
    hists_mmsfalse.h_N_neg->Draw("hist,e,same");
    can->cd(6);
    hists_mmsfalse.h_P_neg->Draw("hist,e,same");
    can->cd(7);
    gPad->SetLogy(0);
    hists_mmsfalse.h_R_neg->Draw("hist,e,same");
    can->cd(8);
    gPad->SetLogy(0);
    hists_mmsfalse.h_Z_neg->Draw("hist,e,same");

    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("figure/resid_vsP_from3D_atR%d_Z%d.pdf",(int)selectR,(int)selectZ));

    TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetHeader("ana573_2026p003_v001, Run 79516, poly seed");
    legend->AddEntry(hists_mmsfalse.h_P_pos, Form("No TPOT"), "F");
    legend->Draw();
    TCanvas* can_leg = new TCanvas("can_leg","",2000,1200);
    legend->Draw();
    can_leg->SaveAs(Form("figure/resid_vsP_from3D_leg.pdf"));

    delete can;
    delete can_leg;

  }
  }

}
