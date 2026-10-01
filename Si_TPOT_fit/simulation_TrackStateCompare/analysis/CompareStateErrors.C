#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
  constexpr unsigned int mvtxLastLayer = 2;
  constexpr unsigned int inttFirstLayer = 3;
  constexpr unsigned int inttLastLayer = 6;
  constexpr unsigned int tpcFirstLayer = 7;
  constexpr unsigned int tpcLastLayer = 54;

  struct ErrorHistograms
  {
    TH1D* mvtxRphi = nullptr;
    TH1D* mvtxZ = nullptr;
    TH1D* inttRphi = nullptr;
    TH1D* inttZ = nullptr;
    TH1D* tpcRphi = nullptr;
    TH1D* tpcZ = nullptr;
    Long64_t acceptedTracks = 0;
    Long64_t mvtxStates = 0;
    Long64_t inttStates = 0;
    Long64_t tpcStates = 0;
  };

  TH1D* makeHistogram(const TString& name, const TString& xTitle)
  {
    auto* histogram = new TH1D(name, TString::Format(";%s;Normalized states", xTitle.Data()),
                               1000, 0., 0.1);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    return histogram;
  }

  ErrorHistograms readErrors(
      const char* fileName,
      const char* tag,
      Long64_t maxEntries)
  {
    ErrorHistograms output;
    output.mvtxRphi = makeHistogram(
        TString::Format("h_%s_mvtx_rphi", tag), "#sigma_{local X} [cm]");
    output.mvtxZ = makeHistogram(
        TString::Format("h_%s_mvtx_z", tag), "#sigma_{local Y} [cm]");
    output.inttRphi = makeHistogram(
        TString::Format("h_%s_intt_rphi", tag), "#sigma_{local X} [cm]");
    output.inttZ = makeHistogram(
        TString::Format("h_%s_intt_z", tag), "#sigma_{local Y} [cm]");
    output.tpcRphi = makeHistogram(
        TString::Format("h_%s_tpc_rphi", tag), "#sigma_{local X} [cm]");
    output.tpcZ = makeHistogram(
        TString::Format("h_%s_tpc_z", tag), "#sigma_{local Y} [cm]");

    TFile input(fileName, "READ");
    auto* tree = input.IsZombie() ? nullptr :
        dynamic_cast<TTree*>(input.Get("trackStateComparison"));
    if (!tree)
    {
      std::cerr << "Could not read trackStateComparison from " << fileName << std::endl;
      return output;
    }

    Float_t recoTrackPt = std::numeric_limits<float>::quiet_NaN();
    Float_t quality = std::numeric_limits<float>::quiet_NaN();
    UInt_t nstateMvtx = 0;
    UInt_t nstateIntt = 0;
    UInt_t nstateTpc = 0;
    UInt_t nclusMvtx = 0;
    UInt_t nclusIntt = 0;
    UInt_t nclusTpc = 0;
    std::vector<unsigned int>* layers = nullptr;
    std::vector<float>* localXErrors = nullptr;
    std::vector<float>* localYErrors = nullptr;

    const std::array<const char*, 11> branches = {
        "reco_track_pt", "quality",
        "nstate_mvtx", "nstate_intt", "nstate_tpc",
        "nclus_mvtx", "nclus_intt", "nclus_tpc",
        "reco_state_layer", "reco_state_local_x_error", "reco_state_local_y_error"};
    tree->SetBranchStatus("*", 0);
    for (const auto* branch : branches)
    {
      if (!tree->GetBranch(branch))
      {
        std::cerr << "Missing branch " << branch << " in " << fileName << std::endl;
        return output;
      }
      tree->SetBranchStatus(branch, 1);
    }

    tree->SetBranchAddress("reco_track_pt", &recoTrackPt);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("nstate_mvtx", &nstateMvtx);
    tree->SetBranchAddress("nstate_intt", &nstateIntt);
    tree->SetBranchAddress("nstate_tpc", &nstateTpc);
    tree->SetBranchAddress("nclus_mvtx", &nclusMvtx);
    tree->SetBranchAddress("nclus_intt", &nclusIntt);
    tree->SetBranchAddress("nclus_tpc", &nclusTpc);
    tree->SetBranchAddress("reco_state_layer", &layers);
    tree->SetBranchAddress("reco_state_local_x_error", &localXErrors);
    tree->SetBranchAddress("reco_state_local_y_error", &localYErrors);

    const Long64_t entries = maxEntries < 0 ? tree->GetEntries() :
        std::min(maxEntries, tree->GetEntries());
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
      if (tree->GetEntry(entry) <= 0)
      {
        continue;
      }

      const bool passTrackSelection =
          std::isfinite(recoTrackPt) && recoTrackPt >= 0.2F &&
          std::isfinite(quality) && quality <= 100.F &&
          nclusMvtx >= 3U && nclusIntt >= 2U && nclusTpc >= 35U &&
          nstateMvtx >= 3U && nstateIntt >= 2U && nstateTpc >= 35U;
      if (!passTrackSelection || !layers || !localXErrors || !localYErrors ||
          layers->size() != localXErrors->size() ||
          layers->size() != localYErrors->size())
      {
        continue;
      }
      ++output.acceptedTracks;

      for (std::size_t state = 0; state < layers->size(); ++state)
      {
        const auto layer = layers->at(state);
        const double rphiError = localXErrors->at(state);
        const double zError = localYErrors->at(state);
        if (layer <= mvtxLastLayer)
        {
          if (std::isfinite(rphiError) && rphiError >= 0.) output.mvtxRphi->Fill(rphiError);
          if (std::isfinite(zError) && zError >= 0.) output.mvtxZ->Fill(zError);
          ++output.mvtxStates;
        }
        else if (layer >= inttFirstLayer && layer <= inttLastLayer)
        {
          if (std::isfinite(rphiError) && rphiError >= 0.) output.inttRphi->Fill(rphiError);
          if (std::isfinite(zError) && zError >= 0.) output.inttZ->Fill(zError);
          ++output.inttStates;
        }
        else if (layer >= tpcFirstLayer && layer <= tpcLastLayer)
        {
          if (std::isfinite(rphiError) && rphiError >= 0.) output.tpcRphi->Fill(rphiError);
          if (std::isfinite(zError) && zError >= 0.) output.tpcZ->Fill(zError);
          ++output.tpcStates;
        }
      }
    }

    std::cout << tag << ": accepted tracks " << output.acceptedTracks
              << ", MVTX states " << output.mvtxStates
              << ", INTT states " << output.inttStates
              << ", TPC states " << output.tpcStates << std::endl;
    return output;
  }

  void normalize(TH1D* histogram)
  {
    const double integral = histogram->Integral(0, histogram->GetNbinsX() + 1);
    if (integral > 0.) histogram->Scale(1. / integral);
  }

  void style(TH1D* histogram, Color_t color, Style_t style = 1)
  {
    histogram->SetStats(false);
    histogram->SetLineColor(color);
    histogram->SetMarkerColor(color);
    histogram->SetLineStyle(style);
    histogram->SetLineWidth(3);
    histogram->GetXaxis()->SetTitleOffset(1.15);
    histogram->GetYaxis()->SetTitleOffset(1.35);
  }

  double upperDisplayEdge(TH1D* first, TH1D* second)
  {
    const double probability = 0.8;
    double firstQuantile = 0.;
    double secondQuantile = 0.;
    first->GetQuantiles(1, &firstQuantile, &probability);
    second->GetQuantiles(1, &secondQuantile, &probability);
    return std::min(0.1, std::max(0.005, 1.25 * std::max(firstQuantile, secondQuantile)));
  }

  void drawPair(
      TCanvas* canvas,
      int pad,
      TH1D* first,
      TH1D* second,
      const char* firstLabel,
      const char* secondLabel,
      Color_t firstColor,
      Color_t secondColor)
  {
    canvas->cd(pad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);
    normalize(first);
    normalize(second);
    style(first, firstColor);
    style(second, secondColor);
    first->GetXaxis()->SetRangeUser(0., upperDisplayEdge(first, second));
    second->GetXaxis()->SetRangeUser(0., upperDisplayEdge(first, second));
    first->SetMaximum(1.18 * std::max(first->GetMaximum(), second->GetMaximum()));
    first->Draw("HIST");
    second->Draw("HIST SAME");
    auto* legend = new TLegend(0.48, 0.72, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(first, firstLabel, "l");
    legend->AddEntry(second, secondLabel, "l");
    legend->Draw();
  }

  void drawTriple(
      TCanvas* canvas,
      int pad,
      TH1D* first,
      TH1D* second,
      TH1D* third,
      const char* firstLabel,
      const char* secondLabel,
      const char* thirdLabel)
  {
    canvas->cd(pad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);
    normalize(first);
    normalize(second);
    normalize(third);
    style(first, kBlue + 1);
    style(second, kGreen + 2);
    style(third, kRed + 1);

    const double firstTwoEdge = upperDisplayEdge(first, second);
    const double displayEdge = std::max(firstTwoEdge, upperDisplayEdge(first, third));
    first->GetXaxis()->SetRangeUser(0., displayEdge);
    second->GetXaxis()->SetRangeUser(0., displayEdge);
    third->GetXaxis()->SetRangeUser(0., displayEdge);
    first->GetXaxis()->SetNdivisions(505);
    second->GetXaxis()->SetNdivisions(505);
    third->GetXaxis()->SetNdivisions(505);
    first->SetMaximum(1.18 * std::max({
        first->GetMaximum(), second->GetMaximum(), third->GetMaximum()}));
    first->Draw("HIST");
    second->Draw("HIST SAME");
    third->Draw("HIST SAME");
    auto* legend = new TLegend(0.63, 0.68, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(first, firstLabel, "l");
    legend->AddEntry(second, secondLabel, "l");
    legend->AddEntry(third, thirdLabel, "l");
    legend->Draw();
  }

  void saveCanvas(TCanvas* canvas, const TString& directory, const char* name)
  {
    canvas->SaveAs(TString::Format("%s/%s.png", directory.Data(), name));
    canvas->SaveAs(TString::Format("%s/%s.pdf", directory.Data(), name));
  }
}  // namespace

void CompareStateErrors(
    const char* actsRawFile = "all_FITacts_CLUSERRraw_tpcinfit.root",
    const char* actsSimFile = "all_FITacts_CLUSERRsim_tpcinfit.root",
    const char* genfitFile = "all_FITgenfit_tpcinfit.root",
    const char* outputFileName = "hist_state_error_comparison.root",
    const char* figureDirectory = "figure_state_error_comparison",
    Long64_t maxEntries = -1)
{
  gSystem->mkdir(figureDirectory, true);

  auto actsRaw = readErrors(actsRawFile, "acts_raw", maxEntries);
  auto actsSim = readErrors(actsSimFile, "acts_sim", maxEntries);
  auto genfit = readErrors(genfitFile, "genfit", maxEntries);

  TFile output(outputFileName, "RECREATE");
  const std::array<TH1D*, 18> allHistograms = {
      actsRaw.mvtxRphi, actsRaw.mvtxZ, actsRaw.inttRphi, actsRaw.inttZ,
      actsRaw.tpcRphi, actsRaw.tpcZ,
      actsSim.mvtxRphi, actsSim.mvtxZ, actsSim.inttRphi, actsSim.inttZ,
      actsSim.tpcRphi, actsSim.tpcZ,
      genfit.mvtxRphi, genfit.mvtxZ, genfit.inttRphi, genfit.inttZ,
      genfit.tpcRphi, genfit.tpcZ};
  for (auto* histogram : allHistograms) histogram->Write();

  auto* fitterCanvas = new TCanvas(
      "c_tpc_acts_raw_vs_genfit", "TPC local state errors: ACTS and GENFIT", 1400, 600);
  fitterCanvas->Divide(2, 1);
  drawPair(fitterCanvas, 1, actsRaw.tpcRphi, genfit.tpcRphi,
           "ACTS, cluster error raw", "GENFIT", kRed + 1, kBlue + 1);
  drawPair(fitterCanvas, 2, actsRaw.tpcZ, genfit.tpcZ,
           "ACTS, cluster error raw", "GENFIT", kRed + 1, kBlue + 1);
  fitterCanvas->Write();
  saveCanvas(fitterCanvas, figureDirectory, "tpc_state_error_acts_raw_vs_genfit");

  auto* detectorCanvas = new TCanvas(
      "c_acts_raw_mvtx_intt_tpc", "ACTS local state errors: MVTX, INTT, and TPC", 1400, 600);
  detectorCanvas->Divide(2, 1);
  drawTriple(detectorCanvas, 1, actsRaw.mvtxRphi, actsRaw.inttRphi, actsRaw.tpcRphi,
             "MVTX states", "INTT states", "TPC states");
  drawTriple(detectorCanvas, 2, actsRaw.mvtxZ, actsRaw.inttZ, actsRaw.tpcZ,
             "MVTX states", "INTT states", "TPC states");
  detectorCanvas->Write();
  saveCanvas(detectorCanvas, figureDirectory, "acts_raw_state_error_mvtx_intt_tpc");

  auto* clusterErrorCanvas = new TCanvas(
      "c_tpc_acts_raw_vs_sim", "TPC local state errors: ACTS cluster errors", 1400, 600);
  clusterErrorCanvas->Divide(2, 1);
  drawPair(clusterErrorCanvas, 1, actsRaw.tpcRphi, actsSim.tpcRphi,
           "Cluster error raw", "Cluster error simulation", kBlack, kMagenta + 1);
  drawPair(clusterErrorCanvas, 2, actsRaw.tpcZ, actsSim.tpcZ,
           "Cluster error raw", "Cluster error simulation", kBlack, kMagenta + 1);
  clusterErrorCanvas->Write();
  saveCanvas(clusterErrorCanvas, figureDirectory, "acts_tpc_state_error_cluster_raw_vs_sim");

  output.Close();
  std::cout << "Wrote " << outputFileName << " and figures under "
            << figureDirectory << std::endl;
}
