#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace
{
  TH1* cloneOneDimensionalHistogram(
      TFile& inputFile, const std::string& histogramName, const char* suffix)
  {
    auto* source = dynamic_cast<TH1*>(inputFile.Get(histogramName.c_str()));
    if (!source || source->GetDimension() != 1 || dynamic_cast<TProfile*>(source))
    {
      return nullptr;
    }

    auto* histogram = dynamic_cast<TH1*>(
        source->Clone(TString::Format("%s_%s", histogramName.c_str(), suffix)));
    if (histogram)
    {
      histogram->SetDirectory(nullptr);
    }
    return histogram;
  }

  void normalizeHistogram(TH1* histogram, bool normalize)
  {
    if (!histogram || !normalize)
    {
      return;
    }

    const double integral = histogram->Integral(0, histogram->GetNbinsX() + 1);
    if (integral > 0.)
    {
      histogram->Scale(1. / integral);
      histogram->GetYaxis()->SetTitle("Normalized entries");
    }
  }

  void styleHistogram(TH1* histogram, Color_t color, Style_t lineStyle = 1)
  {
    histogram->SetStats(false);
    histogram->SetLineColor(color);
    histogram->SetMarkerColor(color);
    histogram->SetLineStyle(lineStyle);
    histogram->SetLineWidth(3);
  }

  void prepareCanvas(TCanvas* canvas)
  {
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.15);
    canvas->SetTopMargin(0.10);
  }

  void saveCanvas(TCanvas* canvas, const TString& figureDirectory, const TString& tag)
  {
    canvas->SaveAs(TString::Format("%s/%s.png", figureDirectory.Data(), tag.Data()));
    canvas->SaveAs(TString::Format("%s/%s.pdf", figureDirectory.Data(), tag.Data()));
  }
}  // namespace

void CompareRecoAlgorithms(
    const char* actsFileName = "hist_acts_includesecondaries_fullphi_minpt0p2.root",
    const char* genfitFileName = "hist_genfit_includesecondaries_fullphi_minpt0p2.root",
    const char* outputFileName = "hist_reco_comparison.root",
    const char* figureDirectory = "figure_compare",
    bool normalize = true)
{
  TFile actsFile(actsFileName, "READ");
  TFile genfitFile(genfitFileName, "READ");
  if (actsFile.IsZombie() || genfitFile.IsZombie())
  {
    std::cerr << "Could not open both input files: " << actsFileName
              << " and " << genfitFileName << std::endl;
    return;
  }

  gSystem->mkdir(figureDirectory, true);
  TFile outputFile(outputFileName, "RECREATE");

  // Four-way track-pT comparison. Color identifies the reconstruction,
  // while line style identifies truth (dashed) or reconstructed (solid) pT.
  std::unique_ptr<TH1> actsTruthPt(
      cloneOneDimensionalHistogram(actsFile, "h_truth_track_pt", "acts_compare"));
  std::unique_ptr<TH1> actsRecoPt(
      cloneOneDimensionalHistogram(actsFile, "h_reco_track_pt", "acts_compare"));
  std::unique_ptr<TH1> genfitTruthPt(
      cloneOneDimensionalHistogram(genfitFile, "h_truth_track_pt", "genfit_compare"));
  std::unique_ptr<TH1> genfitRecoPt(
      cloneOneDimensionalHistogram(genfitFile, "h_reco_track_pt", "genfit_compare"));

  if (actsTruthPt && actsRecoPt && genfitTruthPt && genfitRecoPt)
  {
    normalizeHistogram(actsTruthPt.get(), normalize);
    normalizeHistogram(actsRecoPt.get(), normalize);
    normalizeHistogram(genfitTruthPt.get(), normalize);
    normalizeHistogram(genfitRecoPt.get(), normalize);
    styleHistogram(actsTruthPt.get(), kRed + 1, 2);
    styleHistogram(actsRecoPt.get(), kRed + 1, 1);
    styleHistogram(genfitTruthPt.get(), kBlue + 1, 2);
    styleHistogram(genfitRecoPt.get(), kBlue + 1, 1);
    actsTruthPt->GetXaxis()->SetRangeUser(0., 2.);
    actsRecoPt->GetXaxis()->SetRangeUser(0., 2.);
    genfitTruthPt->GetXaxis()->SetRangeUser(0., 2.);
    genfitRecoPt->GetXaxis()->SetRangeUser(0., 2.);

    auto* trackPtCanvas = new TCanvas(
        "c_compare_track_pt", "Track pT: ACTS and GENFIT", 900, 700);
    prepareCanvas(trackPtCanvas);
    actsTruthPt->SetTitle("");
    actsTruthPt->GetXaxis()->SetTitle("Track p_{T} [GeV/c]");
    actsTruthPt->GetXaxis()->SetTitleOffset(1.05);
    actsTruthPt->GetYaxis()->SetTitleOffset(1.20);
    const double maximum = std::max(
        std::max(actsTruthPt->GetMaximum(), actsRecoPt->GetMaximum()),
        std::max(genfitTruthPt->GetMaximum(), genfitRecoPt->GetMaximum()));
    actsTruthPt->SetMaximum(1.20 * maximum);
    actsTruthPt->Draw("HIST");
    actsRecoPt->Draw("HIST SAME");
    genfitTruthPt->Draw("HIST SAME");
    genfitRecoPt->Draw("HIST SAME");

    auto* legend = new TLegend(0.59, 0.65, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(actsTruthPt.get(), "ACTS truth p_{T}", "l");
    legend->AddEntry(actsRecoPt.get(), "ACTS reco p_{T}", "l");
    legend->AddEntry(genfitTruthPt.get(), "GENFIT truth p_{T}", "l");
    legend->AddEntry(genfitRecoPt.get(), "GENFIT reco p_{T}", "l");
    legend->Draw();

    outputFile.cd();
    actsTruthPt->Write();
    actsRecoPt->Write();
    genfitTruthPt->Write();
    genfitRecoPt->Write();
    trackPtCanvas->Write();
    saveCanvas(trackPtCanvas, figureDirectory, "compare_track_pt");
  }
  else
  {
    std::cerr << "Skipping combined pT plot because one or more pT histograms are missing"
              << std::endl;
  }

  // Compare every common one-dimensional distribution produced by
  // TrackStateAnalysis.C. TProfiles are intentionally excluded: they encode
  // mean values rather than one-dimensional event/state distributions.
  std::set<std::string> histogramNames;
  TIter nextKey(actsFile.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(nextKey()))
  {
    histogramNames.insert(key->GetName());
  }

  std::size_t comparedHistograms = 0;
  for (const auto& histogramName : histogramNames)
  {
    if (histogramName == "h_truth_track_pt" || histogramName == "h_reco_track_pt")
    {
      continue;
    }

    std::unique_ptr<TH1> actsHistogram(
        cloneOneDimensionalHistogram(actsFile, histogramName, "acts_compare"));
    std::unique_ptr<TH1> genfitHistogram(
        cloneOneDimensionalHistogram(genfitFile, histogramName, "genfit_compare"));
    if (!actsHistogram || !genfitHistogram)
    {
      continue;
    }

    normalizeHistogram(actsHistogram.get(), normalize);
    normalizeHistogram(genfitHistogram.get(), normalize);
    styleHistogram(actsHistogram.get(), kRed + 1);
    styleHistogram(genfitHistogram.get(), kBlue + 1);
    actsHistogram->SetTitle("");
    actsHistogram->GetXaxis()->SetTitleOffset(1.05);
    actsHistogram->GetYaxis()->SetTitleOffset(1.20);
    actsHistogram->SetMaximum(
        1.20 * std::max(actsHistogram->GetMaximum(), genfitHistogram->GetMaximum()));

    const TString canvasName = TString::Format("c_compare_%s", histogramName.c_str());
    auto* canvas = new TCanvas(canvasName, canvasName, 900, 700);
    prepareCanvas(canvas);
    actsHistogram->Draw("HIST");
    genfitHistogram->Draw("HIST SAME");

    auto* legend = new TLegend(0.67, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(actsHistogram.get(), "ACTS", "l");
    legend->AddEntry(genfitHistogram.get(), "GENFIT", "l");
    legend->Draw();

    outputFile.cd();
    actsHistogram->Write();
    genfitHistogram->Write();
    canvas->Write();
    saveCanvas(canvas, figureDirectory, TString::Format("compare_%s", histogramName.c_str()));
    ++comparedHistograms;
  }

  outputFile.Close();
  std::cout << "Compared " << comparedHistograms
            << " common one-dimensional histograms" << std::endl;
  std::cout << "Wrote " << outputFileName << " and figures under "
            << figureDirectory << std::endl;
}
