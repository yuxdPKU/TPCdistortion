#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
//#include <RooDoubleCB.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooStats/SPlot.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TBranch.h>

using namespace RooFit;
using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 2)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

void Ks_fitter_nopull()
{

  string path, inputFile, inputFile2, plotTitle, saveName, xAxisTitle, branch;
  double minMass = 0.4;
  double maxMass = 0.6;
  int nBins = 50;

  path = "/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/staticCorrOn_scale1_BenMap_SiTPCMatchingWindow_newwindow/defaultWindow_newcode_ChangeZ/analysis/";
  inputFile = path + "all_pipi_reco.root";
  inputFile2 = path + "all_pipi_like_reco.root";

  branch = "K_S0_mass"; 
  plotTitle = "";
  saveName = path + "pipi_nopull"; 
  xAxisTitle = "m(#pi^{+}#pi^{-}) [GeV]"; 

  TCut masscut = "track_1_MVTX_nStates>=3 && track_2_MVTX_nStates>=3 && track_1_INTT_nStates>=2 && track_2_INTT_nStates>=2 && track_1_TPC_nStates>=30 && track_2_TPC_nStates>=30";

  /*
  stringstream cutStream;
  cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << "&& (track_1_MVTX_nStates>1 && track_2_MVTX_nStates>1 && track_1_track_2_DCA<0.04 && fabs(K_S0_decayLength)<0.05 && max(track_1_p,track_2_p)<0.7)";  // This draws the member "variable_name" in the tree
  TCut masscut = cutStream.str().c_str();
  cout<<masscut<<endl;
  */

  /*
   * Get files and data sets
   */

  TFile* dataFile = new TFile(inputFile.c_str());
  TTree* dataTree = (TTree*)dataFile->Get("DecayTree");

  RooRealVar mass(branch.c_str(), "mass", minMass, maxMass);

  //unbinned fit
  TFile* tempfile = new TFile("temp.root","recreate");
  TTree* dataTree_cut = (TTree*) dataTree->CopyTree(masscut);
  RooDataSet dataSet(branch.c_str(), "data", dataTree_cut, mass);
  //RooDataSet dataSet(branch.c_str(), "data", dataTree, mass);

  //binned fit
  //TH1F *h_data = new TH1F("h_data","", nBins, minMass, maxMass);
  //dataTree->Project("h_data",branch.c_str(),masscut);
  //RooDataHist dataSet(branch.c_str(),"data",mass,h_data);

  /*
   * Signal Model
   */


  RooRealVar  meanGauss("K_S0_mean", "mean", 0.50, 0.40, 0.60);
  RooRealVar  sigmaGauss("sigmaGauss", "sigmaGauss", 0.009, 0.003, 0.03);
  RooGaussian Gauss("Gauss", "Voigt Profile", mass, meanGauss, sigmaGauss);

  //end of K_S0 Signal model

  RooRealVar nSig("nSig", "nSig", 0.9*dataSet.sumEntries(), 0, 1e8);

  /*
   * Background Model
   */

  RooRealVar p1("p1","coeff #1", 0, -100., 100.);
  RooRealVar p2("p2","coeff #2", 0, -100., 100.);
  RooRealVar p3("p3","coeff #3", 0, -100., 100.);
  RooRealVar p4("p4","coeff #4", 10, -100., 100.);
  RooRealVar p5("p5","coeff #5", 10, -100., 100.);
  RooRealVar p6("p6","coeff #6", 10, -100., 100.);
  RooRealVar p7("p7","coeff #7", 10, -100., 100.);
  RooRealVar p8("p8","coeff #8", 10, -100., 100.);
  RooPolynomial background("background","background", mass, RooArgList(p1, p2, p3));

  RooRealVar nBkg("nBkg", "nBkg", 0.1*dataSet.sumEntries(), 0, 1e8);

  /*
   * Fitting to the data
   */

  RooAddPdf model("model", "model", RooArgList(Gauss, background), RooArgList(nSig, nBkg));
  RooFitResult *m_fitres(0);
  m_fitres = model.fitTo(dataSet, Save(kTRUE), Extended(kTRUE));
  m_fitres->Print("v");

  //Print signal yield
  string fittedSigYield = to_string_with_precision(nSig.getValV(), 0);
  string fittedSigError = to_string_with_precision(nSig.getError(), 0);
  int nEvents = dataTree->GetEntries(masscut);
  cout << "*\n*\n* nEvents = " << nEvents << "\n* nSig = " << fittedSigYield << " +/- " << fittedSigError << "\n*\n*" <<endl;

  RooPlot* frame = mass.frame(Title(plotTitle.c_str()));       //creating the frame

  RooBinning bins(minMass, maxMass);
  bins.addUniform(nBins, minMass, maxMass);

  dataSet.plotOn(frame, Binning(bins), XErrorSize(0), DataError(RooAbsData::SumW2));          //plotting the raw unbinned data on the fram
  dataSet.plotOn(frame, Binning(bins), XErrorSize(0), DataError(RooAbsData::SumW2));          //plotting the raw unbinned data on the fram
  model.plotOn(frame, Components(background), LineColor(kGray), DrawOption("F"), FillColor(kGray));
  model.plotOn(frame, Components(RooArgSet(Gauss, background)), LineColor(kAzure+8), DrawOption("F"), FillColor(kAzure+8), MoveToBack());
  model.plotOn(frame, LineColor(kBlack));        //plotting the fit onto the frame
  dataSet.plotOn(frame, DrawOption("PE1"), Binning(bins),XErrorSize(0), DataError(RooAbsData::SumW2));          //plotting the raw unbinned data on the fram

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  //RooHist* pull = frame->pullHist();
  //RooPlot* frame2 = mass.frame(Title(""));
  //frame2->addPlotable(pull,"PE1");
  int maxDigits = 3;
  TGaxis::SetMaxDigits(maxDigits);

  TCanvas* c = new TCanvas("massFitCanvas", "massFitCanvas",800, 800);  

  TPad mainPad("mainPad", "mainPad", 0., 0.3, 1., 1.);
  //mainPad.SetBottomMargin(0);
  mainPad.Draw();
  frame->SetMarkerStyle(kCircle);
  frame->SetMarkerSize(0.02);
  frame->SetLineWidth(1);
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetYaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetLabelFont(42);
  //frame->SetMinimum(85);
  //frame->SetMaximum(200);
  //frame->GetYaxis()->SetRangeUser(405, 1700); //For 2 GeV
  frame->GetXaxis()->SetNdivisions(5, 6, 0, kTRUE);
  float binWidth = 1000.*(maxMass - minMass) / nBins;
  //string yAxisTitle = "Candidates / (" + to_string_with_precision(binWidth) + " MeV)";
  //frame->GetYaxis()->SetTitle(yAxisTitle.c_str());
  TString yAxisTitle = Form("Candidates / (%.0f MeV)",binWidth);
  frame->GetYaxis()->SetTitle(yAxisTitle);
  //frame->GetYaxis()->SetTitleOffset(1.0);
  frame->Draw();
  //gPad->RedrawAxis();

  string fittedSigYieldString = "Yield = " + to_string_with_precision(nSig.getValV(), 0) + "#kern[-0.3]{ #pm} " + fittedSigError;

  string fittedK_S0_mean = to_string_with_precision(meanGauss.getValV() * 1000, 0);
  string fittedK_S0_meanError = to_string_with_precision(meanGauss.getError() * 1000, 0);

  string fittedK_S0_meanString = "#mu = " + fittedK_S0_mean + "#kern[-0.4]{ #pm} " + fittedK_S0_meanError + " MeV";

  string fittedsigmaGauss = to_string_with_precision(sigmaGauss.getValV() * 1000, 1);
  string fittedsigmaGaussError = to_string_with_precision(sigmaGauss.getError() * 1000, 1);

  string fittedsigmaGaussString = "#sigma = " + fittedsigmaGauss + "#kern[-0.4]{ #pm} " + fittedsigmaGaussError + " MeV";

  mass.setRange("sigwindow",0.47,0.51);

  RooAbsReal *sigInte = Gauss.createIntegral(mass,NormSet(mass),Range("sigwindow"));
  RooAbsReal *bkgInte = background.createIntegral(mass,NormSet(mass),Range("sigwindow"));

  Double_t fintsig = sigInte->getVal();
  Double_t fintbkg = bkgInte->getVal();

  Double_t nSigint = (nSig.getValV())*fintsig;
  Double_t nBkgint = (nBkg.getValV())*fintbkg;
  Double_t purity = nSigint / (nSigint+nBkgint);
  std::cout<<"[0.47, 0.51]: nsignal = "<<nSigint<<" , nbkg = "<<nBkgint<<" , purity = "<<purity<<std::endl;

  TPaveText *pt;
  pt = new TPaveText(0.60,0.66,0.90,0.90, "NDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  //TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Preliminary");
  TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt_LaTex = pt->AddText("#it{p}+#it{p} #sqrt{s} = 200 GeV");
  //pt_LaTex = pt->AddText("~1 hour");
  pt_LaTex = pt->AddText(fittedK_S0_meanString.c_str());
  pt_LaTex = pt->AddText(fittedsigmaGaussString.c_str());
  pt_LaTex = pt->AddText(fittedSigYieldString.c_str());
  pt->SetBorderSize(0);
  pt->SetTextSize(0.04);
  pt->DrawClone();

  //date timestamp
  TPaveText *ptDate;
  ptDate = new TPaveText(0.67,0.92,1.0,1.00, "NDC");
  ptDate->SetFillColor(0);
  ptDate->SetFillStyle(0);
  ptDate->SetTextFont(42);
  TText *pt_LaTexDate = ptDate->AddText("05/22/2025");
  ptDate->SetBorderSize(0);
  ptDate->SetTextSize(0.05);
  ptDate->Draw();
  gPad->Modified();

  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.05);
  TLegend *legend = new TLegend(0.20,0.66,0.49,0.90);
  legend->AddEntry(frame->findObject("h_K_S0_mass"),"Data","PE2");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]"),"Fit","L");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[Gauss,background]"),"K_{S}^{0}#rightarrow#pi^{+}#pi^{-}","f");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[background]"),"Comb. Bkg.","f");
  
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->Draw();

  gPad->Modified();

  vector<string> extensions = {".C", ".pdf", ".png"};
  for (auto &extension : extensions)
  {
    string outputFile = saveName + extension;
    c->SaveAs(outputFile.c_str());
  }
 
  gSystem->Exec("rm -f temp*.root");
}
