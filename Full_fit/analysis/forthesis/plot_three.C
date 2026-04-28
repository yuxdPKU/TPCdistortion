void plot1D(TH1* h_noavg, TH1* h_Ben, TH1* h_Xudong, TH1* h_DC, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TChain* chain_DC, TString name, TString output);
void plot2D(TH2* h_noavg, TH2* h_Ben, TH2* h_Xudong, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TString name, TString output);
void plot2D(TH2* h_Ben, TChain* chain_Ben, TString name, TString output);

void plot_three()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 2;
  int runs[nrun] = {79507, 79516};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain_noavg = new TChain("residualtree");
    chain_noavg->Add(Form("../../staticCorrOn_scale1_noavgCorr_5/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Ben = new TChain("residualtree");
    chain_Ben->Add(Form("../../staticCorrOn_scale1_BenMap_5/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Xudong = new TChain("residualtree");
    chain_Xudong->Add(Form("../../staticCorrOn_scale1_XudongMap_7/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_DC = new TChain("residualtree");
    chain_DC->Add(Form("../../staticCorrOn_scale1_DCMap_1/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));

    TH1* h_quality_noavg = new TH1F("h_quality_noavg","quality;quality;Counts",100,0,10);
    TH1* h_quality_Ben = new TH1F("h_quality_Ben","quality;quality;Counts",100,0,10);
    TH1* h_quality_Xudong = new TH1F("h_quality_Xudong","quality;quality;Counts",100,0,10);
    TH1* h_quality_DC = new TH1F("h_quality_DC","quality;quality;Counts",100,0,10);
    plot1D(h_quality_noavg, h_quality_Ben, h_quality_Xudong, h_quality_DC, chain_noavg, chain_Ben, chain_Xudong, chain_DC, Form("quality"), Form("%d_quality.pdf",runs[k]));
  }
}

void plot1D(TH1* h_noavg, TH1* h_Ben, TH1* h_Xudong, TH1* h_DC, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TChain* chain_DC, TString name, TString output)
{
    TCut cut = "m_nmmsstate>0";
    //TCut cut = "1";
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), cut);
    chain_Ben->Draw((name + " >>" + h_Ben->GetName()), cut);
    chain_Xudong->Draw((name + " >>" + h_Xudong->GetName()), cut);
    chain_DC->Draw((name + " >>" + h_DC->GetName()), cut);
    h_Ben->Scale(h_noavg->Integral() / h_Ben->Integral());
    h_Xudong->Scale(h_noavg->Integral() / h_Xudong->Integral());
    h_DC->Scale(h_noavg->Integral() / h_DC->Integral());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_noavg->SetLineColor(kBlack);
    h_Ben->SetLineColor(kRed);
    h_Xudong->SetLineColor(kBlue);
    h_DC->SetLineColor(kViolet);
    h_Xudong->Draw("hist");
    h_Ben->Draw("hist,same");
    h_noavg->Draw("hist,same");
    h_DV->Draw("hist,same");
    TLegend *legend = new TLegend(0.45, 0.70, 0.72, 0.9);
    legend->SetHeader("nTPOTstate>0");
    legend->AddEntry(h_noavg, "No average correction", "l");
    legend->AddEntry(h_Ben, "Lamination correction", "l");
    legend->AddEntry(h_Xudong, "Si-TPOT correction", "l");
    legend->AddEntry(h_DC, "DC correction", "l");
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}

void plot2D(TH2* h_noavg, TH2* h_Ben, TH2* h_Xudong, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TString name, TString output)
{
    TCut cut = "cluslayer<3";
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), cut, "colz");
    chain_Ben->Draw((name + " >>" + h_Ben->GetName()), cut, "colz");
    chain_Xudong->Draw((name + " >>" + h_Xudong->GetName()), cut, "colz");
    h_noavg->Scale(h_Ben->Integral() / h_noavg->Integral());
    h_Xudong->Scale(h_Ben->Integral() / h_Xudong->Integral());

    TCanvas *can_noavg = new TCanvas("can_noavg","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_noavg->Draw("colz");
    can_noavg->Update();
    can_noavg->SaveAs(output+Form("noavg.pdf"));

    TCanvas *can_Ben = new TCanvas("can_Ben","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Ben->Draw("colz");
    can_Ben->Update();
    can_Ben->SaveAs(output+Form("Ben.pdf"));

    TCanvas *can_Xudong = new TCanvas("can_Xudong","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Xudong->Draw("colz");
    can_Xudong->Update();
    can_Xudong->SaveAs(output+Form("Xudong.pdf"));
}

void plot2D(TH2* h_Ben, TChain* chain_Ben, TString name, TString output)
{
    //TCut cut = "cluslayer<3"; //mvtx
    //TCut cut = "cluslayer>2 && cluslayer<7"; //intt
    TCut cut = "cluslayer>6 && cluslayer<55"; //tpc
    //TCut cut = "cluslayer>54"; //tpot
    chain_Ben->Draw((name + " >>" + h_Ben->GetName()), cut, "colz");

    TCanvas *can_Ben = new TCanvas("can_Ben","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Ben->Draw("colz");
    can_Ben->Update();
    can_Ben->SaveAs(output+Form("Ben.pdf"));
}
