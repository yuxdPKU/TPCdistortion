void plot1D(TH1* h_100, TH1* h_1, TH1* h_1_tpc, TChain* chain_100, TChain* chain_1, TChain* chain_1_tpc, TString name, TString output);

void plot_dca_misalignment()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 1;
  //int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
  //int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
  //int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
  //int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
  int mbdrates[nrun] = {250};
  int runs[nrun] = {53534};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain_Misalignment100 = new TChain("residualtree");
    chain_Misalignment100->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment100/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Misalignment1 = new TChain("residualtree");
    chain_Misalignment1->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment1/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Misalignment1_blowupTPC = new TChain("residualtree");
    chain_Misalignment1_blowupTPC->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment1_blowupTPC/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));

    TH1* h_dcaxy_Misalignment100 = new TH1F("h_dcaxy_Misalignment100","dcaxy;dcaxy (cm);Counts",100,-1,1);
    TH1* h_dcaxy_Misalignment1 = new TH1F("h_dcaxy_Misalignment1","dcaxy;dcaxy (cm);Counts",100,-1,1);
    TH1* h_dcaxy_Misalignment1_blowupTPC = new TH1F("h_dcaxy_Misalignment1_blowupTPC","dcaxy;dcaxy (cm);Counts",100,-1,1);
    plot1D(h_dcaxy_Misalignment100, h_dcaxy_Misalignment1, h_dcaxy_Misalignment1_blowupTPC, chain_Misalignment100, chain_Misalignment1, chain_Misalignment1, Form("dcaxy"), Form("figure/%d_dcaxy_Misalignment.pdf",runs[k]));

    TH1* h_dcaz_Misalignment100 = new TH1F("h_dcaz_Misalignment100","dcaz;dcaz (cm);Counts",100,-50,50);
    TH1* h_dcaz_Misalignment1 = new TH1F("h_dcaz_Misalignment1","dcaz;dcaz (cm);Counts",100,-50,50);
    TH1* h_dcaz_Misalignment1_blowupTPC = new TH1F("h_dcaz_Misalignment1_blowupTPC","dcaz;dcaz (cm);Counts",100,-50,50);
    plot1D(h_dcaz_Misalignment100, h_dcaz_Misalignment1, h_dcaz_Misalignment1_blowupTPC, chain_Misalignment100, chain_Misalignment1, chain_Misalignment1, Form("dcaz"), Form("figure/%d_dcaz_Misalignment.pdf",runs[k]));

  }
}

void plot1D(TH1* h_100, TH1* h_1, TH1* h_1_tpc, TChain* chain_100, TChain* chain_1, TChain* chain_1_tpc, TString name, TString output)
{
    TCut cut="ntpc>20 && nmaps>=2 && nintt>=1 && pt>=0.5";
    chain_100->Draw((name + " >>" + h_100->GetName()), cut);
    chain_1->Draw((name + " >>" + h_1->GetName()), cut);
    chain_1_tpc->Draw((name + " >>" + h_1_tpc->GetName()), cut);
    h_1->Scale(h_100->Integral() / h_1->Integral());
    h_1_tpc->Scale(h_100->Integral() / h_1_tpc->Integral());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_100->SetLineColor(kBlack);
    h_1->SetLineColor(kRed);
    h_1_tpc->SetLineColor(kBlue);
    h_100->Draw("hist");
    h_1->Draw("hist,same");
    h_1_tpc->Draw("hist,same");
    TLegend *legend = new TLegend(0.55, 0.65, 0.9, 0.9);
    legend->AddEntry(h_100, "Si+TPOTMisalignment==100 , tpcMisalignment==1", "l");
    legend->AddEntry(h_1, "Si+TPOTMisalignment==1, tpcMisalignment==1", "l");
    legend->AddEntry(h_1_tpc, "Si+TPOTMisalignment==1,tpcMisalignment==2", "l");
    legend->SetTextSize(0.02);
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}
