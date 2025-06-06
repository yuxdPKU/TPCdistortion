void plot1D(TH1* h_nofluc, TH1* h_fluc, TChain* chain_nofluc, TChain* chain_fluc, TString name, TString output);
void plot2D(TH2* h_nofluc, TH2* h_fluc, TChain* chain_nofluc, TChain* chain_fluc, TString name, TString output);
void plot2D(TH2* h_nofluc, TChain* chain_nofluc, TString name, TString output);

void plot()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 1;
  //int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
  //int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
  //int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
  //int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
  int mbdrates[nrun] = {400};
  int runs[nrun] = {53877};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain_nofluc = new TChain("residualtree");
    chain_nofluc->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/staticCorrOn_scale1_BenMap_2/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    TChain* chain_fluc = new TChain("residualtree");
    chain_fluc->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/fluctuation/Reconstructed/%d/clusters_seeds_%d-*_resid.root",runs[k],runs[k]));

    TH1* h_quality_nofluc = new TH1F("h_quality_nofluc","quality;quality;Counts",30, 0,300);
    TH1* h_quality_fluc = new TH1F("h_quality_fluc","quality;quality;Counts",30, 0,300);
    plot1D(h_quality_nofluc, h_quality_fluc, chain_nofluc, chain_fluc, Form("quality"), Form("figure/%d_quality.pdf",runs[k]));
  }
}

void plot1D(TH1* h_nofluc, TH1* h_fluc, TChain* chain_nofluc, TChain* chain_fluc, TString name, TString output)
{
    chain_nofluc->Draw((name + " >>" + h_nofluc->GetName()), "");
    chain_fluc->Draw((name + " >>" + h_fluc->GetName()), "");
    h_fluc->Scale(h_nofluc->Integral() / h_fluc->Integral());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_nofluc->SetLineColor(kRed);
    h_fluc->SetLineColor(kBlue);
    h_nofluc->Draw("hist,e");
    h_fluc->Draw("hist,e,same");
    TLegend *legend = new TLegend(0.45, 0.75, 0.72, 0.9);
    legend->AddEntry(h_nofluc, "w/o fluctuation corr", "l");
    legend->AddEntry(h_fluc, "w/ fluctuation corr", "l");
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}

void plot2D(TH2* h_nofluc, TH2* h_fluc, TChain* chain_nofluc, TChain* chain_fluc, TString name, TString output)
{
    TCut cut = "cluslayer<3";
    chain_nofluc->Draw((name + " >>" + h_nofluc->GetName()), cut, "colz");
    chain_fluc->Draw((name + " >>" + h_fluc->GetName()), cut, "colz");
    h_fluc->Scale(h_nofluc->Integral() / h_fluc->Integral());

    TCanvas *can_nofluc = new TCanvas("can_nofluc","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_nofluc->Draw("colz");
    can_nofluc->Update();
    can_nofluc->SaveAs(output+Form("nofluc.pdf"));

    TCanvas *can_fluc = new TCanvas("can_fluc","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_fluc->Draw("colz");
    can_fluc->Update();
    can_fluc->SaveAs(output+Form("fluc.pdf"));
}

void plot2D(TH2* h_nofluc, TChain* chain_nofluc, TString name, TString output)
{
    //TCut cut = "cluslayer<3"; //mvtx
    //TCut cut = "cluslayer>2 && cluslayer<7"; //intt
    TCut cut = "cluslayer>6 && cluslayer<55"; //tpc
    //TCut cut = "cluslayer>54"; //tpot
    chain_nofluc->Draw((name + " >>" + h_nofluc->GetName()), cut, "colz");

    TCanvas *can_nofluc = new TCanvas("can_nofluc","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_nofluc->Draw("colz");
    can_nofluc->Update();
    can_nofluc->SaveAs(output+Form("nofluc.pdf"));
}
