void plot1D(TH1* h, TChain* chain, TCut cut, TString name, TString output);

void plot_residual_1D()
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
    TChain* chain = new TChain("residualtree");
    chain->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment1/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));

    TH1* h_dlx_mvtx_l0 = new TH1F("h_dlx_mvtx_l0","statelx-cluslx in MVTX L0;statelx-cluslx (cm);Counts",100,-0.005,0.005);
    TH1* h_dlx_mvtx_l1 = new TH1F("h_dlx_mvtx_l1","statelx-cluslx in MVTX L1;statelx-cluslx (cm);Counts",100,-0.005,0.005);
    TH1* h_dlx_mvtx_l2 = new TH1F("h_dlx_mvtx_l2","statelx-cluslx in MVTX L2;statelx-cluslx (cm);Counts",100,-0.005,0.005);
    TH1* h_dlx_intt_l3 = new TH1F("h_dlx_intt_l3","statelx-cluslx in INTT L3;statelx-cluslx (cm);Counts",100,-0.02,0.02);
    TH1* h_dlx_intt_l4 = new TH1F("h_dlx_intt_l4","statelx-cluslx in INTT L4;statelx-cluslx (cm);Counts",100,-0.02,0.02);
    TH1* h_dlx_intt_l5 = new TH1F("h_dlx_intt_l5","statelx-cluslx in INTT L5;statelx-cluslx (cm);Counts",100,-0.02,0.02);
    TH1* h_dlx_intt_l6 = new TH1F("h_dlx_intt_l6","statelx-cluslx in INTT L5;statelx-cluslx (cm);Counts",100,-0.02,0.02);

    plot1D(h_dlx_mvtx_l0, chain, Form("cluslayer==0"), Form("statelx-cluslx"), Form("figure/%d_dlx_mvtx_l0.pdf",runs[k]));
    plot1D(h_dlx_mvtx_l1, chain, Form("cluslayer==1"), Form("statelx-cluslx"), Form("figure/%d_dlx_mvtx_l1.pdf",runs[k]));
    plot1D(h_dlx_mvtx_l2, chain, Form("cluslayer==2"), Form("statelx-cluslx"), Form("figure/%d_dlx_mvtx_l2.pdf",runs[k]));
    plot1D(h_dlx_intt_l3, chain, Form("cluslayer==3"), Form("statelx-cluslx"), Form("figure/%d_dlx_intt_l3.pdf",runs[k]));
    plot1D(h_dlx_intt_l4, chain, Form("cluslayer==4"), Form("statelx-cluslx"), Form("figure/%d_dlx_intt_l4.pdf",runs[k]));
    plot1D(h_dlx_intt_l5, chain, Form("cluslayer==5"), Form("statelx-cluslx"), Form("figure/%d_dlx_intt_l5.pdf",runs[k]));
    plot1D(h_dlx_intt_l6, chain, Form("cluslayer==6"), Form("statelx-cluslx"), Form("figure/%d_dlx_intt_l6.pdf",runs[k]));

  }
}

void plot1D(TH1* h, TChain* chain, TCut cut, TString name, TString output)
{
    chain->Draw((name + " >>" + h->GetName()), cut);
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h->SetLineColor(kBlack);
    h->Draw("hist");
    can->Update();
    can->SaveAs(output);
}
