void plot1D(TH1* h_100, TH1* h_1, TH1* h_1_blowupTPC, TChain* chain_100, TChain* chain_1, TChain* chain_1_blowupTPC, TString name, TString output);
void fit(TH1* h, float fitxmin, float fitxmax, TString output);

void plot_vertex_misalignment()
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
    TChain* chain_Misalignment100 = new TChain("vertextree");
    chain_Misalignment100->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment100/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Misalignment1 = new TChain("vertextree");
    chain_Misalignment1->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment1/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_Misalignment1_blowupTPC = new TChain("vertextree");
    chain_Misalignment1_blowupTPC->Add(Form("../staticCorrOn_scale1_BenMap_Misalignment1_blowupTPC/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));

    TH1* h_vx_Misalignment100 = new TH1F("h_vx_Misalignment100","vx;vx (cm);Counts",100,-0.2,0.2);
    TH1* h_vx_Misalignment1 = new TH1F("h_vx_Misalignment1","vx;vx (cm);Counts",100,-0.2,0.2);
    TH1* h_vx_Misalignment1_blowupTPC = new TH1F("h_vx_Misalignment1_blowupTPC","vx;vx (cm);Counts",100,-0.2,0.2);
    plot1D(h_vx_Misalignment100, h_vx_Misalignment1, h_vx_Misalignment1_blowupTPC, chain_Misalignment100, chain_Misalignment1, chain_Misalignment1_blowupTPC, Form("vx"), Form("figure/%d_vx_Misalignment.pdf",runs[k]));

    TH1* h_vy_Misalignment100 = new TH1F("h_vy_Misalignment100","vy;vy (cm);Counts",100,-0.2,0.5);
    TH1* h_vy_Misalignment1 = new TH1F("h_vy_Misalignment1","vy;vy (cm);Counts",100,-0.2,0.5);
    TH1* h_vy_Misalignment1_blowupTPC = new TH1F("h_vy_Misalignment1_blowupTPC","vy;vy (cm);Counts",100,-0.2,0.5);
    plot1D(h_vy_Misalignment100, h_vy_Misalignment1, h_vy_Misalignment1_blowupTPC, chain_Misalignment100, chain_Misalignment1, chain_Misalignment1_blowupTPC, Form("vy"), Form("figure/%d_vy_Misalignment.pdf",runs[k]));

    TH1* h_vz_Misalignment100 = new TH1F("h_vz_Misalignment100","vz;vz (cm);Counts",100,-20,20);
    TH1* h_vz_Misalignment1 = new TH1F("h_vz_Misalignment1","vz;vz (cm);Counts",100,-20,20);
    TH1* h_vz_Misalignment1_blowupTPC = new TH1F("h_vz_Misalignment1_blowupTPC","vz;vz (cm);Counts",100,-20,20);
    plot1D(h_vz_Misalignment100, h_vz_Misalignment1, h_vz_Misalignment1_blowupTPC, chain_Misalignment100, chain_Misalignment1, chain_Misalignment1_blowupTPC, Form("vz"), Form("figure/%d_vz_Misalignment.pdf",runs[k]));

    fit(h_vx_Misalignment1, -0.05,0.01, Form("figure/%d_vx_fit.pdf",runs[k]));
    fit(h_vy_Misalignment1, 0.08, 0.16, Form("figure/%d_vy_fit.pdf",runs[k]));
  }
}

void plot1D(TH1* h_100, TH1* h_1, TH1* h_1_blowupTPC, TChain* chain_100, TChain* chain_1, TChain* chain_1_blowupTPC, TString name, TString output)
{
    TCut cut="1";
    chain_100->Draw((name + " >>" + h_100->GetName()), cut);
    chain_1->Draw((name + " >>" + h_1->GetName()), cut);
    chain_1_blowupTPC->Draw((name + " >>" + h_1_blowupTPC->GetName()), cut);
    h_1->Scale(h_100->Integral() / h_1->Integral());
    h_1_blowupTPC->Scale(h_100->Integral() / h_1_blowupTPC->Integral());
    h_100->SetMaximum(1.5 * h_100->GetMaximum());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_100->SetLineColor(kBlack);
    h_1->SetLineColor(kRed);
    h_1_blowupTPC->SetLineColor(kBlue);
    h_100->Draw("hist");
    h_1->Draw("hist,same");
    h_1_blowupTPC->Draw("hist,same");
    TLegend *legend = new TLegend(0.55, 0.65, 0.9, 0.9);
    legend->AddEntry(h_100, "Si+TPOTMisalignment==100 , tpcMisalignment==1", "l");
    legend->AddEntry(h_1, "Si+TPOTMisalignment==1, tpcMisalignment==1", "l");
    legend->AddEntry(h_1_blowupTPC, "Si+TPOTMisalignment==1,tpcMisalignment==2", "l");
    legend->SetTextSize(0.02);
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}

void fit(TH1* h, float fitxmin, float fitxmax, TString output)
{
    TF1 *gausFit = new TF1("gausFit", "gaus", fitxmin, fitxmax);
    h->Fit(gausFit, "Q", "", fitxmin, fitxmax);

    Double_t mean = gausFit->GetParameter(1);
    Double_t mean_error = gausFit->GetParError(1);
    Double_t sigma = gausFit->GetParameter(2);
    Double_t sigma_error = gausFit->GetParError(2);

    TPaveText *pt = new TPaveText(.55, .70, .85, .85, "NDC");
    pt->SetFillColor(0);
    //pt->SetFillStyle(0);//transparent
    pt->SetLineColor(0);
    pt->SetBorderSize(0);
    pt->SetTextColor(kBlack);
    pt->AddText(Form("Mean = %.4f#pm%.4f mm", 10*mean, 10*mean_error));
    pt->AddText(Form("Sigma = %.4f#pm%.4f mm", 10*sigma, 10*sigma_error));

    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h->SetLineColor(kBlack);
    gausFit->SetLineColor(kRed);
    h->Draw("hist");
    gausFit->Draw("hist,same");
    pt->Draw("same");

    can->Update();
    can->SaveAs(output);
}
