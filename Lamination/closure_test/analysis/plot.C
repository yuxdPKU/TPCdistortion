void plot1D(int run, TH1* h_noavg, TH1* h_Ben, TChain* chain_noavg, TChain* chain_avg, TString name, TString output);

void plot()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 24;
  int runs[nrun] = {53018,53046,53079,53080,53081,53194,53195,53196,53197,53494,53513,53517,53530,53531,53532,53571,53578,53579,53580,53581,53586,53587,53590,53686};
  //const int nrun = 1;
  //int runs[nrun] = {53877};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain_noavg = new TChain("residualtree");
    chain_noavg->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Lamination/closure_test/disable_avg_distortion/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    TChain* chain_avg = new TChain("residualtree");
    chain_avg->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Lamination/closure_test/enable_avg_distortion/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));

    //TH1* h_quality_noavg = new TH1F("h_quality_noavg","quality;quality;Counts",90, 0,300);
    //TH1* h_quality_avg = new TH1F("h_quality_avg","quality;quality;Counts",90, 0,300);
    //plot1D(runs[k],h_quality_noavg, h_quality_avg, chain_noavg, chain_avg, Form("quality"), Form("figure/%d_quality.pdf",runs[k]));

    //TH1* h_dx_noavg = new TH1F("h_dx_noavg",";tpcseedx-silseedx;Counts",100,-10,10);
    //TH1* h_dx_avg = new TH1F("h_dx_avg",";tpcseedx-silseedx;Counts",100,-10,10);
    //plot1D(runs[k],h_dx_noavg, h_dx_avg, chain_noavg, chain_avg, Form("tpcseedx-silseedx"), Form("figure/%d_dx.pdf",runs[k]));

    //TH1* h_dy_noavg = new TH1F("h_dy_noavg",";tpcseedy-silseedy;Counts",100,-10,10);
    //TH1* h_dy_avg = new TH1F("h_dy_avg",";tpcseedy-silseedy;Counts",100,-10,10);
    //plot1D(runs[k],h_dy_noavg, h_dy_avg, chain_noavg, chain_avg, Form("tpcseedy-silseedy"), Form("figure/%d_dy.pdf",runs[k]));

    TH1* h_dz_noavg = new TH1F("h_dz_noavg",";tpcseedz-silseedz;Counts",100,-200,200);
    TH1* h_dz_avg = new TH1F("h_dz_avg",";tpcseedz-silseedz;Counts",100,-200,200);
    plot1D(runs[k],h_dz_noavg, h_dz_avg, chain_noavg, chain_avg, Form("tpcseedz-silseedz"), Form("figure/%d_dz.pdf",runs[k]));

    TH1* h_dphi_noavg = new TH1F("h_dphi_noavg",";tpcseedphi-silseedphi;Counts",100,-0.3,0.3);
    TH1* h_dphi_avg = new TH1F("h_dphi_avg",";tpcseedphi-silseedphi;Counts",100,-0.3,0.3);
    plot1D(runs[k],h_dphi_noavg, h_dphi_avg, chain_noavg, chain_avg, Form("tpcseedphi-silseedphi"), Form("figure/%d_dphi.pdf",runs[k]));

    TH1* h_deta_noavg = new TH1F("h_deta_noavg",";tpcseedeta-silseedeta;Counts",100,-0.3,0.3);
    TH1* h_deta_avg = new TH1F("h_deta_avg",";tpcseedeta-silseedeta;Counts",100,-0.3,0.3);
    plot1D(runs[k],h_deta_noavg, h_deta_avg, chain_noavg, chain_avg, Form("tpcseedeta-silseedeta"), Form("figure/%d_deta.pdf",runs[k]));

  }
}

void plot1D(int run, TH1* h_noavg, TH1* h_avg, TChain* chain_noavg, TChain* chain_avg, TString name, TString output)
{
    //TCut cut = "nmaps>=3 && nintt>=2 && ntpc>=30";
    TCut cut = "ntpc>=30";
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), cut);
    chain_avg->Draw((name + " >>" + h_avg->GetName()), cut);
    //h_noavg->Scale(h_avg->Integral() / h_noavg->Integral());
    int ymax = h_noavg->GetMaximum() > h_avg->GetMaximum() ? 1.1 *h_noavg->GetMaximum() : 1.1 * h_avg->GetMaximum();
    h_noavg->SetMaximum(ymax);
    h_noavg->SetMinimum(0);
    h_avg->SetMaximum(ymax);
    h_avg->SetMinimum(0);
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    h_noavg->SetLineColor(kBlack);
    h_avg->SetLineColor(kRed);
    h_avg->Draw("hist");
    h_noavg->Draw("hist,same");
    //TLegend *legend = new TLegend(0.40, 0.70, 0.70, 0.90);
    TLegend *legend = new TLegend(0.60, 0.70, 0.90, 0.90);
    legend->SetHeader(Form("Run %d, 10 thousand events",run));
    legend->AddEntry(h_noavg, "W/o average correction", "l");
    legend->AddEntry(h_avg, "W/ average correction", "l");
    legend->SetTextSize(0.03);
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}
