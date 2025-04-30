void plot1D(TH1* h_noavg, TH1* h_avg, TChain* chain_noavg, TChain* chain_avg, TString name, TString output);
void plot2D(TH2* h_noavg, TH2* h_avg, TChain* chain_noavg, TChain* chain_avg, TString name, TString output);
void plot2D(TH2* h_avg, TChain* chain_avg, TString name, TString output);

void plot()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int nrun = 4;
  int runs[nrun] = {54966,54967,54968,54969};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain_noavg = new TChain("residualtree");
    chain_noavg->Add(Form("../noavgCorr/Reconstructed/%d/clusters_seeds_%d-*_resid.root",runs[k],runs[k]));
    TChain* chain_avg = new TChain("residualtree");
    chain_avg->Add(Form("../laminationMap/Reconstructed/%d/clusters_seeds_%d-*_resid.root",runs[k],runs[k]));

    TH1* h_quality_noavg = new TH1F("h_quality_noavg","quality;quality;Counts",30, 0,300);
    TH1* h_quality_avg = new TH1F("h_quality_avg","quality;quality;Counts",30, 0,300);
    plot1D(h_quality_noavg, h_quality_avg, chain_noavg, chain_avg, Form("quality"), Form("figure/%d_quality.pdf",runs[k]));

    TH1* h_tpcseedx_noavg = new TH1F("h_tpcseedx_noavg","tpcseedx;tpcseedx (cm);Counts",100, -5,5);
    TH1* h_tpcseedx_avg = new TH1F("h_tpcseedx_avg","tpcseedx;tpcseedx (cm);Counts",100, -5,5);
    plot1D(h_tpcseedx_noavg, h_tpcseedx_avg, chain_noavg, chain_avg, Form("tpcseedx"), Form("figure/%d_tpcseedx.pdf",runs[k]));

    TH1* h_tpcseedy_noavg = new TH1F("h_tpcseedy_noavg","tpcseedy;tpcseedy (cm);Counts",100, -5,5);
    TH1* h_tpcseedy_avg = new TH1F("h_tpcseedy_avg","tpcseedy;tpcseedy (cm);Counts",100, -5,5);
    plot1D(h_tpcseedy_noavg, h_tpcseedy_avg, chain_noavg, chain_avg, Form("tpcseedy"), Form("figure/%d_tpcseedy.pdf",runs[k]));


  }
}

void plot1D(TH1* h_noavg, TH1* h_avg, TChain* chain_noavg, TChain* chain_avg, TString name, TString output)
{
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), "");
    chain_avg->Draw((name + " >>" + h_avg->GetName()), "");
    //chain_noavg->Draw((name + " >>" + h_noavg->GetName()), "phi>-1.72816 && phi<-1.42667 && fabs(eta)<0.88");
    //chain_avg->Draw((name + " >>" + h_avg->GetName()), "phi>-1.72816 && phi<-1.42667 && fabs(eta)<0.88");
    cout<<"noavg nIntegral = "<<h_noavg->Integral()<<" , avg nIntegral = "<<h_avg->Integral()<<endl;
    h_noavg->Scale(h_avg->Integral() / h_noavg->Integral());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_noavg->SetLineColor(kBlack);
    h_avg->SetLineColor(kRed);
    double ymax = h_avg->GetMaximum() > h_noavg->GetMaximum() ? 1.1*(h_avg->GetMaximum()) : 1.1*(h_noavg->GetMaximum());
    h_avg->SetMaximum(ymax);
    h_avg->Draw("hist,e");
    h_noavg->Draw("hist,e,same");
    TLegend *legend = new TLegend(0.60, 0.75, 0.72, 0.9);
    legend->AddEntry(h_noavg, "No average correction", "l");
    legend->AddEntry(h_avg, "Lamination Map", "l");
    legend->SetTextSize(0.03);
    legend->Draw("same");
    can->Update();
    can->SaveAs(output);
    delete can;
}

void plot2D(TH2* h_noavg, TH2* h_avg, TChain* chain_noavg, TChain* chain_avg, TString name, TString output)
{
    TCut cut = "cluslayer<3";
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), cut, "colz");
    chain_avg->Draw((name + " >>" + h_avg->GetName()), cut, "colz");
    h_noavg->Scale(h_avg->Integral() / h_noavg->Integral());

    TCanvas *can_noavg = new TCanvas("can_noavg","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_noavg->Draw("colz");
    can_noavg->Update();
    can_noavg->SaveAs(output+Form("noavg.pdf"));

    TCanvas *can_avg = new TCanvas("can_avg","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_avg->Draw("colz");
    can_avg->Update();
    can_avg->SaveAs(output+Form("avg.pdf"));
    delete can_avg;
    delete can_noavg;
}

void plot2D(TH2* h_avg, TChain* chain_avg, TString name, TString output)
{
    //TCut cut = "cluslayer<3"; //mvtx
    //TCut cut = "cluslayer>2 && cluslayer<7"; //intt
    TCut cut = "cluslayer>6 && cluslayer<55"; //tpc
    //TCut cut = "cluslayer>54"; //tpot
    chain_avg->Draw((name + " >>" + h_avg->GetName()), cut, "colz");

    TCanvas *can_avg = new TCanvas("can_avg","",800,600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_avg->Draw("colz");
    can_avg->Update();
    can_avg->SaveAs(output+Form("avg.pdf"));
    delete can_avg;
}
