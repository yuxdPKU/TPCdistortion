void plot1D(TH1* h_noavg, TH1* h_Ben, TH1* h_Xudong, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TString name, TString output);
void plot2D(TH2* h_noavg, TH2* h_Ben, TH2* h_Xudong, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TString name, TString output);
void plot2D(TH2* h_Ben, TChain* chain_Ben, TString name, TString output);

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
    TChain* chain_noavg = new TChain("residualtree");
    chain_noavg->Add(Form("../staticCorrOn_scale1_noavgCorr_3/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    //chain_noavg->Add(Form("../staticCorrOn_scale1_noavgCorr_2/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    TChain* chain_Ben = new TChain("residualtree");
    chain_Ben->Add(Form("../staticCorrOn_scale1_BenMap_3/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    //chain_Ben->Add(Form("../staticCorrOn_scale1_BenMap_2/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    TChain* chain_Xudong = new TChain("residualtree");
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_5/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_5/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_6/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_4/root/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_4/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_3/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));
    //chain_Xudong->Add(Form("../staticCorrOn_scale1_XudongMap_2/Reconstructed/%d/clusters_seeds_%d-0.root_resid.root",runs[k],runs[k]));

    TH1* h_quality_noavg = new TH1F("h_quality_noavg","quality;quality;Counts",90, 0,300);
    TH1* h_quality_Ben = new TH1F("h_quality_Ben","quality;quality;Counts",90, 0,300);
    TH1* h_quality_Xudong = new TH1F("h_quality_Xudong","quality;quality;Counts",90, 0,300);
    plot1D(h_quality_noavg, h_quality_Ben, h_quality_Xudong, chain_noavg, chain_Ben, chain_Xudong, Form("quality"), Form("figure/%d_quality.pdf",runs[k]));

    /*
    TH2* h_quality_phi_noavg = new TH2F("h_quality_phi_noavg","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,-M_PI,M_PI);
    TH2* h_quality_phi_Ben = new TH2F("h_quality_phi_Ben","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,-M_PI,M_PI);
    TH2* h_quality_phi_Xudong = new TH2F("h_quality_phi_Xudong","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,-M_PI,M_PI);
    //TH2* h_quality_phi_noavg = new TH2F("h_quality_phi_noavg","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,1,2);
    //TH2* h_quality_phi_Ben = new TH2F("h_quality_phi_Ben","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,1,2);
    //TH2* h_quality_phi_Xudong = new TH2F("h_quality_phi_Xudong","quality vs. phi;quality;phi (rad);Counts",100, 0,300,100,1,2);
    plot2D(h_quality_phi_noavg, h_quality_phi_Ben, h_quality_phi_Xudong, chain_noavg, chain_Ben, chain_Xudong, Form("phi:quality"), Form("figure/%d_quality_phi.pdf",runs[k]));
    */

    //TH2* h_dlx_phi_Ben = new TH2F("h_dlx_phi_Ben","dlx vs. phi;phi (rad);dlx;Counts",500,-M_PI,M_PI,500,-0.02,0.02); // mvtx
    //TH2* h_dlx_phi_Ben = new TH2F("h_dlx_phi_Ben","dlx vs. phi;phi (rad);dlx;Counts",500,-M_PI,M_PI,500,-0.3,0.3); // intt
    //TH2* h_dlx_phi_Ben = new TH2F("h_dlx_phi_Ben","dlx vs. phi;phi (rad);dlx;Counts",500,-M_PI,M_PI,500,-0.3,0.3); // tpc
    //h_dlx_phi_Ben->GetXaxis()->SetNoExponent(true);
    //plot2D(h_dlx_phi_Ben, chain_Ben, Form("statelx-cluslx:atan2(clusgy,clusgx)"), Form("figure/%d_dlx_phi.pdf",runs[k]));

    //TH2* h_cluslx_phi_Ben = new TH2F("h_cluslx_phi_Ben","cluster local x vs. phi;phi (rad);cluster local x (cm);Counts",500,-M_PI,M_PI,500,-1,1); // mvtx
    //TH2* h_cluslx_phi_Ben = new TH2F("h_cluslx_phi_Ben","cluster local x vs. phi;phi (rad);cluster local x (cm);Counts",500,-M_PI,M_PI,500,-2,2); // intt
    //TH2* h_cluslx_phi_Ben = new TH2F("h_cluslx_phi_Ben","cluster local x vs. phi;phi (rad);cluster local x (cm);Counts",500,-M_PI,M_PI,500,-2,2); // tpc
    //h_cluslx_phi_Ben->GetXaxis()->SetNoExponent(true);
    //plot2D(h_cluslx_phi_Ben, chain_Ben, Form("cluslx:atan2(clusgy,clusgx)"), Form("figure/%d_cluslx_phi.pdf",runs[k]));

    //TH2* h_statelx_phi_Ben = new TH2F("h_statelx_phi_Ben","state local x vs. phi;phi (rad);state local x (cm);Counts",500,-M_PI,M_PI,500,-1,1); // mvtx
    //TH2* h_statelx_phi_Ben = new TH2F("h_statelx_phi_Ben","state local x vs. phi;phi (rad);state local x (cm);Counts",500,-M_PI,M_PI,500,-2,2); // intt
    //TH2* h_statelx_phi_Ben = new TH2F("h_statelx_phi_Ben","state local x vs. phi;phi (rad);state local x (cm);Counts",500,-M_PI,M_PI,500,-2,2); // tpc
    //h_statelx_phi_Ben->GetXaxis()->SetNoExponent(true);
    //plot2D(h_statelx_phi_Ben, chain_Ben, Form("statelx:atan2(stategy,stategx)"), Form("figure/%d_statelx_phi.pdf",runs[k]));

  }
}

void plot1D(TH1* h_noavg, TH1* h_Ben, TH1* h_Xudong, TChain* chain_noavg, TChain* chain_Ben, TChain* chain_Xudong, TString name, TString output)
{
    chain_noavg->Draw((name + " >>" + h_noavg->GetName()), "");
    chain_Ben->Draw((name + " >>" + h_Ben->GetName()), "");
    chain_Xudong->Draw((name + " >>" + h_Xudong->GetName()), "");
    //chain_noavg->Draw((name + " >>" + h_noavg->GetName()), "phi>-1.72816 && phi<-1.42667 && fabs(eta)<0.88");
    //chain_Ben->Draw((name + " >>" + h_Ben->GetName()), "phi>-1.72816 && phi<-1.42667 && fabs(eta)<0.88");
    //chain_Xudong->Draw((name + " >>" + h_Xudong->GetName()), "phi>-1.72816 && phi<-1.42667 && fabs(eta)<0.88");
    h_noavg->Scale(h_Ben->Integral() / h_noavg->Integral());
    h_Xudong->Scale(h_Ben->Integral() / h_Xudong->Integral());
    TCanvas *can = new TCanvas("can","",800,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h_noavg->SetLineColor(kBlack);
    h_Ben->SetLineColor(kRed);
    h_Xudong->SetLineColor(kBlue);
    h_Ben->Draw("hist");
    h_noavg->Draw("e,same");
    h_Xudong->Draw("hist,same");
    TLegend *legend = new TLegend(0.45, 0.75, 0.72, 0.9);
    legend->AddEntry(h_noavg, "No average correction", "l");
    legend->AddEntry(h_Ben, "Ben's Map", "l");
    legend->AddEntry(h_Xudong, "Xudong's Map", "l");
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
