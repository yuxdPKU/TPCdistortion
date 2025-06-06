void plot1D(TH1* h, TChain* chain, TString name, TCut cut, TString output, bool logy=0);
void plot1D_QoverPt(TH1* h[], TChain* chain, TString name, TCut cut, TString output, double fit_xmin, double fit_xmax, int nsigmawindow=3, bool logy=0, bool doexpfit=0);
void plot2D(TH2* h, TChain* chain, TString name, TCut cut, TString output, bool logz=0);
void plot2D(TH2* h_Ben, TChain* chain_Ben, TString name, TString output);

const int nQoverPtbin=44;
double QoverPt_min[nQoverPtbin]={-8.00,-6.02,-5.53,-5.15,-4.82,-4.51,-4.25,-4.00,-3.76,-3.53,-3.32,-3.11,-2.92,-2.72,-2.52,-2.32,-2.11,-1.90,-1.66,-1.38,-1.00,-0.33,0.00,0.20,0.33,0.81,1.11,1.34,1.55,1.75,1.94,2.13,2.32,2.52,2.72,2.94,3.17,3.42,3.68,3.97,4.27,4.65,5.08,5.67};
double QoverPt_max[nQoverPtbin]={-6.02,-5.53,-5.15,-4.82,-4.51,-4.25,-4.00,-3.76,-3.53,-3.32,-3.11,-2.92,-2.72,-2.52,-2.32,-2.11,-1.90,-1.66,-1.38,-1.00,-0.33,0.00,0.20,0.33,0.81,1.11,1.34,1.55,1.75,1.94,2.13,2.32,2.52,2.72,2.94,3.17,3.42,3.68,3.97,4.27,4.65,5.08,5.67,8.00};
double QoverPt_mid[nQoverPtbin];
double QoverPt_interval[nQoverPtbin];

void plot()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  for (int i=0; i<nQoverPtbin; i++)
  {
    QoverPt_mid[i] = (QoverPt_min[i]+QoverPt_max[i])/2;
    QoverPt_interval[i] = (QoverPt_max[i] - QoverPt_min[i])/2;
  }
  const int nrun = 1;
  //int mbdrates[7] = {70, 250, 300, 380, 400, 430, 550};
  //int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
  //int mbdrates[nrun] = {250, 300, 380, 400, 430, 550};
  //int runs[nrun] = {53534, 53744, 53756, 53877, 53876, 53630};
  int mbdrates[nrun] = {400};
  int runs[nrun] = {53877};

  for (int k=0; k<nrun; k++)
  {
    TChain* chain = new TChain("residualtree");
    //chain->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/staticCorrOn_scale1_BenMap_SiTPCMatchingWindow/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runs[k],runs[k]));
    for (int i=0; i<100; i++)
    {
     chain->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/staticCorrOn_scale1_BenMap_SiTPCMatchingWindow/Reconstructed/%d/clusters_seeds_%d-%d.root_resid.root",runs[k],runs[k],i));
    }

    TCut cut = "crossing==0 && ntpc>34 && !std::isnan(silseedeta)";

    /*
    TH1* h1_QoverPt = new TH1F("h1_QoverPt","TPC seed Q/pT;Q_{TPC}/pT_{TPC};Counts",100,-8,8);
    plot1D(h1_QoverPt, chain, Form("m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpx)+TMath::Sq(m_tpcseedpy))"), cut, Form("figure/%d_QoverPt.pdf",runs[k]), 1);

    TH1* h1_dx = new TH1F("h1_dx","Inclusive TPCx - Six;x_{TPC}-x_{Si};Counts",100,-5,5);
    plot1D(h1_dx, chain, Form("m_tpcseedx-m_silseedx"), cut, Form("figure/%d_dx_inclusive.pdf",runs[k]), 1);
    TH1* h1_dy = new TH1F("h1_dy","Inclusive TPCy - Siy;y_{TPC}-y_{Si};Counts",100,-5,5);
    plot1D(h1_dy, chain, Form("m_tpcseedy-m_silseedy"), cut, Form("figure/%d_dy_inclusive.pdf",runs[k]), 1);
    TH1* h1_dz = new TH1F("h1_dz","Inclusive TPCz - Siz;z_{TPC}-z_{Si};Counts",100,-10,10);
    plot1D(h1_dz, chain, Form("m_tpcseedz-m_silseedz"), cut, Form("figure/%d_dz_inclusive.pdf",runs[k]), 1);
    TH1* h1_dz_etapos = new TH1F("h1_dz_etapos","Inclusive TPCz - Siz #eta>0;z_{TPC}-z_{Si};Counts",100,-10,10);
    plot1D(h1_dz_etapos, chain, Form("m_tpcseedz-m_silseedz"), cut+"tpcseedeta>0", Form("figure/%d_dz_inclusive_etapos.pdf",runs[k]), 1);
    TH1* h1_dz_etaneg = new TH1F("h1_dz_etaneg","Inclusive TPCz - Siz #eta<0;z_{TPC}-z_{Si};Counts",100,-10,10);
    plot1D(h1_dz_etaneg, chain, Form("m_tpcseedz-m_silseedz"), cut+"tpcseedeta<0", Form("figure/%d_dz_inclusive_etaneg.pdf",runs[k]), 1);
    TH1* h1_deta = new TH1F("h1_deta","Inclusive TPC#eta - Si#eta;#eta_{TPC}-#eta_{Si};Counts",100,-0.2,0.2);
    plot1D(h1_deta, chain, Form("m_tpcseedeta-m_silseedeta"), cut, Form("figure/%d_deta_inclusive.pdf",runs[k]), 1);
    TH1* h1_dphi = new TH1F("h1_dphi","Inclusive TPC#phi - Si#phi;#phi_{TPC}-#phi_{Si};Counts",100,-0.2,0.2);
    plot1D(h1_dphi, chain, Form("m_tpcseedphi-m_silseedphi"), cut, Form("figure/%d_dphi_inclusive.pdf",runs[k]), 1);

    TH2* h2_dx_QoverPt = new TH2F("h2_dx_QoverPt","Inclusive TPCx - Six vs Q_{TPC}/pT_{TPC};TPCx - Six;Q_{TPC}/pT_{TPC}",100,-5,5,100,-8,8);
    plot2D(h2_dx_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpx)+TMath::Sq(m_tpcseedpy)):m_tpcseedx-m_silseedx", cut, Form("figure/%d_dx_vs_QoverPt.pdf",runs[k]), 1);
    TH2* h2_dy_QoverPt = new TH2F("h2_dy_QoverPt","Inclusive TPCy - Siy vs Q_{TPC}/pT_{TPC};TPCy - Siy;Q_{TPC}/pT_{TPC}",100,-5,5,100,-8,8);
    plot2D(h2_dy_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpx)+TMath::Sq(m_tpcseedpy)):m_tpcseedy-m_silseedy", cut, Form("figure/%d_dy_vs_QoverPt.pdf",runs[k]), 1);
    TH2* h2_dz_QoverPt = new TH2F("h2_dz_QoverPt","Inclusive TPCz - Siz vs Q_{TPC}/pT_{TPC};TPCz - Siz;Q_{TPC}/pT_{TPC}",100,-10,10,100,-8,8);
    plot2D(h2_dz_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpz)+TMath::Sq(m_tpcseedpy)):m_tpcseedz-m_silseedz", cut, Form("figure/%d_dz_vs_QoverPt.pdf",runs[k]), 1);
    TH2* h2_dz_etapos_QoverPt = new TH2F("h2_dz_etapos_QoverPt","Inclusive TPCz - Siz vs Q_{TPC}/pT_{TPC} #eta>0;TPCz - Siz;Q_{TPC}/pT_{TPC}",100,-10,10,100,-8,8);
    plot2D(h2_dz_etapos_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpz)+TMath::Sq(m_tpcseedpy)):m_tpcseedz-m_silseedz", cut+"tpcseedeta>0", Form("figure/%d_dz_vs_QoverPt_etapos.pdf",runs[k]), 1);
    TH2* h2_dz_etaneg_QoverPt = new TH2F("h2_dz_etaneg_QoverPt","Inclusive TPCz - Siz vs Q_{TPC}/pT_{TPC} #eta<0;TPCz - Siz;Q_{TPC}/pT_{TPC}",100,-10,10,100,-8,8);
    plot2D(h2_dz_etaneg_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpz)+TMath::Sq(m_tpcseedpy)):m_tpcseedz-m_silseedz", cut+"tpcseedeta<0", Form("figure/%d_dz_vs_QoverPt_etaneg.pdf",runs[k]), 1);
    TH2* h2_deta_QoverPt = new TH2F("h2_deta_QoverPt","Inclusive TPC#eta - Si#eta vs Q_{TPC}/pT_{TPC};TPC#eta - Si#eta;Q_{TPC}/pT_{TPC}",100,-0.2,0.2,100,-8,8);
    plot2D(h2_deta_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpz)+TMath::Sq(m_tpcseedpy)):m_tpcseedeta-m_silseedeta", cut, Form("figure/%d_deta_vs_QoverPt.pdf",runs[k]), 1);
    TH2* h2_dphi_QoverPt = new TH2F("h2_dphi_QoverPt","Inclusive TPC#phi - Si#phi vs Q_{TPC}/pT_{TPC};TPC#phi - Si#phi;Q_{TPC}/pT_{TPC}",100,-0.2,0.2,100,-8,8);
    plot2D(h2_dphi_QoverPt, chain, "m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpz)+TMath::Sq(m_tpcseedpy)):m_tpcseedphi-m_silseedphi", cut, Form("figure/%d_dphi_vs_QoverPt.pdf",runs[k]), 1);
    */

    TH1* h1_dx_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dx_QoverPt[i] = new TH1F(Form("h1_dx_QoverPt_%d",i),"TPCx - Six;x_{TPC}-x_{Si};Counts",100,-5,5);
    }
    plot1D_QoverPt(h1_dx_QoverPt, chain, Form("m_tpcseedx-m_silseedx"), cut, Form("figure/%d_dx_allQoverPt",runs[k]), -5, 5, 4, 0);

    TH1* h1_dy_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dy_QoverPt[i] = new TH1F(Form("h1_dy_QoverPt_%d",i),"TPCy - Siy;y_{TPC}-y_{Si};Counts",100,-5,5);
    }
    plot1D_QoverPt(h1_dy_QoverPt, chain, Form("m_tpcseedy-m_silseedy"), cut, Form("figure/%d_dy_allQoverPt",runs[k]), -5, 5, 4, 0);

    TH1* h1_dz_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dz_QoverPt[i] = new TH1F(Form("h1_dz_QoverPt_%d",i),"TPCz - Siz;z_{TPC}-z_{Si};Counts",250,-10,10);
    }
    plot1D_QoverPt(h1_dz_QoverPt, chain, Form("m_tpcseedz-m_silseedz"), cut, Form("figure/%d_dz_allQoverPt",runs[k]), -10, 10, 4, 0, 1);

    /*
    TH1* h1_dz_etapos_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dz_etapos_QoverPt[i] = new TH1F(Form("h1_dz_etapos_QoverPt_%d",i),"TPCz - Siz #eta>0;z_{TPC}-z_{Si};Counts",250,-10,10);
    }
    plot1D_QoverPt(h1_dz_etapos_QoverPt, chain, Form("m_tpcseedz-m_silseedz"), cut+"tpcseedeta>0", Form("figure/%d_dz_allQoverPt_etapos",runs[k]), -10, 10, 4, 0);

    TH1* h1_dz_etaneg_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dz_etaneg_QoverPt[i] = new TH1F(Form("h1_dz_etaneg_QoverPt_%d",i),"TPCz - Siz #eta<0;z_{TPC}-z_{Si};Counts",250,-10,10);
    }
    plot1D_QoverPt(h1_dz_etaneg_QoverPt, chain, Form("m_tpcseedz-m_silseedz"), cut+"tpcseedeta<0", Form("figure/%d_dz_allQoverPt_etaneg",runs[k]), -10, 10, 4, 0);

    TH1* h1_deta_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_deta_QoverPt[i] = new TH1F(Form("h1_deta_QoverPt_%d",i),"TPC#eta - Si#eta;#eta_{TPC}-#eta_{Si};Counts",300,-0.3,0.3);
    }
    plot1D_QoverPt(h1_deta_QoverPt, chain, Form("m_tpcseedeta-m_silseedeta"), cut, Form("figure/%d_deta_allQoverPt",runs[k]), -0.3, 0.3, 4, 0, 1);

    TH1* h1_dphi_QoverPt[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      h1_dphi_QoverPt[i] = new TH1F(Form("h1_dphi_QoverPt_%d",i),"TPC#phi - Si#phi;#phi_{TPC}-#phi_{Si};Counts",100,-0.2,0.2);
    }
    plot1D_QoverPt(h1_dphi_QoverPt, chain, Form("m_tpcseedphi-m_silseedphi"), cut, Form("figure/%d_dphi_allQoverPt",runs[k]), -0.2, 0.2, 4, 0);
    */
  }
}

void plot1D(TH1* h, TChain* chain, TString name, TCut cut, TString output, bool logy=0)
{
    chain->Draw((name + " >>" + h->GetName()), cut);
    TCanvas *can = new TCanvas("can","",800,600);
    gPad->SetLogy(logy);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h->SetLineColor(kBlack);
    h->Draw("hist");
    TLegend *legend = new TLegend(0.45, 0.75, 0.72, 0.9);
    legend->AddEntry(h, "", "l");
    //legend->Draw("same");
    can->Update();
    can->SaveAs(output);
}

void plot1D_QoverPt(TH1* h[], TChain* chain, TString name, TCut cut, TString output, double fit_xmin, double fit_xmax, int nsigmawindow=3, bool logy=0, bool doexpfit=0)
{
    int ymax=0;
    Double_t arr_mean[nQoverPtbin], arr_mean_error[nQoverPtbin];
    Double_t arr_sigma[nQoverPtbin], arr_sigma_error[nQoverPtbin];
    Double_t arr_nsigma[nQoverPtbin];
    Double_t arr_nsigmawidnow[nQoverPtbin], arr_nsigmawindow_error[nQoverPtbin];
    TF1 *gausFit[nQoverPtbin];
    for (int i=0; i<nQoverPtbin; i++)
    {
      std::cout<<"Fit to Q_{TPC}/pT_{TPC} in ["<<QoverPt_min[i]<<" , "<<QoverPt_max[i]<<"]"<<std::endl;
      TCut cut_QoverPt = Form("m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpx)+TMath::Sq(m_tpcseedpy))>%f && m_tpcseedcharge/TMath::Sqrt(TMath::Sq(m_tpcseedpx)+TMath::Sq(m_tpcseedpy))<%f",QoverPt_min[i],QoverPt_max[i]);
      chain->Draw((name + " >>" + h[i]->GetName()), cut+cut_QoverPt);
      if (h[i]->GetMaximum() > ymax)
      {
        ymax = h[i]->GetMaximum();
      }

      gausFit[i] = new TF1(Form("gausFit_%d",i), "gaus", h[i]->GetXaxis()->GetXmin(), h[i]->GetXaxis()->GetXmax());
      h[i]->Fit(gausFit[i], "Q", "", fit_xmin, fit_xmax);
      gausFit[i]->SetLineColor(kRed);
      //h[i]->GetXaxis()->SetRangeUser(-0.05, 0.05);

      Double_t mean = gausFit[i]->GetParameter(1);
      Double_t mean_error = gausFit[i]->GetParError(1);
      Double_t sigma = gausFit[i]->GetParameter(2);
      Double_t sigma_error = gausFit[i]->GetParError(2);
      arr_mean[i] = mean;
      arr_mean_error[i] = mean_error;
      arr_sigma[i] = sigma;
      arr_nsigma[i] = nsigmawindow*sigma;
      arr_sigma_error[i] = sigma_error;
      if (mean>0)
      {
        arr_nsigmawidnow[i] = mean + nsigmawindow*sigma;
      }
      else if (mean<0)
      {
        arr_nsigmawidnow[i] = fabs(mean - nsigmawindow*sigma);
      }
      arr_nsigmawindow_error[i] = sqrt(pow(mean_error,2)+pow(nsigmawindow*sigma_error,2));
    }

    if (doexpfit)
    {
      TGraphErrors* gr_nsigmawindow = new TGraphErrors(nQoverPtbin, QoverPt_mid, arr_nsigmawidnow, QoverPt_interval, arr_nsigmawindow_error);
      gr_nsigmawindow->GetXaxis()->SetTitle("Q_{TPC}/pT_{TPC}");
      gr_nsigmawindow->GetYaxis()->SetTitle(Form("Max Gaussian #mu #pm %d#sigma",nsigmawindow));
      //gr_nsigmawindow->SetMaximum(1);

      TF1 *fexp_neg = new TF1("fexp_neg", "[0]+[1]*exp(-[2]*x)", QoverPt_mid[0], 0);
      fexp_neg->SetParameter(0, arr_nsigmawidnow[nQoverPtbin/2-1]);
      //fexp_neg->SetParameter(1, 0.005);
      //fexp_neg->SetParameter(2, 1);
      fexp_neg->SetParameter(1, 0.02);//for dz
      fexp_neg->SetParameter(2, 0.5);//for dz
      fexp_neg->SetParLimits(0, 0.0, 2.0);//for dz
      fexp_neg->SetParLimits(1, 0.0, 1.0);//for dz
      fexp_neg->SetParLimits(2, 0.0, 10.0);//for dz
      gr_nsigmawindow->Fit(fexp_neg, "Q", "", QoverPt_mid[0], 0);
      fexp_neg->SetLineColor(kRed);
      Double_t a0_fexp_neg = fexp_neg->GetParameter(0);
      Double_t a0_fexp_neg_err = fexp_neg->GetParError(0);
      Double_t a1_fexp_neg = fexp_neg->GetParameter(1);
      Double_t a1_fexp_neg_err = fexp_neg->GetParError(1);
      Double_t a2_fexp_neg = fexp_neg->GetParameter(2);
      Double_t a2_fexp_neg_err = fexp_neg->GetParError(2);

      TF1 *fexp_pos = new TF1("fexp_pos", "[0]+[1]*exp([2]*x)", 0, QoverPt_mid[nQoverPtbin-1]);
      fexp_pos->SetParameter(0, arr_nsigmawidnow[nQoverPtbin/2-1]);
      //fexp_pos->SetParameter(1, 0.005);
      //fexp_pos->SetParameter(2, 1);
      fexp_pos->SetParameter(1, 0.02);//for dz
      fexp_pos->SetParameter(2, 0.5);//for dz
      fexp_pos->SetParLimits(0, 0.0, 2.0);//for dz
      fexp_pos->SetParLimits(1, 0.0, 1.0);//for dz
      fexp_pos->SetParLimits(2, 0.0, 10.0);//for dz
      gr_nsigmawindow->Fit(fexp_pos, "Q", "", 0, QoverPt_mid[nQoverPtbin-1]);
      fexp_pos->SetLineColor(kBlue);
      Double_t a0_fexp_pos = fexp_pos->GetParameter(0);
      Double_t a0_fexp_pos_err = fexp_pos->GetParError(0);
      Double_t a1_fexp_pos = fexp_pos->GetParameter(1);
      Double_t a1_fexp_pos_err = fexp_pos->GetParError(1);
      Double_t a2_fexp_pos = fexp_pos->GetParameter(2);
      Double_t a2_fexp_pos_err = fexp_pos->GetParError(2);

      TCanvas* can_fitwindow = new TCanvas("can_fitwindow","",800,600);

      can_fitwindow->cd(1);
      gr_nsigmawindow->Draw("AP");
      fexp_neg->Draw("same");
      fexp_pos->Draw("same");

      myText(0.25,0.85,kBlack,Form("Negative: %.3f+%.4f#timese^{%.2fQ/pT}",a0_fexp_neg,a1_fexp_neg,a2_fexp_neg));
      myText(0.25,0.75,kBlack,Form("Positive: %.3f+%.4f#timese^{%.2fQ/pT}",a0_fexp_pos,a1_fexp_pos,a2_fexp_pos));

      can_fitwindow->Update();
      can_fitwindow->SaveAs(output+"_fitwindow.pdf");

      TCanvas* can_fitwindow_sub = new TCanvas("can_fitwindow_sub","",800,600);

      can_fitwindow_sub->cd(1);
      //gr_nsigmawindow->SetMaximum(0.08);//for deta
      //gr_nsigmawindow->SetMinimum(0.01);//for deta
      gr_nsigmawindow->SetMaximum(2.7);//for dz
      gr_nsigmawindow->SetMinimum(1.1);//for dz
      gr_nsigmawindow->GetXaxis()->SetRangeUser(-2.5, 2.5);

      gr_nsigmawindow->Draw("AP");
      fexp_neg->Draw("same");
      fexp_pos->Draw("same");

      myText(0.25,0.85,kBlack,Form("Negative: %.3f+%.4f#timese^{%.2fQ/pT}",a0_fexp_neg,a1_fexp_neg,a2_fexp_neg));
      myText(0.25,0.75,kBlack,Form("Positive: %.3f+%.4f#timese^{%.2fQ/pT}",a0_fexp_pos,a1_fexp_pos,a2_fexp_pos));

      can_fitwindow_sub->Update();
      can_fitwindow_sub->SaveAs(output+"_fitwindow_sub.pdf");
    }

    TGraphErrors* gr_mean = new TGraphErrors(nQoverPtbin, QoverPt_mid, arr_mean, QoverPt_interval, arr_mean_error);
    TGraphErrors* gr_sigma = new TGraphErrors(nQoverPtbin, QoverPt_mid, arr_sigma, QoverPt_interval, arr_sigma_error);
    TGraphErrors* gr_mean_1sigma = new TGraphErrors(nQoverPtbin, QoverPt_mid, arr_mean, QoverPt_interval, arr_sigma);
    TGraphErrors* gr_mean_nsigma = new TGraphErrors(nQoverPtbin, QoverPt_mid, arr_mean, QoverPt_interval, arr_nsigma);

    gr_mean->GetXaxis()->SetTitle("Q_{TPC}/pT_{TPC}");
    gr_mean->GetYaxis()->SetTitle("Gaussian #mu #pm #mu error");
    gr_sigma->GetXaxis()->SetTitle("Q_{TPC}/pT_{TPC}");
    gr_sigma->GetYaxis()->SetTitle("Gaussian #sigma #pm #sigma error");
    gr_mean_1sigma->GetXaxis()->SetTitle("Q_{TPC}/pT_{TPC}");
    gr_mean_1sigma->GetYaxis()->SetTitle("Gaussian #mu #pm 1#sigma");
    gr_mean_nsigma->GetXaxis()->SetTitle("Q_{TPC}/pT_{TPC}");
    gr_mean_nsigma->GetYaxis()->SetTitle(Form("Gaussian #mu #pm %d#sigma",nsigmawindow));

    TCanvas* can_fit = new TCanvas("can_fit","",800,600);
    can_fit->Divide(2,2);

    can_fit->cd(1);
    gr_mean->Draw("AP");

    can_fit->cd(2);
    gr_sigma->Draw("AP");

    can_fit->cd(3);
    gr_mean_1sigma->Draw("AP");

    can_fit->cd(4);
    gr_mean_nsigma->Draw("AP");

    can_fit->Update();
    can_fit->SaveAs(output+"_fit.pdf");

    for (int i=0; i<nQoverPtbin; i++)
    {
      TCanvas *c = new TCanvas("c","",800,600);
      gPad->SetLogy(logy);
      gPad->SetRightMargin(0.05);
      h[i]->Draw("hist,e");
      gausFit[i]->Draw("same");

      TLegend *leg = new TLegend(0.62, 0.80, 0.82, 0.90);
      TString legstr = Form("Q_{TPC}/pT_{TPC} #in [%.2f,%.2f]",QoverPt_min[i],QoverPt_max[i]);
      leg->AddEntry(h[i], legstr, "l");
      leg->SetTextSize(0.03);
      leg->Draw();

      myText(0.65,0.75,kBlack,Form("#mu=%.4f#pm%.4f",arr_mean[i],arr_mean_error[i]));
      myText(0.65,0.70,kBlack,Form("#sigma=%.4f#pm%.4f",arr_sigma[i],arr_sigma_error[i]));

      c->Update();
      c->SaveAs(output+Form("_%d.pdf",i));
    }

    TCanvas *can = new TCanvas("can","",800,600);
    gPad->SetLogy(logy);

    for (int i=0; i<nQoverPtbin; i++)
    {
      h[i]->SetLineColor(i+1);
      h[i]->SetMaximum(ymax);
      if (i==0)
      {
        h[i]->Draw("hist");
      }
      else
      {
        h[i]->Draw("hist,same");
      }
    }
    can->Update();
    can->SaveAs(output+".pdf");

    TCanvas *can_leg = new TCanvas("can_leg","",500,600);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);

    TLegend *legend = new TLegend(0,0.03,1,0.97);
    for (int i=0; i<nQoverPtbin; i++)
    {
      TString legstr = Form("Q_{TPC}/pT_{TPC} #in [%.2f,%.2f]: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",QoverPt_min[i],QoverPt_max[i],arr_mean[i],arr_mean_error[i],arr_sigma[i],arr_sigma_error[i]);
      legend->AddEntry(h[i], legstr, "l");
    }
    legend->SetTextSize(0.03);
    legend->Draw();
    can_leg->Update();
    can_leg->SaveAs(output+"_leg.pdf");
}

void plot2D(TH2* h, TChain* chain, TString name, TCut cut, TString output, bool logz=0)
{
    chain->Draw((name + " >>" + h->GetName()), cut, "colz");
    TCanvas *can = new TCanvas("can","",800,600);
    gPad->SetLogz(logz);
    //gPad->SetLeftMargin(0.15);
    //gPad->SetRightMargin(0.15);
    h->Draw("colz");
    can->Update();
    can->SaveAs(output);
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
