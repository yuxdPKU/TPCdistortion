void makecanvas1d(TH1* h1, TString name, TString xtitle, TString ytitle, bool logy=false);
void makecanvas2d(TH2* h2, TString name, TString xtitle, TString ytitle, bool logz=false);
void makecanvas2d_slicefit(TH2* h2, TString name, TString xtitle, TString ytitle, bool logz=false, bool dofit=false);
void makecanvas2d_slicefit_simul(TH2* h2_pos, TH2* h2_neg, TString name, TString xtitle, TString ytitle, bool logz=false, bool dofit=false);

double fit_function( double* x, double* par )
{
  const int ieta = x[0] >= 0 ? 1:0;
  const double z = x[1];
  const double offset = par[ieta];
  const double slope = par[2];
  return offset+slope*z;
}

double fit_function_1d( double* x, double* par )
{
  const double z = x[0];
  const double offset = par[0];
  const double slope = par[1];
  return offset+slope*z;
}

void si_tpc_z_resid_plot()
{
    TCut cut_nocrossingcut="ntpc>20&&nmaps>=3&&nintt>=2&&pt>-999&&quality<10000";
    TCut cut="ntpc>20&&nmaps>=3&&nintt>=2&&pt>-999&&quality<10000&&crossing==0";
    TCut cut_etapos="ntpc>20&&nmaps>=3&&nintt>=2&&pt>-999&&quality<10000&&eta>0&&crossing==0";
    TCut cut_etaneg="ntpc>20&&nmaps>=3&&nintt>=2&&pt>-999&&quality<10000&&eta<0&&crossing==0";

    //int runnumber=53877;
    int runnumber=52077;
    //TString type="CDBavgmap/root"; double dv=0.00726566;
    //TString type="CDBavgmap_dv743/root"; double dv=0.0074252;
    //TString type="CDBavgmap_tzero300/root"; double dv=0.00726566;
    //TString type="CDBavgmap_dv73/root"; double dv=0.0073;
    //TString type="CDBavgmap_dv75/root"; double dv=0.0074911;
    //TString type="CDBavgmap_dv75_t082/root"; double dv=0.0074911;
    //TString type="CDBavgmap_dv75_t082_sampa5/root"; double dv=0.0074911;
    //TString type="zeroField/root"; double dv=0.00745455;
    //TString type="zeroField_dv_round1/root"; double dv=0.00733896;
    //TString type="zeroField_dv_round2/root"; double dv=0.0073747274;
    //TString type="zeroField_dv_round3/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round4/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round5/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round6/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round7/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round8/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round9/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round10/root"; double dv=0.0073821095;
    //TString type="zeroField_dvt0_round11/root"; double dv=0.0073821095;
    TString type="zeroField_dvt0_round12/root"; double dv=0.0073821095;

    double deltaT=106.65237;

    TString outpath_prefix=Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/") + type + Form("/");

    TChain* chain = new TChain("residualtree");
    chain->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/") + type + Form("/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runnumber,runnumber));

    TH2* h2_dz_z = new TH2F("h2_dz_z","",40,-20,20,100,-4,4);
    chain->Draw("m_tpcseedz-m_silseedz:m_silseedz>>h2_dz_z",cut);
    makecanvas2d(h2_dz_z,outpath_prefix+Form("figure/dz_z.pdf"),Form("z_{Si}"),Form("z_{TPC} - z_{Si}"));

    TH2* h2_dz_z_etapos = new TH2F("h2_dz_z_etapos","",40,-15,15,100,-4,4);
    chain->Draw("m_tpcseedz-m_silseedz:m_silseedz>>h2_dz_z_etapos",cut_etapos);
    //makecanvas2d_slicefit(h2_dz_z_etapos,outpath_prefix+Form("figure/dz_z_etapos.pdf"),Form("z_{Si}"),Form("z_{TPC} - z_{Si}"),false,true);

    TH2* h2_dz_z_etaneg = new TH2F("h2_dz_z_etaneg","",40,-15,15,100,-4,4);
    chain->Draw("m_tpcseedz-m_silseedz:m_silseedz>>h2_dz_z_etaneg",cut_etaneg);
    //makecanvas2d_slicefit(h2_dz_z_etaneg,outpath_prefix+Form("figure/dz_z_etaneg.pdf"),Form("z_{Si}"),Form("z_{TPC} - z_{Si}"),false,true);

    makecanvas2d_slicefit_simul(h2_dz_z_etapos,h2_dz_z_etaneg,outpath_prefix+Form("figure/dz_z_simul.pdf"),Form("z_{Si}"),Form("z_{TPC} - z_{Si}"),false,true);

    TH2* h2_dz_eta = new TH2F("h2_dz_eta","",40,-1.5,1.5,100,-4,4);
    chain->Draw("m_tpcseedz-m_silseedz:m_eta>>h2_dz_eta",cut);
    makecanvas2d(h2_dz_eta,outpath_prefix+Form("figure/dz_eta.pdf"),Form("#eta"),Form("z_{TPC} - z_{Si}"));

    TH2* h2_dz_crossing = new TH2F("h2_dz_crossing","",500,-100,400,100,-4,4);
    TString ts_dz_crossing = Form("(m_tpcseedz-m_silseedz) + (eta>0?1:-1)*(m_crossing*%g*%g):m_crossing>>h2_dz_crossing",dv,deltaT);
    cout<<ts_dz_crossing<<endl;
    chain->Draw(ts_dz_crossing, cut_nocrossingcut);
    makecanvas2d(h2_dz_crossing,outpath_prefix+Form("figure/dz_crossing.pdf"),Form("crossing"),Form("z_{TPC} - z_{Si}"));


}

void makecanvas1d(TH1* h1, TString name, TString xtitle, TString ytitle="Counts", bool logy=false)
{
    TCanvas *can = new TCanvas("can","",800,600);
    can->cd();
    if (logy)
    {
      gPad->SetLogy(1);
      h1->SetMinimum(0.1);
    }
    else
    {
      gPad->SetLogy(0);
      h1->SetMinimum(0);
    }
    h1->GetXaxis()->SetTitle(xtitle);
    h1->GetYaxis()->SetTitle(ytitle);
    h1->Draw();
    can->SaveAs(name);
    delete can;
    gPad->SetLogy(0);
}

void makecanvas2d(TH2* h2, TString name, TString xtitle, TString ytitle, bool logz=false)
{
    TCanvas *can = new TCanvas("can","",800,600);
    can->cd();
    gPad->SetRightMargin(0.15);
    if (logz)
    {
      gPad->SetLogz(1);
      h2->SetMinimum(0.9);
    }
    h2->GetXaxis()->SetTitle(xtitle);
    h2->GetYaxis()->SetTitle(ytitle);
    h2->Draw("colz");

    can->SaveAs(name);
    delete can;
}

void makecanvas2d_slicefit(TH2* h2, TString name, TString xtitle, TString ytitle, bool logz=false, bool dofit=false)
{
    TCanvas *can = new TCanvas("can","",800,600);
    can->cd();
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.20);
    if (logz)
    {
      gPad->SetLogz(1);
      h2->SetMinimum(0.9);
    }
    h2->GetXaxis()->SetTitle(xtitle);
    h2->GetYaxis()->SetTitle(ytitle);
    h2->Draw("colz");

    /*
    Int_t nBinsX = h2->GetXaxis()->GetNbins();
    TH1D* h1_projY = new TH1D(Form("%s_projY",h2->GetName()), "Projection Y", nBinsX, h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());

    for (Int_t i = 1; i <= nBinsX; i++) {
        TH1D* h1_temp = h2->ProjectionY("_py", i, i);
        Double_t meanY = h1_temp->GetMean();
        Double_t meanErrorY = h1_temp->GetMeanError();
	//std::cout<<"sigma for bin "<<i<<" = "<<h1_temp->GetStdDev()<<std::endl;
        h1_projY->SetBinContent(i, meanY);
        h1_projY->SetBinError(i, meanErrorY);
        delete h1_temp;
    }
    h1_projY->Draw("E1,same");
    */

    h2->FitSlicesY();
    TH1 *wid = (TH1*)gDirectory->Get(Form("%s_2",h2->GetName()));
    TH1 *mean = (TH1*)gDirectory->Get(Form("%s_1",h2->GetName()));

    TF1 *fitfunc = new TF1("fitfunc", "[0]+[1]*x", -10, 10);
    if (dofit)
    {
      fitfunc->SetParameter(0,mean->GetBinContent(mean->GetNbinsX() / 2));
      fitfunc->SetParameter(1,0);
      cout<<"input intercept = "<<mean->GetBinContent(mean->GetNbinsX() / 2)<<endl;
      mean->Fit(fitfunc, "Q", "", -10, 10);
    }

    h2->Draw("colz");
    wid->SetLineColor(kRed);
    wid->SetMarkerColor(kRed);
    wid->Draw("same");
    mean->SetLineColor(kBlack);
    mean->Draw("same");
    fitfunc->SetLineColor(kViolet);
    if (dofit)
    {
      fitfunc->Draw("same");

      double intercept_mean = fitfunc->GetParameter(0);
      double intercept_error = fitfunc->GetParError(0);
      double slope_mean = fitfunc->GetParameter(1);
      double slope_error = fitfunc->GetParError(1);
      cout << "fit intercept = " << intercept_mean << " +/- " << intercept_error << endl;
      cout << "fit slope = " << slope_mean << " +/- " << slope_error << endl;

      TPaveText *pt = new TPaveText(.40, .00, .70, .15, "NDC");
      pt->SetFillColor(0);
      //pt->SetFillStyle(0);//transparent
      pt->SetLineColor(0);
      pt->SetBorderSize(0);
      pt->SetTextColor(kBlack);
      pt->AddText(Form("Slope: (%.2f#pm%.2f)#times10^{-2}",100*slope_mean,100*slope_error));
      pt->AddText(Form("Intercept: (%.2f#pm%.2f)#times10^{-2}",100*intercept_mean,100*intercept_error));
      pt->Draw("same");
    }

    for (int i = 1; i <= mean->GetNbinsX(); i++)
    {
    Double_t bin_center = mean->GetBinCenter(i);
    Double_t bin_content = mean->GetBinContent(i);
    Double_t bin_error = mean->GetBinError(i);
    Double_t low_edge = mean->GetBinLowEdge(i);
    Double_t up_edge = mean->GetBinLowEdge(i) + mean->GetBinWidth(i);

    //cout << "bin " << i << " center " << bin_center << " , content " << bin_content << endl;
    }

    can->SaveAs(name);
    delete can;
}

//refer to /direct/phenix+u/workarea/hpereira/sphenix/work/g4simulations/macros/SiliconMatchingDeltaZ_fit.C
void makecanvas2d_slicefit_simul(TH2* h2_pos, TH2* h2_neg, TString name, TString xtitle, TString ytitle, bool logz=false, bool dofit=false)
{
    TCanvas *can = new TCanvas("can","",1600,600);
    can->Divide(2,1);

    h2_pos->FitSlicesY();
    TH1 *wid_pos = (TH1*)gDirectory->Get(Form("%s_2",h2_pos->GetName()));
    TH1 *mean_pos = (TH1*)gDirectory->Get(Form("%s_1",h2_pos->GetName()));

    h2_neg->FitSlicesY();
    TH1 *wid_neg = (TH1*)gDirectory->Get(Form("%s_2",h2_neg->GetName()));
    TH1 *mean_neg = (TH1*)gDirectory->Get(Form("%s_1",h2_neg->GetName()));

    double zmin = h2_neg->GetXaxis()->GetXmin();
    double zmax = h2_neg->GetXaxis()->GetXmax();

    // create 2D histogram with mean and sigma from slice fit
    auto h_fit = new TH2F( "h_fit", "", 2, -1.5, 1.5, h2_neg->GetNbinsX(), zmin, zmax);
    h_fit->GetXaxis()->SetTitle( "#eta (tpc seed)" );
    h_fit->GetYaxis()->SetTitle( "z (silicon seed) (cm)" );
    h_fit->GetZaxis()->SetTitle( "#Deltaz(TPC-silicon) (cm)" );
    h_fit->SetStats(0);

    for( int i = 0; i<mean_neg->GetNbinsX(); ++i )
    {
      int ieta=0;
      const auto entries = h2_neg->Integral(i+1,i+1,1,h2_neg->GetNbinsY());
      if( entries > 0 )
      {
        h_fit->SetBinContent(ieta+1, i+1, mean_neg->GetBinContent(i+1) );
        h_fit->SetBinError(ieta+1, i+1, wid_neg->GetBinContent(i+1)/std::sqrt(entries));
      }
    }

    for( int i = 0; i<mean_pos->GetNbinsX(); ++i )
    {
      int ieta=1;
      const auto entries = h2_pos->Integral(i+1,i+1,1,h2_pos->GetNbinsY());
      if( entries > 0 )
      {
        h_fit->SetBinContent(ieta+1, i+1, mean_pos->GetBinContent(i+1) );
        h_fit->SetBinError(ieta+1, i+1, wid_pos->GetBinContent(i+1)/std::sqrt(entries));
      }
    }

    TF2 *fitfunc = new TF2("fitfunc", fit_function, -1.5, 1.5, zmin, zmax, 3); // 3 parameters
    if (dofit)
    {
      fitfunc->SetParameter(0, mean_pos->GetBinContent(mean_pos->GetNbinsX() / 2)); // hist1 截距初始值
      fitfunc->SetParameter(1, mean_neg->GetBinContent(mean_neg->GetNbinsX() / 2)); // hist2 截距初始值
      fitfunc->SetParameter(2, 0); // 设定初始斜率
      cout<<"pos input intercept = "<<mean_pos->GetBinContent(mean_pos->GetNbinsX() / 2)<<endl;
      cout<<"neg input intercept = "<<mean_neg->GetBinContent(mean_neg->GetNbinsX() / 2)<<endl;
      h_fit->Fit( fitfunc, "0R" );
    }

    can->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.20);
    if (logz)
    {
      gPad->SetLogz(1);
      h2_pos->SetMinimum(0.9);
    }
    h2_pos->GetXaxis()->SetTitle(xtitle);
    h2_pos->GetYaxis()->SetTitle(ytitle);
    h2_pos->Draw("colz");
    if (dofit)
    {
      int ieta=1;

      // project 2D histogram
      auto h_fit_proj = h_fit->ProjectionY( Form( "%s%i", h_fit->GetName(), ieta ), ieta+1, ieta+1);
      h_fit_proj->SetLineColor(kBlack);
      h_fit_proj->Draw("same");

      auto f_loc = new TF1( "f0", fit_function_1d, zmin, zmax, 2 );
      f_loc->SetParameter(0, fitfunc->GetParameter(ieta));
      f_loc->SetParameter(1, fitfunc->GetParameter(2));
      f_loc->SetLineColor(kViolet);
      f_loc->SetLineWidth(2);
      f_loc->Draw("same");

      double intercept_mean = fitfunc->GetParameter(ieta);
      double intercept_error = fitfunc->GetParError(ieta);
      double slope_mean = fitfunc->GetParameter(2);
      double slope_error = fitfunc->GetParError(2);

      TPaveText *pt = new TPaveText(.40, .00, .70, .15, "NDC");
      pt->SetFillColor(0);
      //pt->SetFillStyle(0);//transparent
      pt->SetLineColor(0);
      pt->SetBorderSize(0);
      pt->SetTextColor(kBlack);
      pt->AddText("#eta>0");
      pt->AddText(Form("Slope: (%.2f#pm%.2f)#times10^{-2}",100*slope_mean,100*slope_error));
      pt->AddText(Form("Intercept: (%.2f#pm%.2f)#times10^{-2}",100*intercept_mean,100*intercept_error));
      pt->Draw("same");
    }

    can->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.20);
    if (logz)
    {
      gPad->SetLogz(1);
      h2_neg->SetMinimum(0.9);
    }
    h2_neg->GetXaxis()->SetTitle(xtitle);
    h2_neg->GetYaxis()->SetTitle(ytitle);
    h2_neg->Draw("colz");
    if (dofit)
    {
      int ieta=0;

      // project 2D histogram
      auto h_fit_proj = h_fit->ProjectionY( Form( "%s%i", h_fit->GetName(), ieta ), ieta+1, ieta+1);
      h_fit_proj->SetLineColor(kBlack);
      h_fit_proj->Draw("same");

      auto f_loc = new TF1( "f0", fit_function_1d, zmin, zmax, 2 );
      f_loc->SetParameter(0, fitfunc->GetParameter(ieta));
      f_loc->SetParameter(1, fitfunc->GetParameter(2));
      f_loc->SetLineColor(kViolet);
      f_loc->SetLineWidth(2);
      f_loc->Draw("same");

      double intercept_mean = fitfunc->GetParameter(ieta);
      double intercept_error = fitfunc->GetParError(ieta);
      double slope_mean = fitfunc->GetParameter(2);
      double slope_error = fitfunc->GetParError(2);
      cout << "fit intercept = " << intercept_mean << " +/- " << intercept_error << endl;
      cout << "fit slope = " << slope_mean << " +/- " << slope_error << endl;

      TPaveText *pt = new TPaveText(.40, .00, .70, .15, "NDC");
      pt->SetFillColor(0);
      //pt->SetFillStyle(0);//transparent
      pt->SetLineColor(0);
      pt->SetBorderSize(0);
      pt->SetTextColor(kBlack);
      pt->AddText("#eta<0");
      pt->AddText(Form("Slope: (%.2f#pm%.2f)#times10^{-2}",100*slope_mean,100*slope_error));
      pt->AddText(Form("Intercept: (%.2f#pm%.2f)#times10^{-2}",100*intercept_mean,100*intercept_error));
      pt->Draw("same");
    }

    can->SaveAs(name);
    delete can;
}
