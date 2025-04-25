void makecanvas1d(TH1* h1, TString name, TString xtitle, TString ytitle, bool logy=false);
void makecanvas2d(TH2* h2, TString name, TString xtitle, TString ytitle, bool logz=false);
void residual_plot()
{
    TCut cut="ntpc>20 && nmaps>=2 && nintt>=1 && pt>=0.5";
    //TCut cut="crossing==0 && ntpc>20 && nmaps>=2 && nintt>=1 && pt>=0.5";
    //TCut cut="crossing==0 && ntpc>20 && nmaps>=2 && nintt>=1 && pt>=0.5 && fabs(silseedz-tpcseedz)<30";
    TCut cut_p=cut + "charge>0";
    TCut cut_m=cut + "charge<0";

    int runnumber=53534;
    //TString type="staticCorrOn_scale1_XudongMap";
    //TString type="staticCorrOn_scale1_BenMap";
    //TString type="staticCorrOn_scale1_noavgCorr";
    //TString type="staticCorrOn_scale1_BenMap_Misalignment1";
    //TString type="staticCorrOn_scale1_BenMap_Misalignment100";
    TString type="staticCorrOn_scale1_BenMap_Misalignment1_blowupTPC";
    //TString type="staticCorrOn_scale1_BenMap_newTPOTalignment";

    TString outpath_prefix=Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/") + type + Form("/");

    TChain* chain = new TChain("residualtree");
    chain->Add(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Full_fit/") + type + Form("/Reconstructed/%d/clusters_seeds_%d-*.root_resid.root",runnumber,runnumber));

/*
    TH2* h2_dphi_phi = new TH2F("h2_dphi_phi","",40,-0.2,0.2,100,-3.15,3.15);
    chain->Draw("m_tpcseedphi:m_tpcseedphi-m_silseedphi>>h2_dphi_phi",cut);
    makecanvas2d(h2_dphi_phi,outpath_prefix+Form("figure/dphi_phi.pdf"),Form("#phi_{TPC} - #phi_{Si}"),Form("#phi_{TPC}"));
*/

    /*
    TH1* h1_dcaxy = new TH1F("h1_dcaxy","",100,-5,5);
    chain->Draw("m_dcaxy>>h1_dcaxy",cut);
    makecanvas1d(h1_dcaxy,outpath_prefix+Form("figure/dcaxy.pdf"),Form("DCA_{xy} [cm]"),Form("Counts"),true);
    */

    TH1* h1_tpcseedx = new TH1F("h1_tpcseedx","",100,-10,10);
    chain->Draw("m_tpcseedx>>h1_tpcseedx",cut);
    makecanvas1d(h1_tpcseedx,outpath_prefix+Form("figure/tpcseedx.pdf"),Form("TPCseed x [cm]"),Form("Counts"),false);

    /*
    TH1* h1_dcaxy_tpcseed = new TH1F("h1_dcaxy_tpcseed","",100,0,30);
    chain->Draw("sqrt(m_tpcseedx*m_tpcseedx+m_tpcseedy*m_tpcseedy)>>h1_dcaxy_tpcseed",cut);
    makecanvas1d(h1_dcaxy_tpcseed,outpath_prefix+Form("figure/dcaxy_tpcseed.pdf"),Form("TPCseed DCA_{xy} [cm]"),Form("Counts"),true);

    TH1* h1_dcaxy_tpcseed_p = new TH1F("h1_dcaxy_tpcseed_p","",100,0,30);
    chain->Draw("sqrt(m_tpcseedx*m_tpcseedx+m_tpcseedy*m_tpcseedy)>>h1_dcaxy_tpcseed_p",cut_p);
    makecanvas1d(h1_dcaxy_tpcseed_p,outpath_prefix+Form("figure/dcaxy_tpcseed_p.pdf"),Form("TPCseed DCA_{xy} [cm]"),Form("Counts"),true);

    TH1* h1_dcaxy_tpcseed_m = new TH1F("h1_dcaxy_tpcseed_m","",100,0,30);
    chain->Draw("sqrt(m_tpcseedx*m_tpcseedx+m_tpcseedy*m_tpcseedy)>>h1_dcaxy_tpcseed_m",cut_m);
    makecanvas1d(h1_dcaxy_tpcseed_m,outpath_prefix+Form("figure/dcaxy_tpcseed_m.pdf"),Form("TPCseed DCA_{xy} [cm]"),Form("Counts"),true);
    */

/*
    TH1* h1_clus_dx = new TH1F("h1_clus_dx","",100,-1,1);
    chain->Draw("statelx-cluslx>>h1_clus_dx",cut);
    makecanvas1d(h1_clus_dx,outpath_prefix+Form("figure/clus_dx.pdf"),Form("statelx - cluslx [cm]"),Form("Counts"));

    TH1* h1_clus_dz = new TH1F("h1_clus_dz","",100,-1,1);
    chain->Draw("statelz-cluslz>>h1_clus_dz",cut);
    makecanvas1d(h1_clus_dz,outpath_prefix+Form("figure/clus_dz.pdf"),Form("statelz - cluslz [cm]"),Form("Counts"));
*/

    // cluster local dx versus layer
    TH2* h2_clus_dx_layer = new TH2F("h2_clus_dx_layer","",60,0,60,500,-0.5,0.5);
    //TH2* h2_clus_dx_layer = new TH2F("h2_clus_dx_layer","",60,0,60,500,-0.01,0.01); //for mvtx
    chain->Draw("(statelx-cluslx):cluslayer>>h2_clus_dx_layer",cut);
    makecanvas2d(h2_clus_dx_layer,outpath_prefix+Form("figure/clus_dx_layer.pdf"),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);

    TH2* h2_clus_dx_layer_p = new TH2F("h2_clus_dx_layer_p","",60,0,60,100,-1.5,1.5);
    chain->Draw("(statelx-cluslx):cluslayer>>h2_clus_dx_layer_p",cut_p);
    makecanvas2d(h2_clus_dx_layer_p,outpath_prefix+Form("figure/clus_dx_layer_p.pdf"),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);

    TH2* h2_clus_dx_layer_m = new TH2F("h2_clus_dx_layer_m","",60,0,60,100,-1.5,1.5);
    chain->Draw("(statelx-cluslx):cluslayer>>h2_clus_dx_layer_m",cut_m);
    makecanvas2d(h2_clus_dx_layer_m,outpath_prefix+Form("figure/clus_dx_layer_m.pdf"),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);

    float Zmin=-100;
    float Zmax=100;
    int nz=10;
    for (int i=0; i<nz; i++)
    {
      float zmin=Zmin+i*(Zmax-Zmin)/nz;
      float zmax=Zmin+(i+1)*(Zmax-Zmin)/nz;
      TH2* h2_clus_dx_layer_z = new TH2F(Form("h2_clus_dx_layer_z%d",i),Form("%.0f<Z<%.f cm",zmin,zmax),60,0,60,100,-1.5,1.5);
      chain->Draw(Form("(statelx-cluslx):cluslayer>>h2_clus_dx_layer_z%d",i),cut + Form("clusgz>%f && clusgz<%f",zmin,zmax));
      makecanvas2d(h2_clus_dx_layer_z,outpath_prefix+Form("figure/clus_dx_layer_z%d.pdf",i),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);

      /*
      TH2* h2_clus_dx_layer_p_z = new TH2F(Form("h2_clus_dx_layer_p_z%d",i),Form("%.0f<Z<%.f cm",zmin,zmax),60,0,60,100,-1.5,1.5);
      chain->Draw(Form("(statelx-cluslx):cluslayer>>h2_clus_dx_layer_p_z%d",i),cut_p + Form("clusgz>%f && clusgz<%f",zmin,zmax));
      makecanvas2d(h2_clus_dx_layer_p_z,outpath_prefix+Form("figure/clus_dx_layer_p_z%d.pdf",i),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);

      TH2* h2_clus_dx_layer_m_z = new TH2F(Form("h2_clus_dx_layer_m_z%d",i),Form("%.0f<Z<%.f cm",zmin,zmax),60,0,60,100,-1.5,1.5);
      chain->Draw(Form("(statelx-cluslx):cluslayer>>h2_clus_dx_layer_m_z%d",i),cut_m + Form("clusgz>%f && clusgz<%f",zmin,zmax));
      makecanvas2d(h2_clus_dx_layer_m_z,outpath_prefix+Form("figure/clus_dx_layer_m_z%d.pdf",i),Form("Layer"),Form("(statelx-cluslx) [cm]"),true);
      */
    }

    // cluster local dz versus layer
    TH2* h2_clus_dz_layer = new TH2F("h2_clus_dz_layer","",60,0,60,100,-2,2);
    chain->Draw("(statelz-cluslz):cluslayer>>h2_clus_dz_layer",cut);
    makecanvas2d(h2_clus_dz_layer,outpath_prefix+Form("figure/clus_dz_layer.pdf"),Form("Layer"),Form("(statelz-cluslz) [cm]"),true);

    TH2* h2_clus_dz_layer_p = new TH2F("h2_clus_dz_layer_p","",60,0,60,100,-2,2);
    chain->Draw("(statelz-cluslz):cluslayer>>h2_clus_dz_layer_p",cut_p);
    makecanvas2d(h2_clus_dz_layer_p,outpath_prefix+Form("figure/clus_dz_layer_p.pdf"),Form("Layer"),Form("(statelz-cluslz) [cm]"),true);

    TH2* h2_clus_dz_layer_m = new TH2F("h2_clus_dz_layer_m","",60,0,60,100,-2,2);
    chain->Draw("(statelz-cluslz):cluslayer>>h2_clus_dz_layer_m",cut_m);
    makecanvas2d(h2_clus_dz_layer_m,outpath_prefix+Form("figure/clus_dz_layer_m.pdf"),Form("Layer"),Form("(statelz-cluslz) [cm]"),true);

    /*
    TH2* h2_clus_resid_layer = new TH2F("h2_clus_resid_layer","",60,0,60,100,0,1);
    chain->Draw("sqrt((statelx-cluslx)*(statelx-cluslx)+(statelz-cluslz)*(statelz-cluslz)):cluslayer>>h2_clus_resid_layer",cut);
    makecanvas2d(h2_clus_resid_layer,outpath_prefix+Form("figure/clus_resid_layer.pdf"),Form("Layer"),Form("#sqrt((statelx-cluslx)^2+(statelz-cluslz)^2) [cm]"));
    */

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
    wid->SetLineColor(kRed);
    wid->Draw("same");
    mean->SetLineColor(kBlack);
    mean->Draw("same");

    can->SaveAs(name);
    delete can;
}
