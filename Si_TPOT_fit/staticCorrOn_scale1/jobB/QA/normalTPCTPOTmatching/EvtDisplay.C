namespace
{

  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  template <class T>
  inline constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    else if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }
    else
    {
      return phi;
    }
  }

}

void EvtDisplay()
{
  gStyle->SetOptStat(0);

  int runno=53877;
  int segmentno=7;
  int eventno=2172;

  TFile* infile = new TFile(Form("~/workarea/TPCdistortion/Si_TPOT_fit/staticCorrOn_scale1/Reconstructed/%d/clusters_seeds_%d-%d.root_qa.root",runno,runno,segmentno));
  TTree* intree = (TTree*) infile->Get("h_QAG4SimulationDistortions_residTree");
  float drphi, dz, clusZErr, stateZErr, clusZ, stateZ, clusR, clusRPhiErr, stateRPhiErr, clusPhi;
  float stateR, statePhi;
  float tanAlpha, tanBeta;
  int event;
  std::vector<float> *sclusgx_mvtx=0;
  std::vector<float> *sclusgy_mvtx=0;
  std::vector<float> *sclusgz_mvtx=0;
  std::vector<float> *sclusgr_mvtx=0;
  std::vector<float> *sclusphi_mvtx=0;
  std::vector<float> *sclusgx_intt=0;
  std::vector<float> *sclusgy_intt=0;
  std::vector<float> *sclusgz_intt=0;
  std::vector<float> *sclusgr_intt=0;
  std::vector<float> *sclusphi_intt=0;
  std::vector<float> *sclusgx_tpot=0;
  std::vector<float> *sclusgy_tpot=0;
  std::vector<float> *sclusgz_tpot=0;
  std::vector<float> *sclusgr_tpot=0;
  std::vector<float> *sclusphi_tpot=0;
  std::vector<float> *stategx_mvtx=0;
  std::vector<float> *stategy_mvtx=0;
  std::vector<float> *stategz_mvtx=0;
  std::vector<float> *stategr_mvtx=0;
  std::vector<float> *statephi_mvtx=0;
  std::vector<float> *stategx_intt=0;
  std::vector<float> *stategy_intt=0;
  std::vector<float> *stategz_intt=0;
  std::vector<float> *stategr_intt=0;
  std::vector<float> *statephi_intt=0;
  std::vector<float> *stategx_tpot=0;
  std::vector<float> *stategy_tpot=0;
  std::vector<float> *stategz_tpot=0;
  std::vector<float> *stategr_tpot=0;
  std::vector<float> *statephi_tpot=0;
  intree->SetBranchAddress("dz",&dz);
  intree->SetBranchAddress("drphi",&drphi);
  intree->SetBranchAddress("clusZErr",&clusZErr);
  intree->SetBranchAddress("stateZErr",&stateZErr);
  intree->SetBranchAddress("clusZ",&clusZ);
  intree->SetBranchAddress("stateZ",&stateZ);
  intree->SetBranchAddress("clusR",&clusR);
  intree->SetBranchAddress("clusPhi",&clusPhi);
  intree->SetBranchAddress("stateR",&stateR);
  intree->SetBranchAddress("statePhi",&statePhi);
  intree->SetBranchAddress("clusRPhiErr",&clusRPhiErr);
  intree->SetBranchAddress("stateRPhiErr",&stateRPhiErr);
  intree->SetBranchAddress("tanAlpha",&tanAlpha);
  intree->SetBranchAddress("tanBeta",&tanBeta);
  intree->SetBranchAddress("event",&event);
  intree->SetBranchAddress("sclusgx_mvtx",&sclusgx_mvtx);
  intree->SetBranchAddress("sclusgy_mvtx",&sclusgy_mvtx);
  intree->SetBranchAddress("sclusgz_mvtx",&sclusgz_mvtx);
  intree->SetBranchAddress("sclusgr_mvtx",&sclusgr_mvtx);
  intree->SetBranchAddress("sclusphi_mvtx",&sclusphi_mvtx);
  intree->SetBranchAddress("sclusgx_intt",&sclusgx_intt);
  intree->SetBranchAddress("sclusgy_intt",&sclusgy_intt);
  intree->SetBranchAddress("sclusgz_intt",&sclusgz_intt);
  intree->SetBranchAddress("sclusgr_intt",&sclusgr_intt);
  intree->SetBranchAddress("sclusphi_intt",&sclusphi_intt);
  intree->SetBranchAddress("sclusgx_tpot",&sclusgx_tpot);
  intree->SetBranchAddress("sclusgy_tpot",&sclusgy_tpot);
  intree->SetBranchAddress("sclusgz_tpot",&sclusgz_tpot);
  intree->SetBranchAddress("sclusgr_tpot",&sclusgr_tpot);
  intree->SetBranchAddress("sclusphi_tpot",&sclusphi_tpot);
  intree->SetBranchAddress("stategx_mvtx",&stategx_mvtx);
  intree->SetBranchAddress("stategy_mvtx",&stategy_mvtx);
  intree->SetBranchAddress("stategz_mvtx",&stategz_mvtx);
  intree->SetBranchAddress("stategr_mvtx",&stategr_mvtx);
  intree->SetBranchAddress("statephi_mvtx",&statephi_mvtx);
  intree->SetBranchAddress("stategx_intt",&stategx_intt);
  intree->SetBranchAddress("stategy_intt",&stategy_intt);
  intree->SetBranchAddress("stategz_intt",&stategz_intt);
  intree->SetBranchAddress("stategr_intt",&stategr_intt);
  intree->SetBranchAddress("statephi_intt",&statephi_intt);
  intree->SetBranchAddress("stategx_tpot",&stategx_tpot);
  intree->SetBranchAddress("stategy_tpot",&stategy_tpot);
  intree->SetBranchAddress("stategz_tpot",&stategz_tpot);
  intree->SetBranchAddress("stategr_tpot",&stategr_tpot);
  intree->SetBranchAddress("statephi_tpot",&statephi_tpot);
  std::vector<float> vec_clus_R, vec_clus_RPhi, vec_clus_x, vec_clus_y, vec_clus_z;
  std::vector<float> vec_state_R, vec_state_RPhi, vec_state_x, vec_state_y, vec_state_z;
  std::vector<float> vec_clusR_statePhi, vec_state_x2, vec_state_y2;
  std::vector<float> vec_drphi, vec_dz;

  vec_clus_R.clear();
  vec_clus_RPhi.clear();
  vec_clus_x.clear();
  vec_clus_y.clear();
  vec_clus_z.clear();
  vec_state_R.clear();
  vec_state_RPhi.clear();
  vec_state_x.clear();
  vec_state_y.clear();
  vec_state_z.clear();
  vec_clusR_statePhi.clear();
  vec_state_x2.clear();
  vec_state_y2.clear();
  vec_drphi.clear();
  vec_dz.clear();

  float min_R = 1e3, min_RPhi = 1e3, min_x = 1e3, min_y = 1e3, min_z = 1e3;
  float max_R = -1e3, max_RPhi = -1e3, max_x = -1e3, max_y = -1e3, max_z = -1e3;
  float min_RPhi2 = 1e3, max_RPhi2 = -1e3, min_x2 = 1e3, min_y2 = 1e3, max_x2 = -1e3, max_y2 = -1e3;

  int nevent = intree->GetEntries();
  int flag = true;
  for (int i=0; i<nevent; i++)
  {
    intree->GetEntry(i);
    if (event!=eventno) continue;
    double erp = square(clusRPhiErr) + square(stateRPhiErr);
    double ez = square(clusZErr) + square(stateZErr);
    //if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
    //if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
    //if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;
  
    //if (std::fabs(clusZ) > 20) continue;

    vec_clus_R.push_back(clusR);
    vec_clus_RPhi.push_back(clusR*clusPhi);
    vec_clus_x.push_back(clusR*cos(clusPhi));
    vec_clus_y.push_back(clusR*sin(clusPhi));
    vec_clus_z.push_back(clusZ);
    vec_state_R.push_back(stateR);
    vec_state_RPhi.push_back(stateR*statePhi);
    vec_state_x.push_back(stateR*cos(statePhi));
    vec_state_y.push_back(stateR*sin(statePhi));
    vec_state_z.push_back(stateZ);
    vec_clusR_statePhi.push_back(clusR*statePhi);
    vec_state_x2.push_back(clusR*cos(statePhi));
    vec_state_y2.push_back(clusR*sin(statePhi));
    vec_drphi.push_back(drphi);
    vec_dz.push_back(dz);

    if (flag)
    {
      int nmvtx = sclusgx_mvtx->size();
      for (int imvtx = 0; imvtx < nmvtx; imvtx++)
      {
          vec_clus_R.push_back(sclusgr_mvtx->at(imvtx));
          vec_clus_RPhi.push_back(sclusgr_mvtx->at(imvtx)*sclusphi_mvtx->at(imvtx));
          vec_clus_x.push_back(sclusgx_mvtx->at(imvtx));
          vec_clus_y.push_back(sclusgy_mvtx->at(imvtx));
          vec_clus_z.push_back(sclusgz_mvtx->at(imvtx));
          vec_state_R.push_back(stategr_mvtx->at(imvtx));
          vec_state_RPhi.push_back(stategr_mvtx->at(imvtx)*statephi_mvtx->at(imvtx));
          vec_state_x.push_back(stategx_mvtx->at(imvtx));
          vec_state_y.push_back(stategy_mvtx->at(imvtx));
          vec_state_z.push_back(stategz_mvtx->at(imvtx));
          vec_clusR_statePhi.push_back(sclusgr_mvtx->at(imvtx)*statephi_mvtx->at(imvtx));
          vec_state_x2.push_back(sclusgr_mvtx->at(imvtx)*cos(statephi_mvtx->at(imvtx)));
          vec_state_y2.push_back(sclusgr_mvtx->at(imvtx)*sin(statephi_mvtx->at(imvtx)));
          vec_drphi.push_back(sclusgr_mvtx->at(imvtx)*(sclusphi_mvtx->at(imvtx)-statephi_mvtx->at(imvtx)));
          vec_dz.push_back(sclusgz_mvtx->at(imvtx)-stategz_mvtx->at(imvtx));
      }
      int nintt = sclusgx_intt->size();
      for (int iintt = 0; iintt < nintt; iintt++)
      {
          vec_clus_R.push_back(sclusgr_intt->at(iintt));
          vec_clus_RPhi.push_back(sclusgr_intt->at(iintt)*sclusphi_intt->at(iintt));
          vec_clus_x.push_back(sclusgx_intt->at(iintt));
          vec_clus_y.push_back(sclusgy_intt->at(iintt));
          vec_clus_z.push_back(sclusgz_intt->at(iintt));
          vec_state_R.push_back(stategr_intt->at(iintt));
          vec_state_RPhi.push_back(stategr_intt->at(iintt)*statephi_intt->at(iintt));
          vec_state_x.push_back(stategx_intt->at(iintt));
          vec_state_y.push_back(stategy_intt->at(iintt));
          vec_state_z.push_back(stategz_intt->at(iintt));
          vec_clusR_statePhi.push_back(sclusgr_intt->at(iintt)*statephi_intt->at(iintt));
          vec_state_x2.push_back(sclusgr_intt->at(iintt)*cos(statephi_intt->at(iintt)));
          vec_state_y2.push_back(sclusgr_intt->at(iintt)*sin(statephi_intt->at(iintt)));
          vec_drphi.push_back(sclusgr_intt->at(iintt)*(sclusphi_intt->at(iintt)-statephi_intt->at(iintt)));
          vec_dz.push_back(sclusgz_intt->at(iintt)-stategz_intt->at(iintt));
      }
      int ntpot = sclusgx_tpot->size();
      for (int itpot = 0; itpot < ntpot; itpot++)
      {
          vec_clus_R.push_back(sclusgr_tpot->at(itpot));
          vec_clus_RPhi.push_back(sclusgr_tpot->at(itpot)*sclusphi_tpot->at(itpot));
          vec_clus_x.push_back(sclusgx_tpot->at(itpot));
          vec_clus_y.push_back(sclusgy_tpot->at(itpot));
          vec_clus_z.push_back(sclusgz_tpot->at(itpot));
          vec_state_R.push_back(stategr_tpot->at(itpot));
          vec_state_RPhi.push_back(stategr_tpot->at(itpot)*statephi_tpot->at(itpot));
          vec_state_x.push_back(stategx_tpot->at(itpot));
          vec_state_y.push_back(stategy_tpot->at(itpot));
          vec_state_z.push_back(stategz_tpot->at(itpot));
          vec_clusR_statePhi.push_back(sclusgr_tpot->at(itpot)*statephi_tpot->at(itpot));
          vec_state_x2.push_back(sclusgr_tpot->at(itpot)*cos(statephi_tpot->at(itpot)));
          vec_state_y2.push_back(sclusgr_tpot->at(itpot)*sin(statephi_tpot->at(itpot)));
          vec_drphi.push_back(sclusgr_tpot->at(itpot)*(sclusphi_tpot->at(itpot)-statephi_tpot->at(itpot)));
          vec_dz.push_back(sclusgz_tpot->at(itpot)-stategz_tpot->at(itpot));
      }
      flag=false;
    }
  }

  int n = vec_clus_x.size();
  for (int i=0; i<n; i++)
  {
    if (vec_clus_R[i] < min_R) min_R = vec_clus_R[i];
    if (vec_state_R[i] < min_R) min_R = vec_state_R[i];
    if (vec_clus_R[i] > max_R) max_R = vec_clus_R[i];
    if (vec_state_R[i] > max_R) max_R = vec_state_R[i];

    if (vec_clus_RPhi[i] < min_RPhi) min_RPhi = vec_clus_RPhi[i];
    if (vec_state_RPhi[i] < min_RPhi) min_RPhi = vec_state_RPhi[i];
    if (vec_clus_RPhi[i] > max_RPhi) max_RPhi = vec_clus_RPhi[i];
    if (vec_state_RPhi[i] > max_RPhi) max_RPhi = vec_state_RPhi[i];

    if (vec_clus_x[i] < min_x) min_x = vec_clus_x[i];
    if (vec_state_x[i] < min_x) min_x = vec_state_x[i];
    if (vec_clus_x[i] > max_x) max_x = vec_clus_x[i];
    if (vec_state_x[i] > max_x) max_x = vec_state_x[i];

    if (vec_clus_y[i] < min_y) min_y = vec_clus_y[i];
    if (vec_state_y[i] < min_y) min_y = vec_state_y[i];
    if (vec_clus_y[i] > max_y) max_y = vec_clus_y[i];
    if (vec_state_y[i] > max_y) max_y = vec_state_y[i];

    if (vec_clus_z[i] < min_z) min_z = vec_clus_z[i];
    if (vec_state_z[i] < min_z) min_z = vec_state_z[i];
    if (vec_clus_z[i] > max_z) max_z = vec_clus_z[i];
    if (vec_state_z[i] > max_z) max_z = vec_state_z[i];

    if (vec_clus_RPhi[i] < min_RPhi2) min_RPhi2 = vec_clus_RPhi[i];
    if (vec_clusR_statePhi[i] < min_RPhi2) min_RPhi2 = vec_clusR_statePhi[i];
    if (vec_clus_RPhi[i] > max_RPhi2) max_RPhi2 = vec_clus_RPhi[i];
    if (vec_clusR_statePhi[i] > max_RPhi2) max_RPhi2 = vec_clusR_statePhi[i];

    if (vec_clus_x[i] < min_x2) min_x2 = vec_clus_x[i];
    if (vec_state_x2[i] < min_x2) min_x2 = vec_state_x2[i];
    if (vec_clus_x[i] > max_x2) max_x2 = vec_clus_x[i];
    if (vec_state_x2[i] > max_x2) max_x2 = vec_state_x2[i];

    if (vec_clus_y[i] < min_y2) min_y2 = vec_clus_y[i];
    if (vec_state_y2[i] < min_y2) min_y2 = vec_state_y2[i];
    if (vec_clus_y[i] > max_y2) max_y2 = vec_clus_y[i];
    if (vec_state_y2[i] > max_y2) max_y2 = vec_state_y2[i];
  }

  TCanvas *can = new TCanvas("can","",800*3,600*3);
  can->Divide(3,3);

  can->cd(1);
  can->SetTopMargin(0.12);
  TGraph *gr_clus_RPhi_R = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_clus_RPhi.data());
  gr_clus_RPhi_R->GetXaxis()->SetTitle("R (cm)");
  gr_clus_RPhi_R->GetYaxis()->SetTitle("R#phi (cm)");
  gr_clus_RPhi_R->GetXaxis()->SetLimits(min_R-2,max_R+2);
  gr_clus_RPhi_R->GetYaxis()->SetRangeUser(min_RPhi-2,max_RPhi+2);
  gr_clus_RPhi_R->SetMarkerColor(kBlack);
  gr_clus_RPhi_R->Draw("AP");
  TGraph *gr_state_RPhi_R = new TGraph(vec_state_R.size(), vec_state_R.data(), vec_state_RPhi.data());
  gr_state_RPhi_R->SetMarkerColor(kRed);
  gr_state_RPhi_R->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));
  TLegend *leg = new TLegend(0.6,0.7,0.8,0.9);
  leg->AddEntry(gr_clus_RPhi_R,"Cluster","p");
  leg->AddEntry(gr_state_RPhi_R,"State","p");
  leg->Draw("same");

  can->cd(2);
  can->SetTopMargin(0.12);
  TGraph *gr_clus_x_y = new TGraph(vec_clus_x.size(), vec_clus_x.data(), vec_clus_y.data());
  gr_clus_x_y->GetXaxis()->SetTitle("X (cm)");
  gr_clus_x_y->GetYaxis()->SetTitle("Y (cm)");
  gr_clus_x_y->GetXaxis()->SetLimits(min_x-2,max_x+2);
  gr_clus_x_y->GetYaxis()->SetRangeUser(min_y-2,max_y+2);
  gr_clus_x_y->SetMarkerColor(kBlack);
  gr_clus_x_y->Draw("AP");
  TGraph *gr_state_x_y = new TGraph(vec_state_x.size(), vec_state_x.data(), vec_state_y.data());
  gr_state_x_y->GetXaxis()->SetLimits(min_x-2,max_x+2);
  gr_state_x_y->GetYaxis()->SetRangeUser(min_y-2,max_y+2);
  gr_state_x_y->SetMarkerColor(kRed);
  gr_state_x_y->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));
  //leg->Draw("same");

  can->cd(3);
  can->SetTopMargin(0.12);
  TGraph *gr_clus_z_R = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_clus_z.data());
  gr_clus_z_R->GetXaxis()->SetTitle("R (cm)");
  gr_clus_z_R->GetYaxis()->SetTitle("Z (cm)");
  gr_clus_z_R->GetXaxis()->SetLimits(min_R-2,max_R+2);
  gr_clus_z_R->GetYaxis()->SetRangeUser(min_z-2,max_z+2);
  gr_clus_z_R->SetMarkerColor(kBlack);
  gr_clus_z_R->Draw("AP");
  TGraph *gr_state_z_R = new TGraph(vec_state_R.size(), vec_state_R.data(), vec_state_z.data());
  gr_state_z_R->SetMarkerColor(kRed);
  gr_state_z_R->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));
  //leg->Draw("same");

  can->cd(4);
  can->SetTopMargin(0.12);
  gr_clus_RPhi_R->Draw("AP");
  TGraph *gr_clusR_x_statePhi_R = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_clusR_statePhi.data());
  gr_clusR_x_statePhi_R->GetXaxis()->SetLimits(min_R-2,max_R+2);
  gr_clusR_x_statePhi_R->GetYaxis()->SetRangeUser(min_RPhi2-2,max_RPhi2+2);
  gr_clusR_x_statePhi_R->SetMarkerColor(kRed);
  gr_clusR_x_statePhi_R->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));
  leg->Draw("same");

  can->cd(5);
  can->SetTopMargin(0.12);
  gr_clus_x_y->Draw("AP");
  TGraph *gr_state_x_y2 = new TGraph(vec_state_x2.size(), vec_state_x2.data(), vec_state_y2.data());
  gr_state_x_y2->GetXaxis()->SetLimits(min_x-2,max_x+2);
  gr_state_x_y2->GetYaxis()->SetRangeUser(min_y-2,max_y+2);
  gr_state_x_y2->SetMarkerColor(kRed);
  gr_state_x_y2->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));
  //leg->Draw("same");

  can->cd(6);
  can->SetTopMargin(0.12);
  gr_clus_z_R->Draw("AP");
  TGraph *gr_state_z_clus_R = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_state_z.data());
  gr_state_z_clus_R->GetXaxis()->SetLimits(min_R-2,max_R+2);
  gr_state_z_clus_R->GetYaxis()->SetRangeUser(min_z-2,max_z+2);
  gr_state_z_clus_R->SetMarkerColor(kRed);
  gr_state_z_clus_R->Draw("P*,same");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));

  can->cd(7);
  can->SetTopMargin(0.12);
  TGraph *gr_clusR_drphi = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_drphi.data());
  gr_clusR_drphi->GetXaxis()->SetTitle("R (cm)");
  gr_clusR_drphi->GetYaxis()->SetTitle("Rd#phi (cm)");
  gr_clusR_drphi->SetMarkerColor(kBlack);
  gr_clusR_drphi->Draw("AP");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));

  can->cd(8);
  can->SetTopMargin(0.12);
  TGraph *gr_clusR_dz = new TGraph(vec_clus_R.size(), vec_clus_R.data(), vec_dz.data());
  gr_clusR_dz->GetXaxis()->SetTitle("R (cm)");
  gr_clusR_dz->GetYaxis()->SetTitle("dz (cm)");
  gr_clusR_dz->SetMarkerColor(kBlack);
  gr_clusR_dz->Draw("AP");
  myText(0.3,0.96,kBlack,Form("Run %d Segment %d Event %d",runno,segmentno,eventno));


  can->Update();
  can->SaveAs(Form("figure/EvtDisplay_run%d_segment%d_event%d.pdf",runno,segmentno,eventno));
  
}
