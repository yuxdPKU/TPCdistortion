void fit_gauss(TH2* h, TString name, bool verbose=0);

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

float getPhi(float phi)
{
  float m_phiMin = 0;
  float m_phiMax = 2. * M_PI;
  while (phi < m_phiMin)
  {
    phi += 2. * M_PI;
  }
  while (phi >= m_phiMax)
  {
    phi -= 2. * M_PI;
  }
  return phi;
}

void fit_2d_residual()
{
  int verbosity = 0;
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
    TFile* infile = new TFile(Form("allqa_%d.root",runs[k]));
    TTree* intree = (TTree*) infile->Get("h_QAG4SimulationDistortions_residTree");
    float drphi, dz, clusZErr, stateZErr, clusZ, clusR, clusRPhiErr, stateRPhiErr, clusPhi;
    float tanAlpha, tanBeta;
    float trackPt;
    int layer, charge;
    int track_nmvtx, track_nintt, track_ntpc, track_ntpot;
    int track_nmvtxstate, track_ninttstate, track_ntpcstate, track_ntpotstate;
    intree->SetBranchAddress("dz",&dz);
    intree->SetBranchAddress("drphi",&drphi);
    intree->SetBranchAddress("clusZErr",&clusZErr);
    intree->SetBranchAddress("stateZErr",&stateZErr);
    intree->SetBranchAddress("clusZ",&clusZ);
    intree->SetBranchAddress("clusR",&clusR);
    intree->SetBranchAddress("clusPhi",&clusPhi);
    intree->SetBranchAddress("clusRPhiErr",&clusRPhiErr);
    intree->SetBranchAddress("stateRPhiErr",&stateRPhiErr);
    intree->SetBranchAddress("tanAlpha",&tanAlpha);
    intree->SetBranchAddress("trackPt",&trackPt);
    intree->SetBranchAddress("tanBeta",&tanBeta);
    intree->SetBranchAddress("layer",&layer);
    intree->SetBranchAddress("charge",&charge);
    intree->SetBranchAddress("track_nmvtx",&track_nmvtx);
    intree->SetBranchAddress("track_nintt",&track_nintt);
    intree->SetBranchAddress("track_ntpc",&track_ntpc);
    intree->SetBranchAddress("track_ntpot",&track_ntpot);
    intree->SetBranchAddress("track_nmvtxstate",&track_nmvtxstate);
    intree->SetBranchAddress("track_ninttstate",&track_ninttstate);
    intree->SetBranchAddress("track_ntpcstate",&track_ntpcstate);
    intree->SetBranchAddress("track_ntpotstate",&track_ntpotstate);
  
    const int nRBins = 16;
    const int nZBins = 40;

    const double rMin = 20;
    const double rMax = 78;
    const double zMin = 0;
    const double zMax = 105.5;
    double rWidth = (rMax-rMin) / nRBins;
    double zWidth = (zMax-zMin) / nZBins;

    TGraphErrors* gr_drphi_tanAlpha_posz[nRBins][nZBins];
    TGraphErrors* gr_drphi_tanAlpha_negz[nRBins][nZBins];
    TGraphErrors* gr_dz_tanBeta_posz[nRBins][nZBins];
    TGraphErrors* gr_dz_tanBeta_negz[nRBins][nZBins];

    std::vector<double> vec_drphi_posz[nRBins][nZBins];
    std::vector<double> vec_drphi_negz[nRBins][nZBins];
    std::vector<double> vec_dz_posz[nRBins][nZBins];
    std::vector<double> vec_dz_negz[nRBins][nZBins];
 
    std::vector<double> vec_drphi_err_posz[nRBins][nZBins];
    std::vector<double> vec_drphi_err_negz[nRBins][nZBins];
    std::vector<double> vec_dz_err_posz[nRBins][nZBins];
    std::vector<double> vec_dz_err_negz[nRBins][nZBins];
  
    std::vector<double> vec_tanAlpha_posz[nRBins][nZBins];
    std::vector<double> vec_tanAlpha_negz[nRBins][nZBins];
    std::vector<double> vec_tanBeta_posz[nRBins][nZBins];
    std::vector<double> vec_tanBeta_negz[nRBins][nZBins];

    std::vector<double> vec_tanAlpha_err_posz[nRBins][nZBins];
    std::vector<double> vec_tanAlpha_err_negz[nRBins][nZBins];
    std::vector<double> vec_tanBeta_err_posz[nRBins][nZBins];
    std::vector<double> vec_tanBeta_err_negz[nRBins][nZBins];

    for (int i = 0; i < nRBins; i++) {
        for (int j = 0; j < nZBins; j++) {
            vec_drphi_posz[i][j].clear();
            vec_drphi_negz[i][j].clear();
            vec_dz_posz[i][j].clear();
            vec_dz_negz[i][j].clear();
            vec_drphi_err_posz[i][j].clear();
            vec_drphi_err_negz[i][j].clear();
            vec_dz_err_posz[i][j].clear();
            vec_dz_err_negz[i][j].clear();
            vec_tanAlpha_posz[i][j].clear();
            vec_tanAlpha_negz[i][j].clear();
            vec_tanBeta_posz[i][j].clear();
            vec_tanBeta_negz[i][j].clear();
            vec_tanAlpha_err_posz[i][j].clear();
            vec_tanAlpha_err_negz[i][j].clear();
            vec_tanBeta_err_posz[i][j].clear();
            vec_tanBeta_err_negz[i][j].clear();
        }
    }
  
    int nevent = intree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
      intree->GetEntry(i);
      double erp = square(clusRPhiErr) + square(stateRPhiErr);
      double ez = square(clusZErr) + square(stateZErr);
      if (sqrt(erp) < 0.005) continue; if (sqrt(ez) < 0.01) continue;
      if (std::abs(tanAlpha) > 0.6 || std::abs(drphi) > 2) continue;
      if (std::abs(tanBeta) > 1.5 || std::abs(dz) > 5) continue;
      if (track_nmvtxstate<3) continue;
      if (track_ninttstate<2) continue;
      if (track_ntpc<35) continue;
      if (track_ntpotstate<2) continue;

      int index_r = -1;
      for (int ir = 0; ir < nRBins; ir++)
      {
        if (clusR > (rMin + ir*rWidth) && clusR < (rMin + (ir+1)*rWidth))
        {
          index_r = ir;
        }
      }
      int index_z = -1;
      for (int iz = 0; iz < nZBins; iz++)
      {
        if (fabs(clusZ) > (zMin + iz*zWidth) && fabs(clusZ) < (zMin + (iz+1)*zWidth))
        {
          index_z = iz;
        }
      }
      if (index_r == -1 || index_z == -1)
      {
        continue;
      }

      if (clusZ>0)
      {
        vec_drphi_posz[index_r][index_z].push_back(drphi);
        vec_dz_posz[index_r][index_z].push_back(dz);
        vec_tanAlpha_posz[index_r][index_z].push_back(tanAlpha);
        vec_tanBeta_posz[index_r][index_z].push_back(tanBeta);
        vec_drphi_err_posz[index_r][index_z].push_back(sqrt(erp));
        vec_dz_err_posz[index_r][index_z].push_back(sqrt(ez));
        vec_tanAlpha_err_posz[index_r][index_z].push_back(0);
        vec_tanBeta_err_posz[index_r][index_z].push_back(0);
      }
      else if (clusZ<0)
      {
        vec_drphi_negz[index_r][index_z].push_back(drphi);
        vec_dz_negz[index_r][index_z].push_back(dz);
        vec_tanAlpha_negz[index_r][index_z].push_back(tanAlpha);
        vec_tanBeta_negz[index_r][index_z].push_back(tanBeta);
        vec_drphi_err_negz[index_r][index_z].push_back(sqrt(erp));
        vec_dz_err_negz[index_r][index_z].push_back(sqrt(ez));
        vec_tanAlpha_err_negz[index_r][index_z].push_back(0);
        vec_tanBeta_err_negz[index_r][index_z].push_back(0);
      }
    }

    TF1 *fit_drphi_tanAlpha_posz[nRBins][nZBins];
    TF1 *fit_drphi_tanAlpha_negz[nRBins][nZBins];
    TF1 *fit_dz_tanBeta_posz[nRBins][nZBins];
    TF1 *fit_dz_tanBeta_negz[nRBins][nZBins];
    for (int i = 0; i < nRBins; i++)
    {
        for (int j = 0; j < nZBins; j++)
        {
            gr_drphi_tanAlpha_posz[i][j] = new TGraphErrors(vec_drphi_posz[i][j].size(), &vec_tanAlpha_posz[i][j][0], &vec_drphi_posz[i][j][0], &vec_tanAlpha_err_posz[i][j][0], &vec_drphi_err_posz[i][j][0]);
            gr_drphi_tanAlpha_posz[i][j]->SetNameTitle(Form("gr_drphi_tanAlpha_posz_R%d_Z%d", i, j),
                              Form("R#in[%.1f,%.1f), Z#in[%.1f,%.1f);tan#alpha;drphi (cm)", 
                                  rMin + i*rWidth, rMin + (i+1)*rWidth, 
                                  zMin + j*zWidth, zMin + (j+1)*zWidth));
            gr_drphi_tanAlpha_posz[i][j]->SetMarkerStyle(20);
            gr_drphi_tanAlpha_posz[i][j]->SetMarkerSize(0.2);

            gr_drphi_tanAlpha_negz[i][j] = new TGraphErrors(vec_drphi_negz[i][j].size(), &vec_tanAlpha_negz[i][j][0], &vec_drphi_negz[i][j][0], &vec_tanAlpha_err_negz[i][j][0], &vec_drphi_err_negz[i][j][0]);
            gr_drphi_tanAlpha_negz[i][j]->SetNameTitle(Form("gr_drphi_tanAlpha_negz_R%d_Z%d", i, j),
                              Form("R#in[%.1f,%.1f), Z#in[%.1f,%.1f);tan#alpha;drphi (cm)", 
                                  -1*(rMin + (i+1)*rWidth), -1*(rMin + i*rWidth), 
                                  -1*(zMin + (j+1)*zWidth), -1*(zMin + j*zWidth)));
            gr_drphi_tanAlpha_negz[i][j]->SetMarkerStyle(20);
            gr_drphi_tanAlpha_negz[i][j]->SetMarkerSize(0.2);

            gr_dz_tanBeta_posz[i][j] = new TGraphErrors(vec_dz_posz[i][j].size(), &vec_tanBeta_posz[i][j][0], &vec_dz_posz[i][j][0], &vec_tanBeta_err_posz[i][j][0], &vec_dz_err_posz[i][j][0]);
            gr_dz_tanBeta_posz[i][j]->SetNameTitle(Form("gr_dz_tanBeta_posz_R%d_Z%d", i, j),
                              Form("R#in[%.1f,%.1f), Z#in[%.1f,%.1f);tan#beta;dz (cm)", 
                                  rMin + i*rWidth, rMin + (i+1)*rWidth, 
                                  zMin + j*zWidth, zMin + (j+1)*zWidth));
            gr_dz_tanBeta_posz[i][j]->SetMarkerStyle(20);
            gr_dz_tanBeta_posz[i][j]->SetMarkerSize(0.2);

            gr_dz_tanBeta_negz[i][j] = new TGraphErrors(vec_dz_negz[i][j].size(), &vec_tanBeta_negz[i][j][0], &vec_dz_negz[i][j][0], &vec_tanBeta_err_negz[i][j][0], &vec_dz_err_negz[i][j][0]);
            gr_dz_tanBeta_negz[i][j]->SetNameTitle(Form("gr_dz_tanBeta_negz_R%d_Z%d", i, j),
                              Form("R#in[%.1f,%.1f), Z#in[%.1f,%.1f);tan#beta;dz (cm)", 
                                  -1*(rMin + (i+1)*rWidth), -1*(rMin + i*rWidth),
                                  -1*(zMin + (j+1)*zWidth), -1*(zMin + j*zWidth)));
            gr_dz_tanBeta_negz[i][j]->SetMarkerStyle(20);
            gr_dz_tanBeta_negz[i][j]->SetMarkerSize(0.2);

            // Create a linear function for fitting: y = a*x + b

            fit_drphi_tanAlpha_posz[i][j] = new TF1(Form("fit_drphi_tanAlpha_posz_R%d_Z%d",i,j), "[0]*x + [1]", -0.5, 0.5);
            fit_drphi_tanAlpha_posz[i][j]->SetParNames("Slope", "Intercept");
            fit_drphi_tanAlpha_posz[i][j]->SetParameters(0, 0); // Initial parameter guesses
            fit_drphi_tanAlpha_posz[i][j]->SetLineColor(kRed);

            fit_drphi_tanAlpha_negz[i][j] = new TF1(Form("fit_drphi_tanAlpha_negz_R%d_Z%d",i,j), "[0]*x + [1]", -0.5, 0.5);
            fit_drphi_tanAlpha_negz[i][j]->SetParNames("Slope", "Intercept");
            fit_drphi_tanAlpha_negz[i][j]->SetParameters(0, 0); // Initial parameter guesses
            fit_drphi_tanAlpha_negz[i][j]->SetLineColor(kRed);

            fit_dz_tanBeta_posz[i][j] = new TF1(Form("fit_dz_tanBeta_posz_R%d_Z%d",i,j), "[0]*x + [1]", -0.5, 0.5);
            fit_dz_tanBeta_posz[i][j]->SetParNames("Slope", "Intercept");
            fit_dz_tanBeta_posz[i][j]->SetParameters(0, 0); // Initial parameter guesses
            fit_dz_tanBeta_posz[i][j]->SetLineColor(kRed);

            fit_dz_tanBeta_negz[i][j] = new TF1(Form("fit_dz_tanBeta_negz_R%d_Z%d",i,j), "[0]*x + [1]", -0.5, 0.5);
            fit_dz_tanBeta_negz[i][j]->SetParNames("Slope", "Intercept");
            fit_dz_tanBeta_negz[i][j]->SetParameters(0, 0); // Initial parameter guesses
            fit_dz_tanBeta_negz[i][j]->SetLineColor(kRed);

            // Perform the fit (equal weights for all points)
            gr_drphi_tanAlpha_posz[i][j]->Fit(Form("fit_drphi_tanAlpha_posz_R%d_Z%d",i,j), "Q"); // "Q" for quiet mode
            gr_drphi_tanAlpha_negz[i][j]->Fit(Form("fit_drphi_tanAlpha_negz_R%d_Z%d",i,j), "Q"); // "Q" for quiet mode
            gr_dz_tanBeta_posz[i][j]->Fit(Form("fit_dz_tanBeta_posz_R%d_Z%d",i,j), "Q"); // "Q" for quiet mode
            gr_dz_tanBeta_negz[i][j]->Fit(Form("fit_dz_tanBeta_negz_R%d_Z%d",i,j), "Q"); // "Q" for quiet mode

            if (verbosity>1)
            {
              TCanvas *can_drphi_tanAlpha_posz = new TCanvas("can_drphi_tanAlpha_posz","",800,600);
              gr_drphi_tanAlpha_posz[i][j]->SetMarkerStyle(20);
              gr_drphi_tanAlpha_posz[i][j]->SetMarkerColor(kBlack);
              gr_drphi_tanAlpha_posz[i][j]->Draw("AP");
              fit_drphi_tanAlpha_posz[i][j]->Draw("same");
              can_drphi_tanAlpha_posz->Update();
              can_drphi_tanAlpha_posz->SaveAs(Form("figure/%d_drphi_tanAlpha_posz_R%d_Z%d.pdf",runs[k],i,j));

              TCanvas *can_drphi_tanAlpha_negz = new TCanvas("can_drphi_tanAlpha_negz","",800,600);
              gr_drphi_tanAlpha_negz[i][j]->SetMarkerStyle(20);
              gr_drphi_tanAlpha_negz[i][j]->SetMarkerColor(kBlack);
              gr_drphi_tanAlpha_negz[i][j]->Draw("AP");
              fit_drphi_tanAlpha_negz[i][j]->Draw("same");
              can_drphi_tanAlpha_negz->Update();
              can_drphi_tanAlpha_negz->SaveAs(Form("figure/%d_drphi_tanAlpha_negz_R%d_Z%d.pdf",runs[k],i,j));

              TCanvas *can_dz_tanBeta_posz = new TCanvas("can_dz_tanBeta_posz","",800,600);
              gr_dz_tanBeta_posz[i][j]->SetMarkerStyle(20);
              gr_dz_tanBeta_posz[i][j]->SetMarkerColor(kBlack);
              gr_dz_tanBeta_posz[i][j]->Draw("AP");
              fit_dz_tanBeta_posz[i][j]->Draw("same");
              can_dz_tanBeta_posz->Update();
              can_dz_tanBeta_posz->SaveAs(Form("figure/%d_dz_tanBeta_posz_R%d_Z%d.pdf",runs[k],i,j));

              TCanvas *can_dz_tanBeta_negz = new TCanvas("can_dz_tanBeta_negz","",800,600);
              gr_dz_tanBeta_negz[i][j]->SetMarkerStyle(20);
              gr_dz_tanBeta_negz[i][j]->SetMarkerColor(kBlack);
              gr_dz_tanBeta_negz[i][j]->Draw("AP");
              fit_dz_tanBeta_negz[i][j]->Draw("same");
              can_dz_tanBeta_negz->Update();
              can_dz_tanBeta_negz->SaveAs(Form("figure/%d_dz_tanBeta_negz_R%d_Z%d.pdf",runs[k],i,j));

            }
        }
    }

    TH2* h2_dr_rphi_posz = new TH2D("h2_dr_rphi_posz","dR distortion from R-#phi plane (posz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,zMin,zMax);
    TH2* h2_dr_rphi_negz = new TH2D("h2_dr_rphi_negz","dR distortion from R-#phi plane (negz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,-1*zMax,-1*zMin);
    TH2* h2_rdphi_rphi_posz = new TH2D("h2_rdphi_rphi_posz","Rdphi distortion from R-#phi plane (posz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,zMin,zMax);
    TH2* h2_rdphi_rphi_negz = new TH2D("h2_rdphi_rphi_negz","Rdphi distortion from R-#phi plane (negz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,-1*zMax,-1*zMin);
    TH2* h2_dr_zr_posz = new TH2D("h2_dr_zr_posz","dR distortion from Z-R plane (posz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,zMin,zMax);
    TH2* h2_dr_zr_negz = new TH2D("h2_dr_zr_negz","dR distortion from Z-R plane (negz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,-1*zMax,-1*zMin);
    TH2* h2_dz_zr_posz = new TH2D("h2_dz_zr_posz","dz distortion from Z-R plane (posz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,zMin,zMax);
    TH2* h2_dz_zr_negz = new TH2D("h2_dz_zr_negz","dz distortion from Z-R plane (negz);R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,-1*zMax,-1*zMin);
    TH2* h2_counts_posz = new TH2D("h2_counts_posz","Statistics;R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,zMin,zMax);
    TH2* h2_counts_negz = new TH2D("h2_counts_negz","Statistics;R (cm);Z (cm)",nRBins,rMin,rMax,nZBins,-1*zMax,-1*zMin);
    for (int i = 1; i <= nRBins; i++)
    {
        for (int j = 1; j <= nZBins; j++)
        {
            h2_dr_rphi_posz->SetBinContent(i,j,fit_drphi_tanAlpha_posz[i-1][j-1]->GetParameter(0));
            h2_dr_rphi_posz->SetBinError(i,j,fit_drphi_tanAlpha_posz[i-1][j-1]->GetParError(0));

            h2_dr_rphi_negz->SetBinContent(i,j,fit_drphi_tanAlpha_negz[i-1][nZBins-j]->GetParameter(0));
            h2_dr_rphi_negz->SetBinError(i,j,fit_drphi_tanAlpha_negz[i-1][nZBins-j]->GetParError(0));

            h2_rdphi_rphi_posz->SetBinContent(i,j,fit_drphi_tanAlpha_posz[i-1][j-1]->GetParameter(1));
            h2_rdphi_rphi_posz->SetBinError(i,j,fit_drphi_tanAlpha_posz[i-1][j-1]->GetParError(1));

            h2_rdphi_rphi_negz->SetBinContent(i,j,fit_drphi_tanAlpha_negz[i-1][nZBins-j]->GetParameter(1));
            h2_rdphi_rphi_negz->SetBinError(i,j,fit_drphi_tanAlpha_negz[i-1][nZBins-j]->GetParError(1));

            h2_dr_zr_posz->SetBinContent(i,j,fit_dz_tanBeta_posz[i-1][j-1]->GetParameter(0));
            h2_dr_zr_posz->SetBinError(i,j,fit_dz_tanBeta_posz[i-1][j-1]->GetParError(0));

            h2_dr_zr_negz->SetBinContent(i,j,fit_dz_tanBeta_negz[i-1][nZBins-j]->GetParameter(0));
            h2_dr_zr_negz->SetBinError(i,j,fit_dz_tanBeta_negz[i-1][nZBins-j]->GetParError(0));

            h2_dz_zr_posz->SetBinContent(i,j,fit_dz_tanBeta_posz[i-1][j-1]->GetParameter(1));
            h2_dz_zr_posz->SetBinError(i,j,fit_dz_tanBeta_posz[i-1][j-1]->GetParError(1));

            h2_dz_zr_negz->SetBinContent(i,j,fit_dz_tanBeta_negz[i-1][nZBins-j]->GetParameter(1));
            h2_dz_zr_negz->SetBinError(i,j,fit_dz_tanBeta_negz[i-1][nZBins-j]->GetParError(1));

            h2_counts_posz->SetBinContent(i,j,vec_drphi_posz[i-1][j-1].size());
            h2_counts_negz->SetBinContent(i,j,vec_drphi_negz[i-1][nZBins-j].size());
        }
    }

    TCanvas *can_rphi = new TCanvas("can_rphi","",1600,1200);
    can_rphi->Divide(2,2);
    can_rphi->cd(1);
    h2_dr_rphi_posz->Draw("colz");
    can_rphi->cd(2);
    h2_rdphi_rphi_posz->Draw("colz");
    can_rphi->cd(3);
    h2_dr_rphi_negz->Draw("colz");
    can_rphi->cd(4);
    h2_rdphi_rphi_negz->Draw("colz");
    can_rphi->Update();
    can_rphi->Print(Form("figure/%d_dr_rphi.pdf",runs[k]));

    TCanvas *can_zr = new TCanvas("can_zr","",800,1200);
    can_zr->Divide(1,2);
    can_zr->cd(1);
    h2_dr_zr_posz->Draw("colz");
    can_zr->cd(2);
    h2_dz_zr_posz->Draw("colz");
    can_zr->cd(3);
    h2_dr_zr_negz->Draw("colz");
    can_zr->cd(4);
    h2_dz_zr_negz->Draw("colz");
    can_zr->Update();
    can_zr->Print(Form("figure/%d_dr_zr.pdf",runs[k]));

    TCanvas *can_counts = new TCanvas("can_counts","",800,1200);
    can_counts->Divide(1,2);
    can_counts->cd(1);
    h2_counts_posz->Draw("colz");
    can_counts->cd(2);
    h2_counts_negz->Draw("colz");
    can_counts->Update();
    can_counts->Print(Form("figure/%d_counts.pdf",runs[k]));

    TFile* ofile = new TFile(Form("hist_residual_2Drz_%d.root",runs[k]),"recreate");
    ofile->cd();
    h2_dr_rphi_posz->Write();
    h2_dr_rphi_negz->Write();
    h2_dr_zr_posz->Write();
    h2_dr_zr_negz->Write();
    h2_rdphi_rphi_posz->Write();
    h2_rdphi_rphi_negz->Write();
    h2_dz_zr_posz->Write();
    h2_dz_zr_negz->Write();
    h2_counts_posz->Write();
    h2_counts_negz->Write();
    ofile->Close();

  }
}
