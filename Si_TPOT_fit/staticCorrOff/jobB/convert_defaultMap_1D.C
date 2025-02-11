void loadTpcRadius();
float get_resid_r(float r, TH3* h, TH3* scaler_phi_r_z);
bool check_boundaries(const TAxis* axis, double value);
bool check_boundaries(const TH2* h, double r, double phi);
bool check_boundaries(const TH3* h, double r, double phi, double z);
std::map<int, float> TpcRadiusMap;

void convert_defaultMap_1D()
{
  loadTpcRadius();
  TGaxis::SetMaxDigits(3);

  std::string static_correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/static_only_inverted_10-new.root";
  TFile* file_map = new TFile(static_correction_filename.c_str(),"");
  TH3* hIntDistortionP_negz = (TH3F*) file_map->Get("hIntDistortionP_negz");
  TH3* hIntDistortionR_negz = (TH3F*) file_map->Get("hIntDistortionR_negz");
  TH3* hIntDistortionZ_negz = (TH3F*) file_map->Get("hIntDistortionZ_negz");
  TH3* hIntDistortionP_posz = (TH3F*) file_map->Get("hIntDistortionP_posz");
  TH3* hIntDistortionR_posz = (TH3F*) file_map->Get("hIntDistortionR_posz");
  TH3* hIntDistortionZ_posz = (TH3F*) file_map->Get("hIntDistortionZ_posz");

  // reconstructed distortion grid size (phi, r, z)
  float nbins_phi = hIntDistortionP_negz->GetXaxis()->GetNbins();
  float xmin_phi = hIntDistortionP_negz->GetXaxis()->GetXmin();
  float xmax_phi = hIntDistortionP_negz->GetXaxis()->GetXmax();
  float nbins_r = hIntDistortionP_negz->GetYaxis()->GetNbins();
  float xmin_r = hIntDistortionP_negz->GetYaxis()->GetXmin();
  float xmax_r = hIntDistortionP_negz->GetYaxis()->GetXmax();
  float nbins_z = hIntDistortionP_negz->GetZaxis()->GetNbins();
  float xmin_z = hIntDistortionP_negz->GetZaxis()->GetXmin();
  float xmax_z = hIntDistortionP_posz->GetZaxis()->GetXmax();
  float nbins_layer = 48;
  float xmin_layer = 7;
  float xmax_layer = 55;

  cout<<"X-axis nbins "<<nbins_phi<<" min "<<xmin_phi<<" max "<<xmax_phi<<endl;
  cout<<"Y-axis nbins "<<nbins_r<<" min "<<xmin_r<<" max "<<xmax_r<<endl;
  cout<<"Z-axis nbins "<<nbins_z<<" min "<<xmin_z<<" max "<<xmax_z<<endl;

  bool read_scaler_from_hist = true;
  TH1* h_phi;
  TH1* h_z;
  TH1* h_r;
  TH1* h_layer;
  TH2* h_phi_z;
  TH2* h_phi_r;
  TH2* h_phi_layer;
  TH2* h_r_z;
  TH2* h_layer_z;
  TH3* h_phi_r_z;
  TH3* h_phi_layer_z;

  if (read_scaler_from_hist)
  {
    TFile* file_scale = new TFile("scale.root","");
    h_phi = (TH1*) file_scale->Get("h_phi");
    h_z = (TH1*) file_scale->Get("h_z");
    h_r = (TH1*) file_scale->Get("h_r");
    h_layer = (TH1*) file_scale->Get("h_layer");
    h_phi_z = (TH2*) file_scale->Get("h_phi_z");
    h_phi_r = (TH2*) file_scale->Get("h_phi_r");
    h_phi_layer = (TH2*) file_scale->Get("h_phi_layer");
    h_r_z = (TH2*) file_scale->Get("h_r_z");
    h_layer_z = (TH2*) file_scale->Get("h_layer_z");
    h_phi_r_z = (TH3*) file_scale->Get("h_phi_r_z");
    h_phi_layer_z = (TH3*) file_scale->Get("h_phi_layer_z");
  }
  else
  {
    TFile* file_scale = new TFile("/sphenix/u/xyu3/workarea/TPCana/Distortion_lumi/edge/clusters_seeds_53877_all.root","");
    TTree* tree_scale = (TTree*) file_scale->Get("residualtree");
    h_phi = new TH1F("h_phi", "Run 53877, clusters on track;#phi (rad);counts", nbins_phi, xmin_phi, xmax_phi);
    h_z = new TH1F("h_z", "Run 53877, clusters on track;Z (cm);Counts", nbins_z, xmin_z, xmax_z);
    h_r = new TH1F("h_r", "Run 53877, clusters on track;R (cm);Counts", nbins_r, xmin_r, xmax_r);
    h_layer = new TH1F("h_layer", "Run 53877, clusters on track;Layer;Counts", nbins_layer, xmin_layer, xmax_layer);
    h_phi_z = new TH2F("h_phi_z", "Run 53877, clusters on track;#phi (rad);Z (cm);Counts", nbins_phi, xmin_phi, xmax_phi, nbins_z, xmin_z, xmax_z);
    h_phi_r = new TH2F("h_phi_r", "Run 53877, clusters on track;#phi (rad);R (cm);Counts", nbins_phi, xmin_phi, xmax_phi, nbins_r, xmin_r, xmax_r);
    h_phi_layer = new TH2F("h_phi_layer", "Run 53877, clusters on track;#phi (rad);Layer;Counts", nbins_phi, xmin_phi, xmax_phi, nbins_layer, xmin_layer, xmax_layer);
    h_r_z = new TH2F("h_r_z", "Run 53877, clusters on track;R (cm);Z (cm);Counts", nbins_r, xmin_r, xmax_r, nbins_z, xmin_z, xmax_z);
    h_layer_z = new TH2F("h_layer_z", "Run 53877, clusters on track;Layer;Z (cm);Counts", nbins_layer, xmin_layer, xmax_layer, nbins_z, xmin_z, xmax_z);
    h_phi_r_z = new TH3F("h_phi_r_z", "Run 53877, clusters on track;#phi (rad);R (cm);Z (cm);Counts", nbins_phi, xmin_phi, xmax_phi, nbins_r, xmin_r, xmax_r, nbins_z, xmin_z, xmax_z);
    h_phi_layer_z = new TH3F("h_phi_layer_z", "Run 53877, clusters on track;#phi (rad);Layer;Z (cm);Counts", nbins_phi, xmin_phi, xmax_phi, nbins_layer, xmin_layer, xmax_layer, nbins_z, xmin_z, xmax_z);

    tree_scale->Draw("clusgz>>h_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("fabs(clusgr)>>h_r","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("cluslayer>>h_layer","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("clusgz:(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("fabs(clusgr):(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi_r","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("cluslayer:(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi_layer","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("clusgz:fabs(clusgr)>>h_r_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("clusgz:cluslayer>>h_layer_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("clusgz:fabs(clusgr):(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi_r_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");
    tree_scale->Draw("clusgz:cluslayer:(atan2(clusgy,clusgx) < 0 ? atan2(clusgy,clusgx)+2*TMath::Pi() : atan2(clusgy,clusgx))>>h_phi_layer_z","cluslayer>=7 && cluslayer<55 && m_pt>0.2 && m_nmaps>=2 && m_nintt>=2");

    h_z->Scale(1./h_z->Integral());
    h_phi->Scale(1./h_phi->Integral());
    h_r->Scale(1./h_r->Integral());
    h_layer->Scale(1./h_layer->Integral());
    h_phi_z->Scale(1./h_phi_z->Integral());
    h_phi_r->Scale(1./h_phi_r->Integral());
    h_phi_layer->Scale(1./h_phi_layer->Integral());
    h_r_z->Scale(1./h_r_z->Integral());
    h_layer_z->Scale(1./h_layer_z->Integral());
    h_phi_r_z->Scale(1./h_phi_r_z->Integral());
    h_phi_layer_z->Scale(1./h_phi_layer_z->Integral());

    TFile* ofile_scale = new TFile("scale.root","recreate");
    ofile_scale->cd();
    h_phi->Write();
    h_z->Write();
    h_r->Write();
    h_layer->Write();
    h_phi_z->Write();
    h_phi_r->Write();
    h_phi_layer->Write();
    h_r_z->Write();
    h_layer_z->Write();
    h_phi_r_z->Write();
    h_phi_layer_z->Write();
    ofile_scale->Close();
  }

  TCanvas* can_scale = new TCanvas("can_scale","",2400,1800);
  can_scale->Divide(3,3);
  can_scale->cd(1);
  h_z->Draw("hist");
  can_scale->cd(2);
  h_phi->Draw("hist");
  can_scale->cd(3);
  h_r->Draw("hist");
  can_scale->cd(4);
  h_layer->Draw("hist");
  can_scale->cd(5);
  h_phi_z->Draw("colz");
  can_scale->cd(6);
  h_phi_r->Draw("colz");
  can_scale->cd(7);
  h_phi_layer->Draw("colz");
  can_scale->cd(8);
  h_r_z->Draw("colz");
  can_scale->cd(9);
  h_layer_z->Draw("colz");
  can_scale->SaveAs("./figure/scale.pdf");

  TH1* h_dr_layer_pos = new TH1F("h_dr_layer_pos","Positive Z Simulation Distortion Map dr vs layer;layer;dR (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dphi_layer_pos = new TH1F("h_dphi_layer_pos","Positive Z Simulation Distortion Map dphi vs layer;layer;dphi (rad)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_rdphi_layer_pos = new TH1F("h_rdphi_layer_pos","Positive Z Simulation Distortion Map rdphi vs layer;layer;rdphi (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dz_layer_pos = new TH1F("h_dz_layer_pos","Positive Z Simulation Distortion Map dz vs layer;layer;dz (cm)",nbins_layer,xmin_layer,xmax_layer);

  TH1* h_dr_layer_neg = new TH1F("h_dr_layer_neg","Negative Z Simulation Distortion Map dr vs layer;layer;dR (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dphi_layer_neg = new TH1F("h_dphi_layer_neg","Negative Z Simulation Distortion Map dphi vs layer;layer;dphi (rad)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_rdphi_layer_neg = new TH1F("h_rdphi_layer_neg","Negative Z Simulation Distortion Map rdphi vs layer;layer;rdphi (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dz_layer_neg = new TH1F("h_dz_layer_neg","Negative Z Simulation Distortion Map dz vs layer;layer;dz (cm)",nbins_layer,xmin_layer,xmax_layer);

  TH1* h_dr_r_pos = new TH1F("h_dr_r_pos","Positive Z Simulation Distortion Map dr vs r;r (cm);dR (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_dphi_r_pos = new TH1F("h_dphi_r_pos","Positive Z Simulation Distortion Map dphi vs r;r;dphi (rad)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_rdphi_r_pos = new TH1F("h_rdphi_r_pos","Positive Z Simulation Distortion Map rdphi vs r;r (cm);rdphi (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_dz_r_pos = new TH1F("h_dz_r_pos","Positive Z Simulation Distortion Map dz vs r;r (cm);dz (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);

  TH1* h_dr_r_neg = new TH1F("h_dr_r_neg","Negative Z Simulation Distortion Map dr vs r;r (cm);dR (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_dphi_r_neg = new TH1F("h_dphi_r_neg","Negative Z Simulation Distortion Map dphi vs r;r (cm);dphi (rad)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_rdphi_r_neg = new TH1F("h_rdphi_r_neg","Negative Z Simulation Distortion Map rdphi vs r;r (cm);rdphi (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);
  TH1* h_dz_r_neg = new TH1F("h_dz_r_neg","Negative Z Simulation Distortion Map dz vs r;r (cm);dz (cm)",nbins_r,TpcRadiusMap[7],TpcRadiusMap[54]);

  for (int i = 0; i < nbins_layer; i++)
  {
    int layer = h_dr_layer_pos->GetXaxis()->GetBinCenter(i+1);
    float r = TpcRadiusMap[layer];
    h_dphi_layer_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_posz,h_phi_r_z)/r);
    h_rdphi_layer_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_posz,h_phi_r_z));
    h_dr_layer_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionR_posz,h_phi_r_z));
    h_dz_layer_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionZ_posz,h_phi_r_z));

    h_dphi_layer_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_negz,h_phi_r_z)/r);
    h_rdphi_layer_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_negz,h_phi_r_z));
    h_dr_layer_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionR_negz,h_phi_r_z));
    h_dz_layer_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionZ_negz,h_phi_r_z));
  }

  for (int i = 0; i < nbins_r; i++)
  {
    float r = h_dr_r_pos->GetXaxis()->GetBinCenter(i+1);
    h_dphi_r_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_posz,h_phi_r_z)/r);
    h_rdphi_r_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_posz,h_phi_r_z));
    h_dr_r_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionR_posz,h_phi_r_z));
    h_dz_r_pos->SetBinContent(i+1,get_resid_r(r,hIntDistortionZ_posz,h_phi_r_z));

    h_dphi_r_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_negz,h_phi_r_z)/r);
    h_rdphi_r_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionP_negz,h_phi_r_z));
    h_dr_r_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionR_negz,h_phi_r_z));
    h_dz_r_neg->SetBinContent(i+1,get_resid_r(r,hIntDistortionZ_negz,h_phi_r_z));
  }

  TCanvas* can_map = new TCanvas("can_map","",3200,1200);
  can_map->Divide(4,2);
  can_map->cd(1);
  h_dphi_layer_pos->Draw("hist");
  can_map->cd(2);
  h_rdphi_layer_pos->Draw("hist");
  can_map->cd(3);
  h_dphi_layer_neg->Draw("hist");
  can_map->cd(4);
  h_rdphi_layer_neg->Draw("hist");
  can_map->cd(5);
  h_dr_layer_pos->Draw("hist");
  can_map->cd(6);
  h_dz_layer_pos->Draw("hist");
  can_map->cd(7);
  h_dr_layer_neg->Draw("hist");
  can_map->cd(8);
  h_dz_layer_neg->Draw("hist");
  can_map->SaveAs("./figure/defaultMap1D.pdf");

  TCanvas* can_map_r = new TCanvas("can_map_r","",3200,1200);
  can_map_r->Divide(4,2);
  can_map_r->cd(1);
  h_dphi_r_pos->Draw("hist");
  can_map_r->cd(2);
  h_rdphi_r_pos->Draw("hist");
  can_map_r->cd(3);
  h_dphi_r_neg->Draw("hist");
  can_map_r->cd(4);
  h_rdphi_r_neg->Draw("hist");
  can_map_r->cd(5);
  h_dr_r_pos->Draw("hist");
  can_map_r->cd(6);
  h_dz_r_pos->Draw("hist");
  can_map_r->cd(7);
  h_dr_r_neg->Draw("hist");
  can_map_r->cd(8);
  h_dz_r_neg->Draw("hist");
  can_map_r->SaveAs("./figure/defaultMap1D_r.pdf");

  TFile* ofile_map = new TFile("default_map.root","recreate");
  ofile_map->cd();
  h_dphi_layer_pos->Write();
  h_rdphi_layer_pos->Write();
  h_dr_layer_pos->Write();
  h_dz_layer_pos->Write();
  h_dphi_layer_neg->Write();
  h_rdphi_layer_neg->Write();
  h_dr_layer_neg->Write();
  h_dz_layer_neg->Write();

  h_dphi_r_pos->Write();
  h_rdphi_r_pos->Write();
  h_dr_r_pos->Write();
  h_dz_r_pos->Write();
  h_dphi_r_neg->Write();
  h_rdphi_r_neg->Write();
  h_dr_r_neg->Write();
  h_dz_r_neg->Write();

  ofile_map->Close();

}

float get_resid_r(float r, TH3* h, TH3* scaler_phi_r_z)
{
  float resid = 0;
  float nscale = 0;
  for (int iphi = 0; iphi < (h->GetXaxis()->GetNbins()); iphi++)
  {
    float phi = h->GetXaxis()->GetBinCenter(iphi+1);
    for (int iz = 0; iz < (h->GetZaxis()->GetNbins()); iz++)
    {
      float z = h->GetZaxis()->GetBinCenter(iz+1);
      if (check_boundaries(h, phi, r, z) && check_boundaries(scaler_phi_r_z, phi, r, z))
      {
        resid+=h->Interpolate(phi, r, z)*scaler_phi_r_z->Interpolate(phi,r,z);
        nscale+=scaler_phi_r_z->Interpolate(phi,r,z);
      }
    }
  }

  if (nscale==0) {return 0;}

  resid/=nscale;
  return resid;
}

void loadTpcRadius()
{
  TpcRadiusMap[7] = 31.372;
  TpcRadiusMap[8] = 31.944;
  TpcRadiusMap[9] = 32.5159;
  TpcRadiusMap[10] = 33.0879;
  TpcRadiusMap[11] = 33.6599;
  TpcRadiusMap[12] = 34.2318;
  TpcRadiusMap[13] = 34.8038;
  TpcRadiusMap[14] = 35.3758;
  TpcRadiusMap[15] = 35.9477;
  TpcRadiusMap[16] = 36.5197;
  TpcRadiusMap[17] = 37.0917;
  TpcRadiusMap[18] = 37.6636;
  TpcRadiusMap[19] = 38.2356;
  TpcRadiusMap[20] = 38.8076;
  TpcRadiusMap[21] = 39.3795;
  TpcRadiusMap[22] = 39.9515;
  TpcRadiusMap[23] = 41.659;
  TpcRadiusMap[24] = 42.6797;
  TpcRadiusMap[25] = 43.7003;
  TpcRadiusMap[26] = 44.721;
  TpcRadiusMap[27] = 45.7417;
  TpcRadiusMap[28] = 46.7623;
  TpcRadiusMap[29] = 47.783;
  TpcRadiusMap[30] = 48.8037;
  TpcRadiusMap[31] = 49.8243;
  TpcRadiusMap[32] = 50.845;
  TpcRadiusMap[33] = 51.8657;
  TpcRadiusMap[34] = 52.8863;
  TpcRadiusMap[35] = 53.907;
  TpcRadiusMap[36] = 54.9277;
  TpcRadiusMap[37] = 55.9483;
  TpcRadiusMap[38] = 56.969;
  TpcRadiusMap[39] = 58.911;
  TpcRadiusMap[40] = 60.0081;
  TpcRadiusMap[41] = 61.1051;
  TpcRadiusMap[42] = 62.2022;
  TpcRadiusMap[43] = 63.2993;
  TpcRadiusMap[44] = 64.3963;
  TpcRadiusMap[45] = 65.4934;
  TpcRadiusMap[46] = 66.5905;
  TpcRadiusMap[47] = 67.6875;
  TpcRadiusMap[48] = 68.7846;
  TpcRadiusMap[49] = 69.8817;
  TpcRadiusMap[50] = 70.9787;
  TpcRadiusMap[51] = 72.0758;
  TpcRadiusMap[52] = 73.1729;
  TpcRadiusMap[53] = 74.2699;
  TpcRadiusMap[54] = 75.367;
}

bool check_boundaries(const TAxis* axis, double value)
{
  const auto bin = axis->FindBin(value);
  return (bin >= 2 && bin < axis->GetNbins());
}

bool check_boundaries(const TH2* h, double r, double phi)
{
  return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi);
}

bool check_boundaries(const TH3* h, double r, double phi, double z)
{
  return check_boundaries(h->GetXaxis(), r) && check_boundaries(h->GetYaxis(), phi) && check_boundaries(h->GetZaxis(), z);
}
