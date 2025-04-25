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

  TH1* h_dr_layer_pos = new TH1F("h_dr_layer_pos","Positive Z Simulation Distortion Map dr vs layer;layer;dR (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dphi_layer_pos = new TH1F("h_dphi_layer_pos","Positive Z Simulation Distortion Map dphi vs layer;layer;dphi (rad)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_rdphi_layer_pos = new TH1F("h_rdphi_layer_pos","Positive Z Simulation Distortion Map rdphi vs layer;layer;rdphi (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dz_layer_pos = new TH1F("h_dz_layer_pos","Positive Z Simulation Distortion Map dz vs layer;layer;dz (cm)",nbins_layer,xmin_layer,xmax_layer);

  TH1* h_dr_layer_neg = new TH1F("h_dr_layer_neg","Negative Z Simulation Distortion Map dr vs layer;layer;dR (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dphi_layer_neg = new TH1F("h_dphi_layer_neg","Negative Z Simulation Distortion Map dphi vs layer;layer;dphi (rad)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_rdphi_layer_neg = new TH1F("h_rdphi_layer_neg","Negative Z Simulation Distortion Map rdphi vs layer;layer;rdphi (cm)",nbins_layer,xmin_layer,xmax_layer);
  TH1* h_dz_layer_neg = new TH1F("h_dz_layer_neg","Negative Z Simulation Distortion Map dz vs layer;layer;dz (cm)",nbins_layer,xmin_layer,xmax_layer);

  TH1* h_dr_r_pos = new TH1F("h_dr_r_pos","Positive Z Simulation Distortion Map dr vs r;r (cm);dR (cm)",nbins_r,xmin_r,xmax_r);
  TH1* h_dphi_r_pos = new TH1F("h_dphi_r_pos","Positive Z Simulation Distortion Map dphi vs r;r;dphi (rad)",nbins_r,xmin_r,xmax_r);
  TH1* h_rdphi_r_pos = new TH1F("h_rdphi_r_pos","Positive Z Simulation Distortion Map rdphi vs r;r (cm);rdphi (cm)",nbins_r,xmin_r,xmax_r);
  TH1* h_dz_r_pos = new TH1F("h_dz_r_pos","Positive Z Simulation Distortion Map dz vs r;r (cm);dz (cm)",nbins_r,xmin_r,xmax_r);

  TH1* h_dr_r_neg = new TH1F("h_dr_r_neg","Negative Z Simulation Distortion Map dr vs r;r (cm);dR (cm)",nbins_r,xmin_r,xmax_r);
  TH1* h_dphi_r_neg = new TH1F("h_dphi_r_neg","Negative Z Simulation Distortion Map dphi vs r;r (cm);dphi (rad)",nbins_r,xmin_r,xmax_r);
  TH1* h_rdphi_r_neg = new TH1F("h_rdphi_r_neg","Negative Z Simulation Distortion Map rdphi vs r;r (cm);rdphi (cm)",nbins_r,xmin_r,xmax_r);
  TH1* h_dz_r_neg = new TH1F("h_dz_r_neg","Negative Z Simulation Distortion Map dz vs r;r (cm);dz (cm)",nbins_r,xmin_r,xmax_r);

  for (int i = 0; i < nbins_layer; i++)
  {
    int layer = h_dr_layer_pos->GetXaxis()->GetBinCenter(i+1);
    float r = TpcRadiusMap[layer];
    h_dphi_layer_pos->SetBinContent(i+1,0);
    h_rdphi_layer_pos->SetBinContent(i+1,0);
    h_dr_layer_pos->SetBinContent(i+1,0);
    h_dz_layer_pos->SetBinContent(i+1,0);

    h_dphi_layer_neg->SetBinContent(i+1,0);
    h_rdphi_layer_neg->SetBinContent(i+1,0);
    h_dr_layer_neg->SetBinContent(i+1,0);
    h_dz_layer_neg->SetBinContent(i+1,0);
  }

  for (int i = 0; i < nbins_r; i++)
  {
    float r = h_dr_r_pos->GetXaxis()->GetBinCenter(i+1);
    h_dphi_r_pos->SetBinContent(i+1,0);
    h_rdphi_r_pos->SetBinContent(i+1,0);
    h_dr_r_pos->SetBinContent(i+1,0);
    h_dz_r_pos->SetBinContent(i+1,0);

    h_dphi_r_neg->SetBinContent(i+1,0);
    h_rdphi_r_neg->SetBinContent(i+1,0);
    h_dr_r_neg->SetBinContent(i+1,0);
    h_dz_r_neg->SetBinContent(i+1,0);
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
