void plot1D_Zbin(TH3* h3, TH1* h1, float z, TString phiType)
{
//  TpcSpaceChargeReconstructionHelper::set_phi_range_central( {-1.73246,-1.43608} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_east( {-2.26272,-1.96089} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_west( {-1.21241,-0.909953} );

  int xbin = 0;
  if (phiType.CompareTo("central")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-1.73246)+(-1.43608))/2.+2*TMath::Pi());
  }
  else if (phiType.CompareTo("east")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-2.26272)+(-1.96089))/2.+2*TMath::Pi());
  }
  else if (phiType.CompareTo("west")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-1.21241)+(-0.909953))/2.+2*TMath::Pi());
  }

  int zbin = h3->GetZaxis()->FindBin(z);

  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      h1->SetBinContent(i, h3->GetBinContent(xbin, i, zbin));
      h1->SetBinError(i, h3->GetBinError(xbin, i, zbin));
  }
}

void plot1D_Ybin(TH3* h3, TH1* h1, float y, TString phiType)
{
//  TpcSpaceChargeReconstructionHelper::set_phi_range_central( {-1.73246,-1.43608} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_east( {-2.26272,-1.96089} );
//  TpcSpaceChargeReconstructionHelper::set_phi_range_west( {-1.21241,-0.909953} );

  int xbin = 0;
  if (phiType.CompareTo("central")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-1.73246)+(-1.43608))/2.+2*TMath::Pi());
  }
  else if (phiType.CompareTo("east")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-2.26272)+(-1.96089))/2.+2*TMath::Pi());
  }
  else if (phiType.CompareTo("west")==0)
  {
    xbin = h3->GetXaxis()->FindBin(((-1.21241)+(-0.909953))/2.+2*TMath::Pi());
  }

  int ybin = h3->GetYaxis()->FindBin(y);

  for (int i = 1; i <= h3->GetNbinsZ(); i++)
  {
      h1->SetBinContent(i, h3->GetBinContent(xbin, ybin, i));
      h1->SetBinError(i, h3->GetBinError(xbin, ybin, i));
  }
}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin>0) {yMin*=0.9;}
    else if (yMin<0) {yMin*=1.1;}
    if (yMax>0) {yMax*=1.1;}
    else if (yMax<0) {yMax*=0.9;}

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
