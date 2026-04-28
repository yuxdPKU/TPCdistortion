std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void draw1D_r_runs()
{

const int nruns = 2;
//int runs[nruns] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
int runs[nruns] = {79507,79516};

int zpoint=10;
TFile* infile_P = new TFile(Form("../Rootfiles/hist_dphi_1Dr_Z%d.root",zpoint),"");
TFile* infile_R = new TFile(Form("../Rootfiles/hist_dr_1Dr_Z%d.root",zpoint),"");
TFile* infile_Z = new TFile(Form("../Rootfiles/hist_dz_1Dr_Z%d.root",zpoint),"");
TH1F *inhist_P_posz[nruns], *inhist_P_negz[nruns];
TH1F *inhist_R_posz[nruns], *inhist_R_negz[nruns];
TH1F *inhist_Z_posz[nruns], *inhist_Z_negz[nruns];

for (int i=0; i<nruns; i++)
{
	inhist_P_posz[i] = (TH1F*) infile_P->Get(Form("hIntDistortionP_posz_%d",runs[i]));
	inhist_P_posz[i]->SetLineColor(i+1);
	inhist_P_posz[i]->SetMarkerColor(i+1);
	inhist_P_posz[i]->SetLineWidth(1);
	inhist_P_posz[i]->SetFillColor(0);

	inhist_P_negz[i] = (TH1F*) infile_P->Get(Form("hIntDistortionP_negz_%d",runs[i]));
	inhist_P_negz[i]->SetLineColor(i+1);
	inhist_P_negz[i]->SetMarkerColor(i+1);
	inhist_P_negz[i]->SetLineWidth(1);
	inhist_P_negz[i]->SetFillColor(0);

	inhist_R_posz[i] = (TH1F*) infile_R->Get(Form("hIntDistortionR_posz_%d",runs[i]));
	inhist_R_posz[i]->SetLineColor(i+1);
	inhist_R_posz[i]->SetMarkerColor(i+1);
	inhist_R_posz[i]->SetLineWidth(1);
	inhist_R_posz[i]->SetFillColor(0);

	inhist_R_negz[i] = (TH1F*) infile_R->Get(Form("hIntDistortionR_negz_%d",runs[i]));
	inhist_R_negz[i]->SetLineColor(i+1);
	inhist_R_negz[i]->SetMarkerColor(i+1);
	inhist_R_negz[i]->SetLineWidth(1);
	inhist_R_negz[i]->SetFillColor(0);

	inhist_Z_posz[i] = (TH1F*) infile_Z->Get(Form("hIntDistortionZ_posz_%d",runs[i]));
	inhist_Z_posz[i]->SetLineColor(i+1);
	inhist_Z_posz[i]->SetMarkerColor(i+1);
	inhist_Z_posz[i]->SetLineWidth(1);
	inhist_Z_posz[i]->SetFillColor(0);

	inhist_Z_negz[i] = (TH1F*) infile_Z->Get(Form("hIntDistortionZ_negz_%d",runs[i]));
	inhist_Z_negz[i]->SetLineColor(i+1);
	inhist_Z_negz[i]->SetMarkerColor(i+1);
	inhist_Z_negz[i]->SetLineWidth(1);
	inhist_Z_negz[i]->SetFillColor(0);
}

std::vector<TH1*> hists_P; hists_P.clear();
std::vector<TH1*> hists_R; hists_R.clear();
std::vector<TH1*> hists_Z; hists_Z.clear();
for (int i=0; i<nruns; i++) { hists_P.push_back(inhist_P_posz[i]); hists_P.push_back(inhist_P_negz[i]); }
for (int i=0; i<nruns; i++) { hists_R.push_back(inhist_R_posz[i]); hists_R.push_back(inhist_R_negz[i]); }
for (int i=0; i<nruns; i++) { hists_Z.push_back(inhist_Z_posz[i]); hists_Z.push_back(inhist_Z_negz[i]); }
std::pair<double,double> yrange_P = SetCommonYRange(hists_P);
std::pair<double,double> yrange_R = SetCommonYRange(hists_R);
std::pair<double,double> yrange_Z = SetCommonYRange(hists_Z);


TCanvas* can = new TCanvas("can","",800,600);
for (int i=0; i<nruns; i++) { inhist_P_posz[i]->Draw("hist,e,same"); }
/*
TCanvas* can = new TCanvas("can","",2400,1200);
can->Divide(3,2);
can->cd(1);
for (int i=0; i<nruns; i++) { inhist_P_posz[i]->Draw("hist,e,same"); }
can->cd(2);
for (int i=0; i<nruns; i++) { inhist_R_posz[i]->Draw("hist,e,same"); }
can->cd(3);
for (int i=0; i<nruns; i++) { inhist_Z_posz[i]->Draw("hist,e,same"); }

can->cd(4);
for (int i=0; i<nruns; i++) { inhist_P_negz[i]->Draw("hist,e,same"); }
can->cd(5);
for (int i=0; i<nruns; i++) { inhist_R_negz[i]->Draw("hist,e,same"); }
can->cd(6);
for (int i=0; i<nruns; i++) { inhist_Z_negz[i]->Draw("hist,e,same"); }
*/

//TLegend *legend = new TLegend(0.01, 0.01, 0.99, 0.99);
TLegend *legend = new TLegend(0.5, 0.5, 0.9, 0.9);
legend->SetHeader(Form("Z = +%d cm",zpoint));
for (int i=0; i<nruns; i++)
{
  legend->AddEntry(inhist_P_posz[i], Form("Run %d",runs[i]), "l");
}

//TCanvas* can_leg = new TCanvas("can_leg","",200,600);
legend->Draw();
can->SaveAs(Form("figure/resid_vsR_from3D_atZ10_runs.pdf"));
//can_leg->SaveAs(Form("figure/resid_vsR_from3D_leg_runs.pdf"));
}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
	//recover default max and min
	h->SetMinimum();
	h->SetMaximum();
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    if (yMin > 0) { yMin *= 0.8; }
    else if (yMin < 0) { yMin *= 1.2; }

    if (yMax > 0) { yMax *= 1.2; }
    else if (yMax < 0) { yMax *= 0.8; }

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
