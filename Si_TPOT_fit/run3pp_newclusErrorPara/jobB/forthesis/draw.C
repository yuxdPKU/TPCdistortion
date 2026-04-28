void draw()
{

const int nruns = 2;
int runs[nruns] = {79507,79516};

TFile* infile = new TFile("../Rootfiles/hist_rdphi_1Dr_Z10.root","");
TH1F* inhist[nruns];

for (int i=0; i<nruns; i++)
{
	inhist[i] = (TH1F*) infile->Get(Form("hIntDistortionP_posz_%d",runs[i]));
	inhist[i]->SetLineColor(i+1);
	inhist[i]->SetMarkerColor(i+1);
	inhist[i]->SetLineWidth(1);
	inhist[i]->SetFillColor(0);
}

TCanvas* can = new TCanvas("can","",800,600);
can->cd(1);
gPad->SetLogy(0);
for (int i=0; i<nruns; i++) { inhist[i]->Draw("hist,e,same"); }

TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
legend->SetHeader("Z = +10 cm");
for (int i=0; i<nruns; i++)
{
  legend->AddEntry(inhist[i], Form("Run %d",runs[i]), "l");
}
legend->Draw();

can->SaveAs(Form("resid_vsR_from3D_runs.pdf"));


}
