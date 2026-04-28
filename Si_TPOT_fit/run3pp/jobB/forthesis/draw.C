void draw()
{

int runs[9] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};

TFile* infile = new TFile("../Rootfiles/hist_rdphi_1Dr_Z10.root","");
TH1F* inhist[9];

for (int i=0; i<9; i++)
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
for (int i=0; i<9; i++) { inhist[i]->Draw("hist,e,same"); }
can->SaveAs(Form("resid_vsR_from3D_runs.pdf"));


}
