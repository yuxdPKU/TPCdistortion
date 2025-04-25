void saveTH1asPDF(TString inputFileName, TString outputDir, TString outputPrefix) {
    TFile* inputFile = new TFile(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Unable to open input file " << inputFileName << std::endl;
        return;
    }

    TList* list = inputFile->GetListOfKeys();
    TIter next(list);
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj && obj->InheritsFrom(TH1::Class())) {
            TH1* hist = (TH1*)obj;
            TString histName = hist->GetName();
            TString outputFileName = Form("%s/%s_%s.pdf", outputDir.Data(), outputPrefix.Data(), histName.Data());

            //cout<<"hist First bin "<<hist->GetBinContent(1)<<" +/- "<<hist->GetBinError(1)<<endl;
	    const auto axis = hist->GetXaxis();
	    //std::cout<<hist->GetName()<<": Xbins "<<axis->GetNbins()<<" , Xmin "<<axis->GetXmin()<<" , Xmax "<<axis->GetXmax()<<std::endl;
	    
            TCanvas* canvas = new TCanvas("canvas", histName, 600, 600);
	    gPad->SetRightMargin(0.10);
            hist->Draw("e");
            canvas->SaveAs(outputFileName);

            delete canvas;
        }
    }

    inputFile->Close();
    delete inputFile;
}

void saveTH2asPDF(TString inputFileName, TString outputDir, TString outputPrefix) {
    TFile* inputFile = new TFile(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Unable to open input file " << inputFileName << std::endl;
        return;
    }

    TList* list = inputFile->GetListOfKeys();
    TIter next(list);
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj && obj->InheritsFrom(TH2::Class())) {
            TH2* hist = (TH2*)obj;
            TString histName = hist->GetName();
            TString outputFileName = Form("%s/%s_%s.pdf", outputDir.Data(), outputPrefix.Data(), histName.Data());

	    int xbins = hist->GetXaxis()->GetNbins();
	    float xmin = hist->GetXaxis()->GetXmin();
	    float xmax = hist->GetXaxis()->GetXmax();
	    int ybins = hist->GetYaxis()->GetNbins();
	    float ymin = hist->GetYaxis()->GetXmin();
	    float ymax = hist->GetYaxis()->GetXmax();
	    TH2* hist_flipped = new TH2F(histName+"_flipped", histName, ybins, ymin, ymax, xbins, xmin, xmax);

            hist_flipped->GetXaxis()->SetTitle(hist->GetYaxis()->GetTitle());
            hist_flipped->GetYaxis()->SetTitle(hist->GetXaxis()->GetTitle());
            for (int i = 1; i <= hist->GetNbinsX(); i++)
	    {
              for (int j = 1; j <= hist->GetNbinsY(); j++)
	      {
		if (hist->GetBinContent(i, j)==0.0)
		{
		  hist_flipped->SetBinContent(j, i, -10000);
		}
		else
		{
                  hist_flipped->SetBinContent(j, i, hist->GetBinContent(i, j));
		}
              }
            }

	    if (histName=="hentries_negz" || histName=="hentries_posz")
	    {
	      hist_flipped->SetMinimum(0);
	    }
	    if (histName=="hIntDistortionP_negz" || histName=="hIntDistortionP_posz")
	    {
	      hist_flipped->SetMaximum(0.02);
	      hist_flipped->SetMinimum(-0.02);
	    }
	    if (histName=="hIntDistortionRP_negz" || histName=="hIntDistortionRP_posz")
	    {
	      hist_flipped->SetMaximum(1.5);
	      hist_flipped->SetMinimum(-1.5);
	    }
	    if (histName=="hIntDistortionR_negz" || histName=="hIntDistortionR_posz")
	    {
	      hist_flipped->SetMaximum(2);
	      hist_flipped->SetMinimum(-2);
	    }
	    if (histName=="hIntDistortionZ_negz" || histName=="hIntDistortionZ_posz")
	    {
	      hist_flipped->SetMaximum(2);
	      hist_flipped->SetMinimum(-2);
	    }

            TCanvas* canvas = new TCanvas("canvas", histName, 600, 600);
	    gPad->SetRightMargin(0.15);
            hist_flipped->Draw("colz");
            canvas->SaveAs(outputFileName);

            delete canvas;
        }
    }

    inputFile->Close();
    delete inputFile;
}

void drawSingle() {
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(100);
    TGaxis::SetMaxDigits(3);
    const int nrun=6;
    //int runs[nrun] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};
    int runs[nrun] = {53534, 53630, 53744, 53756, 53877, 53876};

    /*
    for(int i=0; i<nrun; i++)
    {
      TString inputFileName = Form("Rootfiles/Distortions_1D_mm_%d_layer.root", runs[i]);
      TString outputDir = Form("figure_%d",runs[i]);
      gSystem->mkdir(outputDir, kTRUE);
      saveTH1asPDF(inputFileName, outputDir, Form("1D_layer"));
    }
    */

    /*
    for(int i=0; i<nrun; i++)
    {
      TString inputFileName = Form("Rootfiles/Distortions_1D_mm_%d_radius.root", runs[i]);
      TString outputDir = Form("figure_%d",runs[i]);
      gSystem->mkdir(outputDir, kTRUE);
      saveTH1asPDF(inputFileName, outputDir, Form("1D_radius"));
    }
    */

    gStyle->SetHistMinimumZero(kTRUE);
    for(int i=0; i<nrun; i++)
    {
      TString inputFileName = Form("Rootfiles/Distortions_2D_mm_%d_rz.root", runs[i]);
      TString outputDir = Form("figure_%d",runs[i]);
      gSystem->mkdir(outputDir, kTRUE);
      saveTH2asPDF(inputFileName, outputDir, Form("2D_rz"));
    }
}
