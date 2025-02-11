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

            //TString xAxisTitle = hist->GetXaxis()->GetTitle();
            //TString yAxisTitle = hist->GetYaxis()->GetTitle();
            //hist->GetXaxis()->SetTitle(yAxisTitle);
            //hist->GetYaxis()->SetTitle(xAxisTitle);

cout<<"hist First bin "<<hist->GetBinContent(1)<<" +/- "<<hist->GetBinError(1)<<endl;
	    const auto axis = hist->GetXaxis();
	    std::cout<<hist->GetName()<<": Xbins "<<axis->GetNbins()<<" , Xmin "<<axis->GetXmin()<<" , Xmax "<<axis->GetXmax()<<std::endl;
	    
            TCanvas* canvas = new TCanvas("canvas", histName, 800, 600);
            hist->Draw();
            canvas->SaveAs(outputFileName);

            delete canvas;
        }
    }

    inputFile->Close();
    delete inputFile;
}

void drawSingle() {
    int runs[7] = {53285, 53534, 53744, 53756, 53877, 53876, 53630};

    for(int i=0; i<6; i++)
    {
      TString inputFileName = Form("Rootfiles/Distortions_1D_mm_%d_layer.root", runs[i]);
      TString outputDir = Form("figure_%d",runs[i]);
      gSystem->mkdir(outputDir, kTRUE);
      saveTH1asPDF(inputFileName, outputDir, Form("1D_layer"));
    }

    for(int i=0; i<6; i++)
    {
      TString inputFileName = Form("Rootfiles/Distortions_1D_mm_%d_radius.root", runs[i]);
      TString outputDir = Form("figure_%d",runs[i]);
      gSystem->mkdir(outputDir, kTRUE);
      saveTH1asPDF(inputFileName, outputDir, Form("1D_radius"));
    }

}
