void make_plot(TH2* h, TString name);
void fit_gauss(TH2* h, TString name, bool verbose=0);
void residual()
{
  gStyle->SetOptStat(0);

  //TFile* file = new TFile("clusters_seeds_53285_qa_all.root");
  //TTree* tree = (TTree*) file->Get("h_QAG4SimulationDistortions_residTree");

  int run = 53744;
  TChain* tree = new TChain("h_QAG4SimulationDistortions_residTree");
  tree->Add(Form("../Reconstructed/%d/clusters_seeds_%d-0-0.root_qa.root",run,run));
  system(Form("mkdir -p figure/%d",run));

  Float_t tanAlpha, tanBeta, drphi, dz, clusR, clusPhi, clusZ, statePhi, stateZ, stateR, stateRPhiErr, stateZErr, clusRPhiErr, clusZErr;
  Int_t event;
  Float_t clusEta, stateEta, statePt, statePz, trackPt;
  Int_t charge;
  Int_t layer;
  std::vector<int> *segtype_tpot = 0;
  std::vector<int> *tileid_tpot = 0;
  std::vector<float> *sclusgx_tpot = 0;
  std::vector<float> *sclusgy_tpot = 0;
  std::vector<float> *sclusgz_tpot = 0;
  std::vector<float> *sclusgr_tpot = 0;
  std::vector<float> *sclusphi_tpot = 0;
  std::vector<float> *scluseta_tpot = 0;

  tree->SetBranchAddress("tanAlpha",&tanAlpha);
  tree->SetBranchAddress("tanBeta",&tanBeta);
  tree->SetBranchAddress("drphi",&drphi);
  tree->SetBranchAddress("dz",&dz);
  tree->SetBranchAddress("clusR",&clusR);
  tree->SetBranchAddress("clusPhi",&clusPhi);
  tree->SetBranchAddress("clusZ",&clusZ);
  tree->SetBranchAddress("statePhi",&statePhi);
  tree->SetBranchAddress("stateZ",&stateZ);
  tree->SetBranchAddress("stateR",&stateR);
  tree->SetBranchAddress("stateRPhiErr",&stateRPhiErr);
  tree->SetBranchAddress("stateZErr",&stateZErr);
  tree->SetBranchAddress("clusRPhiErr",&clusRPhiErr);
  tree->SetBranchAddress("clusZErr",&clusZErr);
  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("clusEta",&clusEta);
  tree->SetBranchAddress("stateEta",&stateEta);
  tree->SetBranchAddress("layer",&layer);
  tree->SetBranchAddress("statePt",&statePt);
  tree->SetBranchAddress("statePz",&statePz);
  tree->SetBranchAddress("trackPt", &trackPt);
  tree->SetBranchAddress("charge", &charge);

  tree->SetBranchAddress("segtype_tpot",&segtype_tpot);
  tree->SetBranchAddress("tileid_tpot",&tileid_tpot);
  tree->SetBranchAddress("sclusgx_tpot",&sclusgx_tpot);
  tree->SetBranchAddress("sclusgy_tpot",&sclusgy_tpot);
  tree->SetBranchAddress("sclusgz_tpot",&sclusgz_tpot);
  tree->SetBranchAddress("sclusgr_tpot",&sclusgr_tpot);
  tree->SetBranchAddress("sclusphi_tpot",&sclusphi_tpot);
  tree->SetBranchAddress("scluseta_tpot",&scluseta_tpot);

  TString histo_prefix = "h_QAG4SimulationDistortions_";

  TH2F *h_beta = new TH2F(histo_prefix + "betadz", "Static distortion off, Run 53285, All charge;tan#beta; #Deltaz [cm]", 100, -1.5, 1.5, 100, -5, 5);
  TH2F *h_beta_pos = new TH2F(histo_prefix + "betadz_pos", "Static distortion off, Run 53285, Positive charge;tan#beta; #Deltaz [cm]", 100, -1.5, 1.5, 100, -5, 5);
  TH2F *h_beta_neg = new TH2F(histo_prefix + "betadz_neg", "Static distortion off, Run 53285, Negative charge;tan#beta; #Deltaz [cm]", 100, -1.5, 1.5, 100, -5, 5);
  TH2F *h_alpha = new TH2F(histo_prefix + "alphardphi", "Static distortion off, Run 53285, All charge;tan#alpha; r#Delta#phi [cm]", 100, -1, 1, 100, -2, 2);
  TH2F *h_alpha_pos = new TH2F(histo_prefix + "alphardphi_pos", "Static distortion off, Run 53285, Positive charge;tan#alpha; r#Delta#phi [cm]", 100, -1, 1, 100, -2, 2);
  TH2F *h_alpha_neg = new TH2F(histo_prefix + "alphardphi_neg", "Static distortion off, Run 53285, Negative charge;tan#alpha; r#Delta#phi [cm]", 100, -1, 1, 100, -2, 2);
  TH2F *h_rphiResid = new TH2F(histo_prefix + "rphiResid", ";r [cm]; #Deltar#phi [cm]", 60, 20, 80, 500, -2, 2);
  TH2F *h_zResid = new TH2F(histo_prefix + "zResid", ";z [cm]; #Deltaz [cm]", 200, -100, 100, 1000, -2, 2);
  TH2F *h_etaResid = new TH2F(histo_prefix + "etaResid", ";#eta;#Delta#eta", 20, -1, 1, 500, -0.2, 0.2);
  TH2F *h_etaResidLayer = new TH2F(histo_prefix + "etaResidLayer", ";r [cm]; #Delta#eta", 60, 20, 80, 500, -0.2, 0.2);
  TH2F *h_zResidLayer = new TH2F(histo_prefix + "zResidLayer", ";r [cm]; #Deltaz [cm]", 60, 20, 80, 1000, -2, 2);
  TH2F *h_deltarphi_layer = new TH2F(histo_prefix + "deltarphi_layer", "Static distortion off, Run 53285, All charge;layer; r#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 500, -3, 3);
  TH2F *h_deltarphi_layer_pos = new TH2F(histo_prefix + "deltarphi_layer_pos", "Static distortion off, Run 53285, Positive charge;layer; r#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 500, -3, 3);
  TH2F *h_deltarphi_layer_neg = new TH2F(histo_prefix + "deltarphi_layer_neg", "Static distortion off, Run 53285, Negative charge;layer; r#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 500, -3, 3);
  TH2F *h_deltaz_layer = new TH2F(histo_prefix + "deltaz_layer", "Static distortion off, Run 53285, All charge;layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 100, -5, 5);
  TH2F *h_deltaz_layer_pos = new TH2F(histo_prefix + "deltaz_layer_pos", "Static distortion off, Run 53285, Positive charge;layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 100, -5, 5);
  TH2F *h_deltaz_layer_neg = new TH2F(histo_prefix + "deltaz_layer_neg", "Static distortion off, Run 53285, Negative charge;layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 100, -5, 5);
  TH2F *h_statez_pulls_layer = new TH2F(histo_prefix + "statez_pulls_layer", ";layer; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 57, 0, 57, 100, -5, 5);
  TH2F *h_staterphi_pulls_layer = new TH2F(histo_prefix + "staterphi_pulls_layer", ";layer; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 57, 0, 57, 100, -5, 5);
  TH2F *h_clusz_pulls_layer = new TH2F(histo_prefix + "clusz_pulls_layer", ";layer; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 57, 0, 57, 100, -5, 5);
  TH2F *h_clusrphi_pulls_layer = new TH2F(histo_prefix + "clusrphi_pulls_layer", ";layer; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 57, 0, 57, 100, -5, 5);

  TH2F *h_deltarphi_z = new TH2F(histo_prefix + "deltarphi_z", "Static distortion off, Run 53285, All charge;z (cm); r#Delta#phi_{track-cluster} (cm)", 100, -80, 80, 500, -3, 3);
  TH2F *h_deltarphi_z_pos = new TH2F(histo_prefix + "deltarphi_z_pos", "Static distortion off, Run 53285, Positive charge;z (cm); r#Delta#phi_{track-cluster} (cm)", 100, -80, 80, 500, -3, 3);
  TH2F *h_deltarphi_z_neg = new TH2F(histo_prefix + "deltarphi_z_neg", "Static distortion off, Run 53285, Negative charge;z (cm); r#Delta#phi_{track-cluster} (cm)", 100, -80, 80, 500, -3, 3);
  TH2F *h_deltaz_z = new TH2F(histo_prefix + "deltaz_z", "Static distortion off, Run 53285, All charge;z (cm); #Deltaz_{track-cluster} (cm)", 100, -80, 80, 100, -5, 5);
  TH2F *h_deltaz_z_pos = new TH2F(histo_prefix + "deltaz_z_pos", "Static distortion off, Run 53285, Positive charge;z (cm); #Deltaz_{track-cluster} (cm)", 100, -80, 80, 100, -5, 5);
  TH2F *h_deltaz_z_neg = new TH2F(histo_prefix + "deltaz_z_neg", "Static distortion off, Run 53285, Negative charge;z (cm); #Deltaz_{track-cluster} (cm)", 100, -80, 80, 100, -5, 5);
  TH2F *h_statez_pulls_z = new TH2F(histo_prefix + "statez_pulls_z", ";z (cm); #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, -80, 80, 100, -5, 5);
  TH2F *h_staterphi_pulls_z = new TH2F(histo_prefix + "staterphi_pulls_z", ";z (cm); #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, -80, 80, 100, -5, 5);
  TH2F *h_clusz_pulls_z = new TH2F(histo_prefix + "clusz_pulls_z", ";z (cm); #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, -80, 80, 100, -5, 5);
  TH2F *h_clusrphi_pulls_z = new TH2F(histo_prefix + "clusrphi_pulls_z", ";z (cm); #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, -80, 80, 100, -5, 5);

  TH2F *h_deltarphi_pt = new TH2F(histo_prefix + "deltarphi_pt", ";pT (GeV/c); r#Delta#phi_{track-cluster} (cm)", 100, 0, 2, 500, -2, 2);
  TH2F *h_deltaz_pt = new TH2F(histo_prefix + "deltaz_pt", ";pT (GeV/c); #Deltaz_{track-cluster} (cm)", 100, 0, 2, 100, -2, 2);
  TH2F *h_statez_pulls_pt = new TH2F(histo_prefix + "statez_pulls_pt", ";pT (GeV/c); #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, 0, 2, 100, -5, 5);
  TH2F *h_staterphi_pulls_pt = new TH2F(histo_prefix + "staterphi_pulls_pt", ";pT (GeV/c); #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, 0, 2, 100, -5, 5);
  TH2F *h_clusz_pulls_pt = new TH2F(histo_prefix + "clusz_pulls_pt", ";pT (GeV/c); #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, 0, 2, 100, -5, 5);
  TH2F *h_clusrphi_pulls_pt = new TH2F(histo_prefix + "clusrphi_pulls_pt", ";pT (GeV/c); #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, 0, 2, 100, -5, 5);

  TH2F *h_deltarphi_pz = new TH2F(histo_prefix + "deltarphi_pz", ";pz (GeV/c); r#Delta#phi_{track-cluster} (cm)", 100, 0, 2, 500, -2, 2);
  TH2F *h_deltaz_pz = new TH2F(histo_prefix + "deltaz_pz", ";pz (GeV/c); #Deltaz_{track-cluster} (cm)", 100, 0, 2, 100, -2, 2);
  TH2F *h_statez_pulls_pz = new TH2F(histo_prefix + "statez_pulls_pz", ";pz (GeV/c); #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, 0, 2, 100, -5, 5);
  TH2F *h_staterphi_pulls_pz = new TH2F(histo_prefix + "staterphi_pulls_pz", ";pz (GeV/c); #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, 0, 2, 100, -5, 5);
  TH2F *h_clusz_pulls_pz = new TH2F(histo_prefix + "clusz_pulls_pz", ";pz (GeV/c); #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, 0, 2, 100, -5, 5);
  TH2F *h_clusrphi_pulls_pz = new TH2F(histo_prefix + "clusrphi_pulls_pz", ";pz (GeV/c); #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, 0, 2, 100, -5, 5);

  TH2F *h_deltarphi_stateeta = new TH2F(histo_prefix + "deltarphi_stateeta", ";track #eta; r#Delta#phi_{track-cluster} (cm)", 100, -1, 1, 500, -2, 2);
  TH2F *h_deltaz_stateeta = new TH2F(histo_prefix + "deltaz_stateeta", ";track #eta; #Deltaz_{track-cluster} (cm)", 100, -1, 1, 100, -2, 2);
  TH2F *h_statez_pulls_stateeta = new TH2F(histo_prefix + "statez_pulls_stateeta", ";track #eta; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, -1, 1, 100, -5, 5);
  TH2F *h_staterphi_pulls_stateeta = new TH2F(histo_prefix + "staterphi_pulls_stateeta", ";track #eta; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, -1, 1, 100, -5, 5);
  TH2F *h_clusz_pulls_stateeta = new TH2F(histo_prefix + "clusz_pulls_stateeta", ";track #eta; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, -1, 1, 100, -5, 5);
  TH2F *h_clusrphi_pulls_stateeta = new TH2F(histo_prefix + "clusrphi_pulls_stateeta", ";track #eta; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, -1, 1, 100, -5, 5);

  TH2F *h_deltarphi_cluseta = new TH2F(histo_prefix + "deltarphi_cluseta", ";Cluster #eta; r#Delta#phi_{track-cluster} (cm)", 100, -1, 1, 500, -2, 2);
  TH2F *h_deltaz_cluseta = new TH2F(histo_prefix + "deltaz_cluseta", ";Cluster #eta; #Deltaz_{track-cluster} (cm)", 100, -1, 1, 100, -2, 2);
  TH2F *h_statez_pulls_cluseta = new TH2F(histo_prefix + "statez_pulls_cluseta", ";Cluster #eta; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, -1, 1, 100, -5, 5);
  TH2F *h_staterphi_pulls_cluseta = new TH2F(histo_prefix + "staterphi_pulls_cluseta", ";Cluster #eta; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, -1, 1, 100, -5, 5);
  TH2F *h_clusz_pulls_cluseta = new TH2F(histo_prefix + "clusz_pulls_cluseta", ";Cluster #eta; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, -1, 1, 100, -5, 5);
  TH2F *h_clusrphi_pulls_cluseta = new TH2F(histo_prefix + "clusrphi_pulls_cluseta", ";Cluster #eta; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, -1, 1, 100, -5, 5);

  TH2F *h_deltarphi_statephi = new TH2F(histo_prefix + "deltarphi_statephi", "Static distortion off, Run 53285, All charge;state #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltarphi_statephi_pos = new TH2F(histo_prefix + "deltarphi_statephi_pos", "Static distortion off, Run 53285, Positive charge;state #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltarphi_statephi_neg = new TH2F(histo_prefix + "deltarphi_statephi_neg", "Static distortion off, Run 53285, Negative charge;state #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltaz_statephi = new TH2F(histo_prefix + "deltaz_statephi", "Static distortion off, Run 53285, All charge;state #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_deltaz_statephi_pos = new TH2F(histo_prefix + "deltaz_statephi", "Static distortion off, Run 53285, Positive charge;state #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_deltaz_statephi_neg = new TH2F(histo_prefix + "deltaz_statephi", "Static distortion off, Run 53285, Negative charge;state #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_statez_pulls_statephi = new TH2F(histo_prefix + "statez_pulls_statephi", ";state #phi [rad]; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_staterphi_pulls_statephi = new TH2F(histo_prefix + "staterphi_pulls_statephi", ";state #phi [rad]; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_clusz_pulls_statephi = new TH2F(histo_prefix + "clusz_pulls_statephi", ";state #phi [rad]; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_clusrphi_pulls_statephi = new TH2F(histo_prefix + "clusrphi_pulls_statephi", ";state #phi [rad]; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, -3.142, 3.142, 100, -5, 5);

  TH2F *h_deltarphi_clusphi = new TH2F(histo_prefix + "deltarphi_clusphi", "Static distortion off, Run 53285, All charge;cluster #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltarphi_clusphi_pos = new TH2F(histo_prefix + "deltarphi_clusphi_pos", "Static distortion off, Run 53285, Positive charge;cluster #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltarphi_clusphi_neg = new TH2F(histo_prefix + "deltarphi_clusphi_neg", "Static distortion off, Run 53285, Negative charge;cluster #phi [rad]; r#Delta#phi_{track-cluster} (cm)", 100, -2.2, -0.8, 500, -3, 3);
  TH2F *h_deltaz_clusphi = new TH2F(histo_prefix + "deltaz_clusphi", "Static distortion off, Run 53285, All charge;cluster #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_deltaz_clusphi_pos = new TH2F(histo_prefix + "deltaz_clusphi_pos", "Static distortion off, Run 53285, Positive charge;cluster #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_deltaz_clusphi_neg = new TH2F(histo_prefix + "deltaz_clusphi_neg", "Static distortion off, Run 53285, Negative charge;cluster #phi [rad]; #Deltaz_{track-cluster} (cm)", 100, -2.2, -0.8, 100, -5, 5);
  TH2F *h_statez_pulls_clusphi = new TH2F(histo_prefix + "statez_pulls_clusphi", ";cluster #phi [rad]; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_staterphi_pulls_clusphi = new TH2F(histo_prefix + "staterphi_pulls_clusphi", ";cluster #phi [rad]; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_clusz_pulls_clusphi = new TH2F(histo_prefix + "clusz_pulls_clusphi", ";cluster #phi [rad]; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 100, -3.142, 3.142, 100, -5, 5);
  TH2F *h_clusrphi_pulls_clusphi = new TH2F(histo_prefix + "clusrphi_pulls_clusphi", ";cluster #phi [rad]; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 100, -3.142, 3.142, 100, -5, 5);


  int nevent = tree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    tree->GetEntry(i);

    int ntpot = tileid_tpot->size();
    bool reject_tpot=false;
    for (int j=0; j<ntpot; j++)
    {
      if (tileid_tpot->at(j)>3) reject_tpot=true;
    }
    if (reject_tpot) continue;

    h_alpha->Fill(tanAlpha, drphi);
    if (charge>0) h_alpha_pos->Fill(tanAlpha, drphi);
    else if (charge<0) h_alpha_neg->Fill(tanAlpha, drphi);
    h_beta->Fill(tanBeta, dz);
    if (charge>0) h_beta_pos->Fill(tanBeta, dz);
    else if (charge<0) h_beta_neg->Fill(tanBeta, dz);

// rphi residual vs R
    h_rphiResid->Fill(clusR, drphi);
// z residual vs Z
    h_zResid->Fill(stateZ, dz);

// (momentum - positional) eta residual vs state Eta
    h_etaResid->Fill(stateEta, clusEta - stateEta);
// z residual vs R
    h_zResidLayer->Fill(clusR, dz);
// (momentum - positional) eta reisual vs R
    h_etaResidLayer->Fill(clusR, clusEta - stateEta);

// layer vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_layer->Fill(layer, drphi);
    if (charge>0) h_deltarphi_layer_pos->Fill(layer, drphi);
    else if (charge<0) h_deltarphi_layer_neg->Fill(layer, drphi);
    h_deltaz_layer->Fill(layer, dz);
    if (charge>0) h_deltaz_layer_pos->Fill(layer, dz);
    else if (charge<0) h_deltaz_layer_neg->Fill(layer, dz);
    h_statez_pulls_layer->Fill(layer, dz / stateZErr);
    h_staterphi_pulls_layer->Fill(layer, drphi / stateRPhiErr);
    h_clusz_pulls_layer->Fill(layer, dz / clusZErr);
    h_clusrphi_pulls_layer->Fill(layer, drphi / clusRPhiErr);

// statez vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_z->Fill(stateZ, drphi);
    if (charge>0) h_deltarphi_z_pos->Fill(stateZ, drphi);
    else if (charge<0) h_deltarphi_z_neg->Fill(stateZ, drphi);
    h_deltaz_z->Fill(stateZ, dz);
    if (charge>0) h_deltaz_z_pos->Fill(stateZ, dz);
    else if (charge<0) h_deltaz_z_neg->Fill(stateZ, dz);
    h_statez_pulls_z->Fill(stateZ, dz / stateZErr);
    h_staterphi_pulls_z->Fill(stateZ, drphi / stateRPhiErr);
    h_clusz_pulls_z->Fill(stateZ, dz / clusZErr);
    h_clusrphi_pulls_z->Fill(stateZ, drphi / clusRPhiErr);

// pt vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_pt->Fill(statePt, drphi);
    h_deltaz_pt->Fill(statePt, dz);
    h_staterphi_pulls_pt->Fill(statePt, drphi / stateRPhiErr);
    h_statez_pulls_pt->Fill(statePt, dz / stateZErr);
    h_clusrphi_pulls_pt->Fill(statePt, drphi / clusRPhiErr);
    h_clusz_pulls_pt->Fill(statePt, dz / clusZErr);

// pz vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_pz->Fill(statePz, drphi);
    h_deltaz_pz->Fill(statePz, dz);
    h_staterphi_pulls_pz->Fill(statePz, drphi / stateRPhiErr);
    h_statez_pulls_pz->Fill(statePz, dz / stateZErr);
    h_clusrphi_pulls_pz->Fill(statePz, drphi / clusRPhiErr);
    h_clusz_pulls_pz->Fill(statePz, dz / clusZErr);

// state eta (momentum) vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_stateeta->Fill(stateEta, drphi);
    h_deltaz_stateeta->Fill(stateEta, dz);
    h_staterphi_pulls_stateeta->Fill(stateEta, drphi / stateRPhiErr);
    h_statez_pulls_stateeta->Fill(stateEta, dz / stateZErr);
    h_clusrphi_pulls_stateeta->Fill(stateEta, drphi / clusRPhiErr);
    h_clusz_pulls_stateeta->Fill(stateEta, dz / clusZErr);

// clus eta (positional) vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_cluseta->Fill(clusEta, drphi);
    h_deltaz_cluseta->Fill(clusEta, dz);
    h_staterphi_pulls_cluseta->Fill(clusEta, drphi / stateRPhiErr);
    h_statez_pulls_cluseta->Fill(clusEta, dz / stateZErr);
    h_clusrphi_pulls_cluseta->Fill(clusEta, drphi / clusRPhiErr);
    h_clusz_pulls_cluseta->Fill(clusEta, dz / clusZErr);

// state phi (positional) vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_statephi->Fill(statePhi, drphi);
    if (charge>0) h_deltarphi_statephi_pos->Fill(statePhi, drphi);
    else if (charge<0) h_deltarphi_statephi_neg->Fill(statePhi, drphi);
    h_deltaz_statephi->Fill(statePhi, dz);
    if (charge>0) h_deltaz_statephi_pos->Fill(statePhi, dz);
    else if (charge<0) h_deltaz_statephi_neg->Fill(statePhi, dz);
    h_staterphi_pulls_statephi->Fill(statePhi, drphi / stateRPhiErr);
    h_statez_pulls_statephi->Fill(statePhi, dz / stateZErr);
    h_clusrphi_pulls_statephi->Fill(statePhi, drphi / clusRPhiErr);
    h_clusz_pulls_statephi->Fill(statePhi, dz / clusZErr);

// clus phi (positional) vs rphi/z residual or residual pulls (clus error or state error)
    h_deltarphi_clusphi->Fill(clusPhi, drphi);
    if (charge>0) h_deltarphi_clusphi_pos->Fill(clusPhi, drphi);
    else if (charge<0) h_deltarphi_clusphi_neg->Fill(clusPhi, drphi);
    h_deltaz_clusphi->Fill(clusPhi, dz);
    if (charge>0) h_deltaz_clusphi_pos->Fill(clusPhi, dz);
    else if (charge<0) h_deltaz_clusphi_neg->Fill(clusPhi, dz);
    h_staterphi_pulls_clusphi->Fill(clusPhi, drphi / stateRPhiErr);
    h_statez_pulls_clusphi->Fill(clusPhi, dz / stateZErr);
    h_clusrphi_pulls_clusphi->Fill(clusPhi, drphi / clusRPhiErr);
    h_clusz_pulls_clusphi->Fill(clusPhi, dz / clusZErr);

  }

  fit_gauss(h_beta, Form("figure/%d/h_beta.pdf",run));
  fit_gauss(h_beta_pos, Form("figure/%d/h_beta_pos.pdf",run));
  fit_gauss(h_beta_neg, Form("figure/%d/h_beta_neg.pdf",run));
  fit_gauss(h_alpha, Form("figure/%d/h_alpha.pdf",run));
  fit_gauss(h_alpha_pos, Form("figure/%d/h_alpha_pos.pdf",run));
  fit_gauss(h_alpha_neg, Form("figure/%d/h_alpha_neg.pdf",run));
  make_plot(h_rphiResid, Form("figure/%d/h_rphiResid.pdf",run));
  make_plot(h_zResid, Form("figure/%d/h_zResid.pdf",run));
  make_plot(h_etaResid, Form("figure/%d/h_etaResid.pdf",run));
  make_plot(h_etaResidLayer, Form("figure/%d/h_etaResidLayer.pdf",run));
  make_plot(h_zResidLayer, Form("figure/%d/h_zResidLayer.pdf",run));
  fit_gauss(h_deltarphi_layer, Form("figure/%d/h_deltarphi_layer.pdf",run));
  fit_gauss(h_deltarphi_layer_pos, Form("figure/%d/h_deltarphi_layer_pos.pdf",run));
  fit_gauss(h_deltarphi_layer_neg, Form("figure/%d/h_deltarphi_layer_neg.pdf",run));
  fit_gauss(h_deltaz_layer, Form("figure/%d/h_deltaz_layer.pdf",run));
  fit_gauss(h_deltaz_layer_pos, Form("figure/%d/h_deltaz_layer_pos.pdf",run));
  fit_gauss(h_deltaz_layer_neg, Form("figure/%d/h_deltaz_layer_neg.pdf",run));
  make_plot(h_statez_pulls_layer, Form("figure/%d/h_statez_pulls_layer.pdf",run));
  make_plot(h_staterphi_pulls_layer, Form("figure/%d/h_staterphi_pulls_layer.pdf",run));
  make_plot(h_clusz_pulls_layer, Form("figure/%d/h_clusz_pulls_layer.pdf",run));
  make_plot(h_clusrphi_pulls_layer, Form("figure/%d/h_clusrphi_pulls_layer.pdf",run));

  fit_gauss(h_deltarphi_z, Form("figure/%d/h_deltarphi_z.pdf",run));
  fit_gauss(h_deltarphi_z_pos, Form("figure/%d/h_deltarphi_z_pos.pdf",run));
  fit_gauss(h_deltarphi_z_neg, Form("figure/%d/h_deltarphi_z_neg.pdf",run));
  fit_gauss(h_deltaz_z, Form("figure/%d/h_deltaz_z.pdf",run));
  fit_gauss(h_deltaz_z_pos, Form("figure/%d/h_deltaz_z_pos.pdf",run));
  fit_gauss(h_deltaz_z_neg, Form("figure/%d/h_deltaz_z_neg.pdf",run));
  make_plot(h_statez_pulls_z, Form("figure/%d/h_statez_pulls_z.pdf",run));
  make_plot(h_staterphi_pulls_z, Form("figure/%d/h_staterphi_pulls_z.pdf",run));
  make_plot(h_clusz_pulls_z, Form("figure/%d/h_clusz_pulls_z.pdf",run));
  make_plot(h_clusrphi_pulls_z, Form("figure/%d/h_clusrphi_pulls_z.pdf",run));

  make_plot(h_deltarphi_pt, Form("figure/%d/h_deltarphi_pt.pdf",run));
  make_plot(h_deltaz_pt, Form("figure/%d/h_deltaz_pt.pdf",run));
  make_plot(h_statez_pulls_pt, Form("figure/%d/h_statez_pulls_pt.pdf",run));
  make_plot(h_staterphi_pulls_pt, Form("figure/%d/h_staterphi_pulls_pt.pdf",run));
  make_plot(h_clusz_pulls_pt, Form("figure/%d/h_clusz_pulls_pt.pdf",run));
  make_plot(h_clusrphi_pulls_pt, Form("figure/%d/h_clusrphi_pulls_pt.pdf",run));

  make_plot(h_deltarphi_pz, Form("figure/%d/h_deltarphi_pz.pdf",run));
  make_plot(h_deltaz_pz, Form("figure/%d/h_deltaz_pz.pdf",run));
  make_plot(h_statez_pulls_pz, Form("figure/%d/h_statez_pulls_pz.pdf",run));
  make_plot(h_staterphi_pulls_pz, Form("figure/%d/h_staterphi_pulls_pz.pdf",run));
  make_plot(h_clusz_pulls_pz, Form("figure/%d/h_clusz_pulls_pz.pdf",run));
  make_plot(h_clusrphi_pulls_pz, Form("figure/%d/h_clusrphi_pulls_pz.pdf",run));

  make_plot(h_deltarphi_stateeta, Form("figure/%d/h_deltarphi_stateeta.pdf",run));
  make_plot(h_deltaz_stateeta, Form("figure/%d/h_deltaz_stateeta.pdf",run));
  make_plot(h_statez_pulls_stateeta, Form("figure/%d/h_statez_pulls_stateeta.pdf",run));
  make_plot(h_staterphi_pulls_stateeta, Form("figure/%d/h_staterphi_pulls_stateeta.pdf",run));
  make_plot(h_clusz_pulls_stateeta, Form("figure/%d/h_clusz_pulls_stateeta.pdf",run));
  make_plot(h_clusrphi_pulls_stateeta, Form("figure/%d/h_clusrphi_pulls_stateeta.pdf",run));

  make_plot(h_deltarphi_cluseta, Form("figure/%d/h_deltarphi_cluseta.pdf",run));
  make_plot(h_deltaz_cluseta, Form("figure/%d/h_deltaz_cluseta.pdf",run));
  make_plot(h_statez_pulls_cluseta, Form("figure/%d/h_statez_pulls_cluseta.pdf",run));
  make_plot(h_staterphi_pulls_cluseta, Form("figure/%d/h_staterphi_pulls_cluseta.pdf",run));
  make_plot(h_clusz_pulls_cluseta, Form("figure/%d/h_clusz_pulls_cluseta.pdf",run));
  make_plot(h_clusrphi_pulls_cluseta, Form("figure/%d/h_clusrphi_pulls_cluseta.pdf",run));

  fit_gauss(h_deltarphi_statephi, Form("figure/%d/h_deltarphi_statephi.pdf",run));
  fit_gauss(h_deltarphi_statephi_pos, Form("figure/%d/h_deltarphi_statephi_pos.pdf",run));
  fit_gauss(h_deltarphi_statephi_neg, Form("figure/%d/h_deltarphi_statephi_neg.pdf",run));
  fit_gauss(h_deltaz_statephi, Form("figure/%d/h_deltaz_statephi.pdf",run));
  fit_gauss(h_deltaz_statephi_pos, Form("figure/%d/h_deltaz_statephi_pos.pdf",run));
  fit_gauss(h_deltaz_statephi_neg, Form("figure/%d/h_deltaz_statephi_neg.pdf",run));
  make_plot(h_statez_pulls_statephi, Form("figure/%d/h_statez_pulls_statephi.pdf",run));
  make_plot(h_staterphi_pulls_statephi, Form("figure/%d/h_staterphi_pulls_statephi.pdf",run));
  make_plot(h_clusz_pulls_statephi, Form("figure/%d/h_clusz_pulls_statephi.pdf",run));
  make_plot(h_clusrphi_pulls_statephi, Form("figure/%d/h_clusrphi_pulls_statephi.pdf",run));

  fit_gauss(h_deltarphi_clusphi, Form("figure/%d/h_deltarphi_clusphi.pdf",run));
  fit_gauss(h_deltarphi_clusphi_pos, Form("figure/%d/h_deltarphi_clusphi_pos.pdf",run));
  fit_gauss(h_deltarphi_clusphi_neg, Form("figure/%d/h_deltarphi_clusphi_neg.pdf",run));
  fit_gauss(h_deltaz_clusphi, Form("figure/%d/h_deltaz_clusphi.pdf",run));
  fit_gauss(h_deltaz_clusphi_pos, Form("figure/%d/h_deltaz_clusphi_pos.pdf",run));
  fit_gauss(h_deltaz_clusphi_neg, Form("figure/%d/h_deltaz_clusphi_neg.pdf",run));
  make_plot(h_statez_pulls_clusphi, Form("figure/%d/h_statez_pulls_clusphi.pdf",run));
  make_plot(h_staterphi_pulls_clusphi, Form("figure/%d/h_staterphi_pulls_clusphi.pdf",run));
  make_plot(h_clusz_pulls_clusphi, Form("figure/%d/h_clusz_pulls_clusphi.pdf",run));
  make_plot(h_clusrphi_pulls_clusphi, Form("figure/%d/h_clusrphi_pulls_clusphi.pdf",run));

}


void make_plot(TH2* h, TString name)
{
  TCanvas* can = new TCanvas("can","",800,600);

  //h->GetXaxis()->SetTitle("#phi_{TPC} - #phi_{Si} [rad]");
  //h->GetYaxis()->SetTitle("Counts");
  //h->SetTitle("#phi_{TPC} - #phi_{Si} for positive");
  //h->SetMaximum(70);
  h->Draw("colz");

  can->Update();
  can->SaveAs(name);
  delete can;
}

void fit_gauss(TH2* h, TString name, bool verbose=0)
{
  TGraphErrors *graph = new TGraphErrors();

  int n = 1;
  for (int i = 1; i <= h->GetNbinsX(); i+=n) {
    TH1D *projection = h->ProjectionY("projection", i, i+n);

    TF1 *gausFit = new TF1("gausFit", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit, "Q");

    if (verbose>0)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      projection->Draw("hist");
      gausFit->SetLineColor(kRed);
      gausFit->Draw("same");
      can->Update();
      can->SaveAs(name + Form(".gausFit_%d.pdf",i));
      delete can;
    }

    Double_t mean = gausFit->GetParameter(1);
    Double_t error = gausFit->GetParError(1);

    double center = 0;
    double width = 0;
    for (int j=0; j<n; j++)
    {
      center += h->GetXaxis()->GetBinCenter(i+j) / n;
      width += h->GetXaxis()->GetBinWidth(i+j) / 2;
    }
    graph->SetPoint(i-1, center, mean);
    graph->SetPointError(i-1, width, error);
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  h->Draw("COLZ");

  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kRed);
  graph->Draw("P");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  c1->Update();
  c1->SaveAs(name);
  delete c1;

}
