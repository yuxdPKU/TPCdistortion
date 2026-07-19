#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
  struct StateBranches
  {
    std::vector<unsigned long>* cluskey = nullptr;
    std::vector<unsigned int>* layer = nullptr;
    std::vector<float>* pathlength = nullptr;
    std::vector<float>* x = nullptr;
    std::vector<float>* y = nullptr;
    std::vector<float>* z = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<float>* xError = nullptr;
    std::vector<float>* yError = nullptr;
    std::vector<float>* zError = nullptr;
    std::vector<float>* clusterX = nullptr;
    std::vector<float>* clusterY = nullptr;
    std::vector<float>* clusterZ = nullptr;

    bool bind(TTree* tree, const std::string& prefix)
    {
      return tree->SetBranchAddress((prefix + "_state_cluskey").c_str(), &cluskey) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_layer").c_str(), &layer) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_pathlength").c_str(), &pathlength) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_x").c_str(), &x) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_y").c_str(), &y) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_z").c_str(), &z) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_px").c_str(), &px) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_py").c_str(), &py) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_pz").c_str(), &pz) >= 0;
    }

    bool bindPositionErrors(TTree* tree, const std::string& prefix)
    {
      return tree->SetBranchAddress((prefix + "_state_x_error").c_str(), &xError) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_y_error").c_str(), &yError) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_z_error").c_str(), &zError) >= 0;
    }

    bool bindClusterPositions(TTree* tree, const std::string& prefix)
    {
      return tree->SetBranchAddress((prefix + "_state_cluster_x").c_str(), &clusterX) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_cluster_y").c_str(), &clusterY) >= 0 &&
             tree->SetBranchAddress((prefix + "_state_cluster_z").c_str(), &clusterZ) >= 0;
    }

    bool hasConsistentSizes() const
    {
      if (!cluskey || !layer || !pathlength || !x || !y || !z || !px || !py || !pz)
      {
        return false;
      }

      const auto size = cluskey->size();
      return layer->size() == size && pathlength->size() == size &&
             x->size() == size && y->size() == size && z->size() == size &&
             px->size() == size && py->size() == size && pz->size() == size;
    }

    std::size_t size() const { return cluskey ? cluskey->size() : 0; }

    bool hasConsistentPositionErrorSizes() const
    {
      return xError && yError && zError &&
             xError->size() == size() && yError->size() == size() && zError->size() == size();
    }

    bool hasConsistentClusterPositionSizes() const
    {
      return clusterX && clusterY && clusterZ &&
             clusterX->size() == size() && clusterY->size() == size() && clusterZ->size() == size();
    }
  };
}  // namespace

void TrackStateAnalysis(
    const char* inputFileName = "all_acts.root",
    const char* outputFileName = "hist_acts.root",
    Long64_t maxEntries = -1,
    const char* figureDirectory = "figure")
{
  constexpr unsigned int tpcFirstLayer = 7;
  constexpr unsigned int tpcLastLayer = 54;
  constexpr int tpcLayerCount = 48;
  constexpr double tpcMinRadius = 31.105;  // cm
  constexpr double tpcMaxRadius = 75.911;  // cm

  TFile inputFile(inputFileName, "READ");
  if (inputFile.IsZombie())
  {
    std::cerr << "Could not open " << inputFileName << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(inputFile.Get("trackStateComparison"));
  if (!tree)
  {
    std::cerr << "Could not find TTree trackStateComparison in " << inputFileName << std::endl;
    return;
  }

  // Track-pair level branches: one entry is one truth/reco track pair.
  UInt_t event = 0;
  UInt_t truthTrackId = 0;
  UInt_t recoTrackId = 0;
  Float_t truthTrackPt = std::numeric_limits<float>::quiet_NaN();
  Float_t recoTrackPt = std::numeric_limits<float>::quiet_NaN();
  ULong_t siliconSeedId = 0;
  ULong_t tpcSeedId = 0;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("truth_track_id", &truthTrackId);
  tree->SetBranchAddress("reco_track_id", &recoTrackId);
  if (!tree->GetBranch("truth_track_pt") || !tree->GetBranch("reco_track_pt"))
  {
    std::cerr << "Missing required truth_track_pt/reco_track_pt branches in "
              << inputFileName << std::endl;
    return;
  }
  tree->SetBranchAddress("truth_track_pt", &truthTrackPt);
  tree->SetBranchAddress("reco_track_pt", &recoTrackPt);
  tree->SetBranchAddress("silicon_seed_id", &siliconSeedId);
  tree->SetBranchAddress("tpc_seed_id", &tpcSeedId);

  StateBranches truth;
  StateBranches reco;
  if (!truth.bind(tree, "truth") || !reco.bind(tree, "reco"))
  {
    std::cerr << "Failed to bind one or more state branches" << std::endl;
    return;
  }
  if (!reco.bindPositionErrors(tree, "reco"))
  {
    std::cerr << "Failed to bind reco state position-error branches" << std::endl;
    return;
  }
  if (!reco.bindClusterPositions(tree, "reco"))
  {
    std::cerr << "Failed to bind reco state cluster-position branches" << std::endl;
    return;
  }

  TFile outputFile(outputFileName, "RECREATE");

  // Example histograms. Change binning/ranges or add your own histograms here.
  auto* hDx = new TH1F("h_dx", ";x_{reco}-x_{truth} [cm];States", 200, -1., 1.);
  auto* hDy = new TH1F("h_dy", ";y_{reco}-y_{truth} [cm];States", 200, -1., 1.);
  auto* hDz = new TH1F("h_dz", ";z_{reco}-z_{truth} [cm];States", 200, -1., 1.);
  auto* hDpx = new TH1F("h_dpx", ";p_{x,reco}-p_{x,truth} [GeV];States", 200, -0.5, 0.5);
  auto* hDpy = new TH1F("h_dpy", ";p_{y,reco}-p_{y,truth} [GeV];States", 200, -0.5, 0.5);
  auto* hDpz = new TH1F("h_dpz", ";p_{z,reco}-p_{z,truth} [GeV];States", 200, -0.5, 0.5);
  auto* hDphi = new TH1F(
      "h_dphi", ";#phi_{reco}-#phi_{truth} [rad];States", 200, -0.02, 0.02);
  auto* hDr = new TH1F(
      "h_dr", ";r_{reco}-r_{truth} [cm];States", 200, -0.02, 0.02);
  auto* hTruthStateRMinusRecoClusterR = new TH1F(
      "h_truth_state_r_minus_reco_cluster_r",
      ";r_{truth state}-r_{reco cluster} [cm];Matched TPC states",
      400, -0.1, 0.1);

  // Matched truth/reco state positions and reconstructed-state uncertainties.
  auto* hTruthStateX = new TH1F(
      "h_truth_state_x", ";Truth state x [cm];Matched TPC states", 200, -80., 80.);
  auto* hTruthStateY = new TH1F(
      "h_truth_state_y", ";Truth state y [cm];Matched TPC states", 200, -80., 80.);
  auto* hTruthStateZ = new TH1F(
      "h_truth_state_z", ";Truth state z [cm];Matched TPC states", 200, -100., 100.);
  auto* hRecoStateX = new TH1F(
      "h_reco_state_x", ";Reco state x [cm];Matched TPC states", 200, -80., 80.);
  auto* hRecoStateY = new TH1F(
      "h_reco_state_y", ";Reco state y [cm];Matched TPC states", 200, -80., 80.);
  auto* hRecoStateZ = new TH1F(
      "h_reco_state_z", ";Reco state z [cm];Matched TPC states", 200, -100., 100.);
  auto* hRecoStateXError = new TH1F(
      "h_reco_state_x_error", ";Reco state #sigma_{x} [cm];Matched TPC states", 200, 0., 1.);
  auto* hRecoStateYError = new TH1F(
      "h_reco_state_y_error", ";Reco state #sigma_{y} [cm];Matched TPC states", 200, 0., 1.);
  auto* hRecoStateZError = new TH1F(
      "h_reco_state_z_error", ";Reco state #sigma_{z} [cm];Matched TPC states", 200, 0., 1.);

  auto* hRecoStateXErrorVsTruthR = new TH2F(
      "h_reco_state_x_error_vs_truth_r", ";r_{truth} [cm];Reco state #sigma_{x} [cm]",
      48, tpcMinRadius, tpcMaxRadius, 200, 0., 1.);
  auto* hRecoStateXErrorVsTruthPhi = new TH2F(
      "h_reco_state_x_error_vs_truth_phi", ";#phi_{truth} [rad];Reco state #sigma_{x} [cm]",
      72, -TMath::Pi(), TMath::Pi(), 200, 0., 1.);
  auto* hRecoStateXErrorVsTruthZ = new TH2F(
      "h_reco_state_x_error_vs_truth_z", ";z_{truth} [cm];Reco state #sigma_{x} [cm]",
      100, -100., 100., 200, 0., 1.);
  auto* hRecoStateYErrorVsTruthR = new TH2F(
      "h_reco_state_y_error_vs_truth_r", ";r_{truth} [cm];Reco state #sigma_{y} [cm]",
      48, tpcMinRadius, tpcMaxRadius, 200, 0., 1.);
  auto* hRecoStateYErrorVsTruthPhi = new TH2F(
      "h_reco_state_y_error_vs_truth_phi", ";#phi_{truth} [rad];Reco state #sigma_{y} [cm]",
      72, -TMath::Pi(), TMath::Pi(), 200, 0., 1.);
  auto* hRecoStateYErrorVsTruthZ = new TH2F(
      "h_reco_state_y_error_vs_truth_z", ";z_{truth} [cm];Reco state #sigma_{y} [cm]",
      100, -100., 100., 200, 0., 1.);
  auto* hRecoStateZErrorVsTruthR = new TH2F(
      "h_reco_state_z_error_vs_truth_r", ";r_{truth} [cm];Reco state #sigma_{z} [cm]",
      48, tpcMinRadius, tpcMaxRadius, 200, 0., 1.);
  auto* hRecoStateZErrorVsTruthPhi = new TH2F(
      "h_reco_state_z_error_vs_truth_phi", ";#phi_{truth} [rad];Reco state #sigma_{z} [cm]",
      72, -TMath::Pi(), TMath::Pi(), 200, 0., 1.);
  auto* hRecoStateZErrorVsTruthZ = new TH2F(
      "h_reco_state_z_error_vs_truth_z", ";z_{truth} [cm];Reco state #sigma_{z} [cm]",
      100, -100., 100., 200, 0., 1.);

  // Track-pair transverse momentum distributions.
  auto* hTruthTrackPt = new TH1F(
      "h_truth_track_pt", ";Track p_{T} [GeV/c];Track pairs", 200, 0., 2.);
  auto* hRecoTrackPt = new TH1F(
      "h_reco_track_pt", ";Track p_{T} [GeV/c];Track pairs", 200, 0., 2.);
  auto* hDeltaTrackPt = new TH1F(
      "h_delta_track_pt", ";p_{T,reco}-p_{T,truth} [GeV/c];Track pairs", 200, -0.2, 0.2);

  // State residuals versus the truth track pT. Each matched TPC state is one
  // entry, while truthTrackPt is constant for all states in a track pair.
  auto* hDxVsTruthPt = new TH2F(
      "h_dx_vs_truth_pt", ";p_{T,truth} [GeV/c];#Deltax [cm]",
      200, 0., 2., 200, -1., 1.);
  auto* hDyVsTruthPt = new TH2F(
      "h_dy_vs_truth_pt", ";p_{T,truth} [GeV/c];#Deltay [cm]",
      200, 0., 2., 200, -1., 1.);
  auto* hDzVsTruthPt = new TH2F(
      "h_dz_vs_truth_pt", ";p_{T,truth} [GeV/c];#Deltaz [cm]",
      200, 0., 2., 200, -1., 1.);
  auto* hDphiVsTruthPt = new TH2F(
      "h_dphi_vs_truth_pt", ";p_{T,truth} [GeV/c];#Delta#phi [rad]",
      200, 0., 2., 200, -0.02, 0.02);

  // Residual distributions in three truth-track pT intervals:
  // [0.2, 0.6), [0.6, 1.0], and (1.0, infinity) GeV/c.
  std::array<TH1F*, 3> hDxTruthPtBins = {
      new TH1F("h_dx_truth_pt_0p2_0p6", "0.2 #leq p_{T,truth} < 0.6 GeV/c;#Deltax [cm];States", 200, -1., 1.),
      new TH1F("h_dx_truth_pt_0p6_1p0", "0.6 #leq p_{T,truth} #leq 1.0 GeV/c;#Deltax [cm];States", 200, -1., 1.),
      new TH1F("h_dx_truth_pt_gt_1p0", "p_{T,truth} > 1.0 GeV/c;#Deltax [cm];States", 200, -1., 1.)};
  std::array<TH1F*, 3> hDyTruthPtBins = {
      new TH1F("h_dy_truth_pt_0p2_0p6", "0.2 #leq p_{T,truth} < 0.6 GeV/c;#Deltay [cm];States", 200, -1., 1.),
      new TH1F("h_dy_truth_pt_0p6_1p0", "0.6 #leq p_{T,truth} #leq 1.0 GeV/c;#Deltay [cm];States", 200, -1., 1.),
      new TH1F("h_dy_truth_pt_gt_1p0", "p_{T,truth} > 1.0 GeV/c;#Deltay [cm];States", 200, -1., 1.)};
  std::array<TH1F*, 3> hDzTruthPtBins = {
      new TH1F("h_dz_truth_pt_0p2_0p6", "0.2 #leq p_{T,truth} < 0.6 GeV/c;#Deltaz [cm];States", 200, -1., 1.),
      new TH1F("h_dz_truth_pt_0p6_1p0", "0.6 #leq p_{T,truth} #leq 1.0 GeV/c;#Deltaz [cm];States", 200, -1., 1.),
      new TH1F("h_dz_truth_pt_gt_1p0", "p_{T,truth} > 1.0 GeV/c;#Deltaz [cm];States", 200, -1., 1.)};
  std::array<TH1F*, 3> hDphiTruthPtBins = {
      new TH1F("h_dphi_truth_pt_0p2_0p6", "0.2 #leq p_{T,truth} < 0.6 GeV/c;#Delta#phi [rad];States", 200, -0.02, 0.02),
      new TH1F("h_dphi_truth_pt_0p6_1p0", "0.6 #leq p_{T,truth} #leq 1.0 GeV/c;#Delta#phi [rad];States", 200, -0.02, 0.02),
      new TH1F("h_dphi_truth_pt_gt_1p0", "p_{T,truth} > 1.0 GeV/c;#Delta#phi [rad];States", 200, -0.02, 0.02)};

  // Residual versus the 48 TPC layers (7--54).
  auto* hDxVsLayer = new TH2F(
      "h_dx_vs_layer", ";Layer;x_{reco}-x_{truth} [cm]",
      tpcLayerCount, tpcFirstLayer - 0.5, tpcLastLayer + 0.5, 200, -1., 1.);
  auto* hDyVsLayer = new TH2F(
      "h_dy_vs_layer", ";Layer;y_{reco}-y_{truth} [cm]",
      tpcLayerCount, tpcFirstLayer - 0.5, tpcLastLayer + 0.5, 200, -1., 1.);
  auto* hDzVsLayer = new TH2F(
      "h_dz_vs_layer", ";Layer;z_{reco}-z_{truth} [cm]",
      tpcLayerCount, tpcFirstLayer - 0.5, tpcLastLayer + 0.5, 200, -1., 1.);

  // Residuals versus truth-state position. Truth coordinates are used on the
  // x-axis so that the coordinate itself is not shifted by the reco residual.
  auto* hDphiVsTruthR = new TH2F(
      "h_dphi_vs_truth_r", ";r_{truth} [cm];#Delta#phi [rad]",
      48, tpcMinRadius, tpcMaxRadius, 200, -0.005, 0.005);
  auto* hDphiVsTruthZ = new TH2F(
      "h_dphi_vs_truth_z", ";z_{truth} [cm];#Delta#phi [rad]",
      220, -100., 100., 200, -0.005, 0.005);
  auto* hDphiVsTruthPhi = new TH2F(
      "h_dphi_vs_truth_phi", ";#phi_{truth} [rad];#Delta#phi [rad]",
      72, -TMath::Pi(), TMath::Pi(), 200, -0.005, 0.005);
  auto* hDzVsTruthR = new TH2F(
      "h_dz_vs_truth_r", ";r_{truth} [cm];#Deltaz [cm]",
      48, tpcMinRadius, tpcMaxRadius, 200, -1., 1.);
  auto* hDzVsTruthZ = new TH2F(
      "h_dz_vs_truth_z", ";z_{truth} [cm];#Deltaz [cm]",
      220, -100., 100., 200, -1., 1.);
  auto* hDzVsTruthPhi = new TH2F(
      "h_dz_vs_truth_phi", ";#phi_{truth} [rad];#Deltaz [cm]",
      72, -TMath::Pi(), TMath::Pi(), 200, -1., 1.);

  // One-dimensional profiles: each bin stores the mean residual after
  // integrating over the other truth-position coordinates.
  auto* pDphiVsTruthR = new TProfile(
      "p_dphi_vs_truth_r", ";r_{truth} [cm];<#Delta#phi> [rad]",
      48, tpcMinRadius, tpcMaxRadius);
  auto* pDphiVsTruthZ = new TProfile(
      "p_dphi_vs_truth_z", ";z_{truth} [cm];<#Delta#phi> [rad]",
      220, -100., 100.);
  auto* pDphiVsTruthPhi = new TProfile(
      "p_dphi_vs_truth_phi", ";#phi_{truth} [rad];<#Delta#phi> [rad]",
      72, -TMath::Pi(), TMath::Pi());
  auto* pDzVsTruthR = new TProfile(
      "p_dz_vs_truth_r", ";r_{truth} [cm];<#Deltaz> [cm]",
      48, tpcMinRadius, tpcMaxRadius);
  auto* pDzVsTruthZ = new TProfile(
      "p_dz_vs_truth_z", ";z_{truth} [cm];<#Deltaz> [cm]",
      220, -100., 100.);
  auto* pDzVsTruthPhi = new TProfile(
      "p_dz_vs_truth_phi", ";#phi_{truth} [rad];<#Deltaz> [cm]",
      72, -TMath::Pi(), TMath::Pi());

  // Two-dimensional profiles: color is the mean residual in each small
  // truth-coordinate bin.
  auto* pDphiTruthRPhi = new TProfile2D(
      "p_dphi_truth_r_phi",
      ";r_{truth} [cm];#phi_{truth} [rad];<#Delta#phi> [rad]",
      48, tpcMinRadius, tpcMaxRadius, 72, -TMath::Pi(), TMath::Pi());
  auto* pDzTruthRPhi = new TProfile2D(
      "p_dz_truth_r_phi",
      ";r_{truth} [cm];#phi_{truth} [rad];<#Deltaz> [cm]",
      48, tpcMinRadius, tpcMaxRadius, 72, -TMath::Pi(), TMath::Pi());
  auto* pDphiTruthZPhi = new TProfile2D(
      "p_dphi_truth_z_phi",
      ";z_{truth} [cm];#phi_{truth} [rad];<#Delta#phi> [rad]",
      100, -100., 100., 72, -TMath::Pi(), TMath::Pi());
  auto* pDzTruthZPhi = new TProfile2D(
      "p_dz_truth_z_phi",
      ";z_{truth} [cm];#phi_{truth} [rad];<#Deltaz> [cm]",
      100, -100., 100., 72, -TMath::Pi(), TMath::Pi());
  auto* pDphiTruthZR = new TProfile2D(
      "p_dphi_truth_z_r",
      ";z_{truth} [cm];r_{truth} [cm];<#Delta#phi> [rad]",
      100, -100., 100., 48, tpcMinRadius, tpcMaxRadius);
  auto* pDzTruthZR = new TProfile2D(
      "p_dz_truth_z_r",
      ";z_{truth} [cm];r_{truth} [cm];<#Deltaz> [cm]",
      100, -100., 100., 48, tpcMinRadius, tpcMaxRadius);

  const Long64_t entriesToProcess =
      maxEntries < 0 ? tree->GetEntries() : std::min(maxEntries, tree->GetEntries());

  Long64_t badVectorEntries = 0;
  Long64_t matchedStates = 0;
  Long64_t dphiProfileStates = 0;
  Long64_t dzProfileStates = 0;
  Long64_t invalidTrackPtEntries = 0;
  Long64_t invalidProjectionStates = 0;

  for (Long64_t entry = 0; entry < entriesToProcess; ++entry)
  {
    if (tree->GetEntry(entry) <= 0)
    {
      continue;
    }

    if (!truth.hasConsistentSizes() || !reco.hasConsistentSizes() ||
        !reco.hasConsistentPositionErrorSizes() ||
        !reco.hasConsistentClusterPositionSizes())
    {
      ++badVectorEntries;
      continue;
    }

    const bool validTruthTrackPt = std::isfinite(truthTrackPt) && truthTrackPt >= 0.F;
    const bool validRecoTrackPt = std::isfinite(recoTrackPt) && recoTrackPt >= 0.F;
    if (validTruthTrackPt)
    {
      hTruthTrackPt->Fill(truthTrackPt);
    }
    if (validRecoTrackPt)
    {
      hRecoTrackPt->Fill(recoTrackPt);
    }
    if (validTruthTrackPt && validRecoTrackPt)
    {
      hDeltaTrackPt->Fill(recoTrackPt - truthTrackPt);
    }
    else
    {
      ++invalidTrackPtEntries;
    }

    int truthPtBin = -1;
    if (validTruthTrackPt && truthTrackPt >= 0.2F && truthTrackPt < 0.6F)
    {
      truthPtBin = 0;
    }
    else if (validTruthTrackPt && truthTrackPt >= 0.6F && truthTrackPt <= 1.0F)
    {
      truthPtBin = 1;
    }
    else if (validTruthTrackPt && truthTrackPt > 1.0F)
    {
      truthPtBin = 2;
    }

    // ================================================================
    // Track-pair level analysis goes here.
    // Available variables:
    //   event, truthTrackId, recoTrackId, truthTrackPt, recoTrackPt,
    //   siliconSeedId, tpcSeedId
    // Example:
    //   hTrackMultiplicity->Fill(...);
    // ================================================================

    // Loop over all truth states, whether or not a reco state exists.
    for (std::size_t truthIndex = 0; truthIndex < truth.size(); ++truthIndex)
    {
      const auto truthKey = truth.cluskey->at(truthIndex);
      const auto truthLayer = truth.layer->at(truthIndex);

      // Fill truth-only histograms here using, for example:
      // truth.x->at(truthIndex), truth.px->at(truthIndex), truthLayer.
      (void) truthKey;
      (void) truthLayer;
    }

    // Loop over all reco states, whether or not a truth state exists.
    for (std::size_t recoIndex = 0; recoIndex < reco.size(); ++recoIndex)
    {
      const auto recoKey = reco.cluskey->at(recoIndex);
      const auto recoLayer = reco.layer->at(recoIndex);

      // Fill reco-only histograms here.
      // Position uncertainties are available as reco.xError/yError/zError.
      (void) recoKey;
      (void) recoLayer;
    }

    // Match truth and reco states by cluster key. Do not match by vector index:
    // Acts and GenFit can store their states in different orders.
    std::unordered_map<unsigned long, std::size_t> recoIndexByClusterKey;
    recoIndexByClusterKey.reserve(reco.size());

    for (std::size_t recoIndex = 0; recoIndex < reco.size(); ++recoIndex)
    {
      const auto key = reco.cluskey->at(recoIndex);
      if (key == std::numeric_limits<unsigned long>::max())
      {
        continue;  // path-length-zero state without an associated cluster
      }
      recoIndexByClusterKey.emplace(key, recoIndex);
    }

    for (std::size_t truthIndex = 0; truthIndex < truth.size(); ++truthIndex)
    {
      const auto key = truth.cluskey->at(truthIndex);
      if (key == std::numeric_limits<unsigned long>::max())
      {
        continue;
      }

      const auto recoIter = recoIndexByClusterKey.find(key);
      if (recoIter == recoIndexByClusterKey.end())
      {
        continue;
      }

      const std::size_t recoIndex = recoIter->second;
      const unsigned int layer = truth.layer->at(truthIndex);

      // This analysis is restricted to the 48 TPC layers. In particular,
      // exclude silicon (0--6) and TPOT/Micromegas (55--56) states.
      if (layer < tpcFirstLayer || layer > tpcLastLayer)
      {
        continue;
      }

      const double truthX = truth.x->at(truthIndex);
      const double truthY = truth.y->at(truthIndex);
      const double truthZ = truth.z->at(truthIndex);
      const double truthR = std::hypot(truthX, truthY);
      const double truthPhi = std::atan2(truthY, truthX);

      // Following PHTpcResiduals, linearly propagate the reconstructed state
      // along its local momentum direction to the radius of the associated
      // reconstructed cluster. This removes the coupling caused by comparing
      // truth and reco states evaluated at different radii.
      const double recoStateX = reco.x->at(recoIndex);
      const double recoStateY = reco.y->at(recoIndex);
      const double recoStateZ = reco.z->at(recoIndex);
      const double recoStatePx = reco.px->at(recoIndex);
      const double recoStatePy = reco.py->at(recoIndex);
      const double recoStatePz = reco.pz->at(recoIndex);
      const double clusterX = reco.clusterX->at(recoIndex);
      const double clusterY = reco.clusterY->at(recoIndex);
      const double trackR = std::hypot(recoStateX, recoStateY);
      const double clusterR = std::hypot(clusterX, clusterY);

      // Cross-check only: the upstream truth fitter is expected to have
      // placed the truth state at the cluster radius already. Do not
      // extrapolate or otherwise modify the truth state in this analysis.
      if (std::isfinite(truthR) && std::isfinite(clusterR))
      {
        hTruthStateRMinusRecoClusterR->Fill(truthR - clusterR);
      }

      if (!std::isfinite(trackR) || !std::isfinite(clusterR) ||
          !std::isfinite(recoStateZ) || !std::isfinite(recoStatePx) ||
          !std::isfinite(recoStatePy) || !std::isfinite(recoStatePz) ||
          trackR <= 0.)
      {
        ++invalidProjectionStates;
        continue;
      }

      const double radiusStep = clusterR - trackR;
      const double trackDrDt =
          (recoStateX * recoStatePx + recoStateY * recoStatePy) / trackR;
      if (!std::isfinite(trackDrDt) || std::abs(trackDrDt) < 1.e-8)
      {
        ++invalidProjectionStates;
        continue;
      }

      const double trackDxDr = recoStatePx / trackDrDt;
      const double trackDyDr = recoStatePy / trackDrDt;
      const double trackDzDr = recoStatePz / trackDrDt;
      const double recoX = recoStateX + radiusStep * trackDxDr;
      const double recoY = recoStateY + radiusStep * trackDyDr;
      const double recoZ = recoStateZ + radiusStep * trackDzDr;
      const double recoR = std::hypot(recoX, recoY);
      const double recoPhi = std::atan2(recoY, recoX);
      if (!std::isfinite(recoX) || !std::isfinite(recoY) ||
          !std::isfinite(recoZ) || !std::isfinite(recoR) ||
          !std::isfinite(recoPhi))
      {
        ++invalidProjectionStates;
        continue;
      }

      const double dx = recoX - truthX;
      const double dy = recoY - truthY;
      const double dz = recoZ - truthZ;
      const double dr = recoR - truthR;
      const double dphi = std::remainder(recoPhi - truthPhi, 2. * TMath::Pi());
      const double dpx = recoStatePx - truth.px->at(truthIndex);
      const double dpy = recoStatePy - truth.py->at(truthIndex);
      const double dpz = recoStatePz - truth.pz->at(truthIndex);

      hDx->Fill(dx);
      hDy->Fill(dy);
      hDz->Fill(dz);
      hDpx->Fill(dpx);
      hDpy->Fill(dpy);
      hDpz->Fill(dpz);
      hDphi->Fill(dphi);
      hDr->Fill(dr);
      hTruthStateX->Fill(truthX);
      hTruthStateY->Fill(truthY);
      hTruthStateZ->Fill(truthZ);
      hRecoStateX->Fill(recoX);
      hRecoStateY->Fill(recoY);
      hRecoStateZ->Fill(recoZ);

      const float recoXError = reco.xError->at(recoIndex);
      const float recoYError = reco.yError->at(recoIndex);
      const float recoZError = reco.zError->at(recoIndex);
      if (std::isfinite(recoXError) && recoXError >= 0.F)
      {
        hRecoStateXError->Fill(recoXError);
        hRecoStateXErrorVsTruthR->Fill(truthR, recoXError);
        hRecoStateXErrorVsTruthPhi->Fill(truthPhi, recoXError);
        hRecoStateXErrorVsTruthZ->Fill(truthZ, recoXError);
      }
      if (std::isfinite(recoYError) && recoYError >= 0.F)
      {
        hRecoStateYError->Fill(recoYError);
        hRecoStateYErrorVsTruthR->Fill(truthR, recoYError);
        hRecoStateYErrorVsTruthPhi->Fill(truthPhi, recoYError);
        hRecoStateYErrorVsTruthZ->Fill(truthZ, recoYError);
      }
      if (std::isfinite(recoZError) && recoZError >= 0.F)
      {
        hRecoStateZError->Fill(recoZError);
        hRecoStateZErrorVsTruthR->Fill(truthR, recoZError);
        hRecoStateZErrorVsTruthPhi->Fill(truthPhi, recoZError);
        hRecoStateZErrorVsTruthZ->Fill(truthZ, recoZError);
      }
      if (validTruthTrackPt)
      {
        hDxVsTruthPt->Fill(truthTrackPt, dx);
        hDyVsTruthPt->Fill(truthTrackPt, dy);
        hDzVsTruthPt->Fill(truthTrackPt, dz);
        hDphiVsTruthPt->Fill(truthTrackPt, dphi);
      }
      if (truthPtBin >= 0)
      {
        hDxTruthPtBins[truthPtBin]->Fill(dx);
        hDyTruthPtBins[truthPtBin]->Fill(dy);
        hDzTruthPtBins[truthPtBin]->Fill(dz);
        hDphiTruthPtBins[truthPtBin]->Fill(dphi);
      }
      hDxVsLayer->Fill(layer, dx);
      hDyVsLayer->Fill(layer, dy);
      hDzVsLayer->Fill(layer, dz);
      hDphiVsTruthR->Fill(truthR, dphi);
      hDphiVsTruthZ->Fill(truthZ, dphi);
      hDphiVsTruthPhi->Fill(truthPhi, dphi);
      hDzVsTruthR->Fill(truthR, dz);
      hDzVsTruthZ->Fill(truthZ, dz);
      hDzVsTruthPhi->Fill(truthPhi, dz);

      // Protect mean-residual profiles from large fit outliers. The ordinary
      // TH1/TH2 histograms above remain uncut so the tails are still visible.
      if (dphi >= -0.02F && dphi <= 0.02F)
      {
        pDphiVsTruthR->Fill(truthR, dphi);
        pDphiVsTruthZ->Fill(truthZ, dphi);
        pDphiVsTruthPhi->Fill(truthPhi, dphi);
        pDphiTruthRPhi->Fill(truthR, truthPhi, dphi);
        pDphiTruthZPhi->Fill(truthZ, truthPhi, dphi);
        pDphiTruthZR->Fill(truthZ, truthR, dphi);
        ++dphiProfileStates;
      }

      if (dz >= -1.F && dz <= 1.F)
      {
        pDzVsTruthR->Fill(truthR, dz);
        pDzVsTruthZ->Fill(truthZ, dz);
        pDzVsTruthPhi->Fill(truthPhi, dz);
        pDzTruthRPhi->Fill(truthR, truthPhi, dz);
        pDzTruthZPhi->Fill(truthZ, truthPhi, dz);
        pDzTruthZR->Fill(truthZ, truthR, dz);
        ++dzProfileStates;
      }
      ++matchedStates;

      // ==============================================================
      // Add layer selections and residual histograms here. For example:
      //   if (layer == 55) { hDxLayer55->Fill(dx); }
      // Other available values include pathlength, position and momentum
      // from both truthIndex and recoIndex.
      // ==============================================================
      (void) layer;
    }

    if ((entry + 1) % 10000 == 0)
    {
      std::cout << "Processed " << entry + 1 << " / " << entriesToProcess << " entries" << std::endl;
    }
  }

  std::cout << "Input: " << inputFileName << std::endl;
  std::cout << "Processed entries: " << entriesToProcess << std::endl;
  std::cout << "Entries with inconsistent vector sizes: " << badVectorEntries << std::endl;
  std::cout << "Matched truth/reco states: " << matchedStates << std::endl;
  std::cout << "States in dphi profiles (-0.02 <= dphi <= 0.02 rad): "
            << dphiProfileStates << std::endl;
  std::cout << "States in dz profiles (-1 <= dz <= 1 cm): "
            << dzProfileStates << std::endl;
  std::cout << "Track-pair entries with invalid pT: " << invalidTrackPtEntries << std::endl;
  std::cout << "Matched states skipped because reco-to-cluster-radius projection was invalid: "
            << invalidProjectionStates << std::endl;

  outputFile.cd();
  hDx->Write();
  hDy->Write();
  hDz->Write();
  hDpx->Write();
  hDpy->Write();
  hDpz->Write();
  hDphi->Write();
  hDr->Write();
  hTruthStateRMinusRecoClusterR->Write();
  hTruthStateX->Write();
  hTruthStateY->Write();
  hTruthStateZ->Write();
  hRecoStateX->Write();
  hRecoStateY->Write();
  hRecoStateZ->Write();
  hRecoStateXError->Write();
  hRecoStateYError->Write();
  hRecoStateZError->Write();
  hRecoStateXErrorVsTruthR->Write();
  hRecoStateXErrorVsTruthPhi->Write();
  hRecoStateXErrorVsTruthZ->Write();
  hRecoStateYErrorVsTruthR->Write();
  hRecoStateYErrorVsTruthPhi->Write();
  hRecoStateYErrorVsTruthZ->Write();
  hRecoStateZErrorVsTruthR->Write();
  hRecoStateZErrorVsTruthPhi->Write();
  hRecoStateZErrorVsTruthZ->Write();
  hTruthTrackPt->Write();
  hRecoTrackPt->Write();
  hDeltaTrackPt->Write();
  hDxVsTruthPt->Write();
  hDyVsTruthPt->Write();
  hDzVsTruthPt->Write();
  hDphiVsTruthPt->Write();
  for (std::size_t ptBin = 0; ptBin < hDxTruthPtBins.size(); ++ptBin)
  {
    hDxTruthPtBins[ptBin]->Write();
    hDyTruthPtBins[ptBin]->Write();
    hDzTruthPtBins[ptBin]->Write();
    hDphiTruthPtBins[ptBin]->Write();
  }
  hDxVsLayer->Write();
  hDyVsLayer->Write();
  hDzVsLayer->Write();
  hDphiVsTruthR->Write();
  hDphiVsTruthZ->Write();
  hDphiVsTruthPhi->Write();
  hDzVsTruthR->Write();
  hDzVsTruthZ->Write();
  hDzVsTruthPhi->Write();
  pDphiVsTruthR->Write();
  pDphiVsTruthZ->Write();
  pDphiVsTruthPhi->Write();
  pDzVsTruthR->Write();
  pDzVsTruthZ->Write();
  pDzVsTruthPhi->Write();
  pDphiTruthRPhi->Write();
  pDzTruthRPhi->Write();
  pDphiTruthZPhi->Write();
  pDzTruthZPhi->Write();
  pDphiTruthZR->Write();
  pDzTruthZR->Write();

  // Divided canvases need extra room on the right for the COLZ palette,
  // including its numerical labels and Z-axis title.
  const auto drawColz = [](TH1* histogram)
  {
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.20);
    gPad->SetBottomMargin(0.16);
    histogram->GetXaxis()->SetTitleOffset(1.05);
    histogram->GetYaxis()->SetTitleOffset(1.05);
    histogram->GetZaxis()->SetLabelSize(0.04);
    histogram->GetZaxis()->SetTitleOffset(1.25);
    histogram->Draw("COLZ");
  };

  const auto drawOneDimensional = [](TH1* histogram)
  {
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.10);
    histogram->GetXaxis()->SetTitleOffset(1.05);
    histogram->GetYaxis()->SetTitleOffset(1.15);
    histogram->Draw("HIST");
  };

  // Example canvas. Replace this block with your preferred plotting style.
  auto* canvas = new TCanvas("c_position_residuals", "Position residuals", 1200, 1000);
  canvas->Divide(2, 2);
  canvas->cd(1);
  hDx->Draw();
  canvas->cd(2);
  hDy->Draw();
  canvas->cd(3);
  hDz->Draw();
  canvas->cd(4);
  hDphi->Draw();
  canvas->Write();

  // Draw the radial-position residual on its own canvas.
  auto* drCanvas = new TCanvas("c_delta_r", "Radial position residual", 800, 650);
  drCanvas->SetLeftMargin(0.14);
  drCanvas->SetBottomMargin(0.14);
  hDr->GetXaxis()->SetTitleOffset(1.05);
  hDr->GetYaxis()->SetTitleOffset(1.15);
  hDr->Draw();
  drCanvas->Write();

  auto* truthClusterRadiusCanvas = new TCanvas(
      "c_truth_state_r_minus_reco_cluster_r",
      "Truth-state radius minus reconstructed-cluster radius", 800, 650);
  drawOneDimensional(hTruthStateRMinusRecoClusterR);
  truthClusterRadiusCanvas->Write();

  auto* statePositionCanvas = new TCanvas(
      "c_state_positions", "Truth and reconstructed state positions", 1800, 550);
  statePositionCanvas->Divide(3, 1);
  const std::array<TH1*, 3> truthStatePositionHistograms = {
      hTruthStateX, hTruthStateY, hTruthStateZ};
  const std::array<TH1*, 3> recoStatePositionHistograms = {
      hRecoStateX, hRecoStateY, hRecoStateZ};
  const std::array<const char*, 3> statePositionAxisTitles = {
      "State x [cm]", "State y [cm]", "State z [cm]"};
  for (std::size_t coordinate = 0; coordinate < truthStatePositionHistograms.size(); ++coordinate)
  {
    statePositionCanvas->cd(coordinate + 1);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.10);
    auto* truthHistogram = truthStatePositionHistograms[coordinate];
    auto* recoHistogram = recoStatePositionHistograms[coordinate];
    truthHistogram->SetStats(false);
    recoHistogram->SetStats(false);
    truthHistogram->SetLineColor(kRed + 1);
    recoHistogram->SetLineColor(kBlue + 1);
    truthHistogram->SetLineWidth(2);
    recoHistogram->SetLineWidth(2);
    truthHistogram->GetXaxis()->SetTitle(statePositionAxisTitles[coordinate]);
    truthHistogram->GetYaxis()->SetTitle("Matched TPC states");
    truthHistogram->GetXaxis()->SetTitleOffset(1.05);
    truthHistogram->GetYaxis()->SetTitleOffset(1.15);
    truthHistogram->SetMaximum(
        1.18 * std::max(truthHistogram->GetMaximum(), recoHistogram->GetMaximum()));
    truthHistogram->Draw("HIST");
    recoHistogram->Draw("HIST SAME");
    auto* legend = new TLegend(0.65, 0.73, 0.87, 0.87);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(truthHistogram, "Truth state", "l");
    legend->AddEntry(recoHistogram, "Reco state", "l");
    legend->Draw();
  }
  statePositionCanvas->Write();

  auto* statePositionErrorCanvas = new TCanvas(
      "c_reco_state_position_errors", "Reconstructed state position uncertainties", 1800, 550);
  statePositionErrorCanvas->Divide(3, 1);
  statePositionErrorCanvas->cd(1);
  drawOneDimensional(hRecoStateXError);
  statePositionErrorCanvas->cd(2);
  drawOneDimensional(hRecoStateYError);
  statePositionErrorCanvas->cd(3);
  drawOneDimensional(hRecoStateZError);
  statePositionErrorCanvas->Write();

  auto* statePositionErrorVsTruthPositionCanvas = new TCanvas(
      "c_reco_state_position_errors_vs_truth_position",
      "Reconstructed state position uncertainties versus truth position", 1800, 1350);
  statePositionErrorVsTruthPositionCanvas->Divide(3, 3);
  const std::array<TH2*, 9> statePositionErrorMaps = {
      hRecoStateXErrorVsTruthR, hRecoStateXErrorVsTruthPhi, hRecoStateXErrorVsTruthZ,
      hRecoStateYErrorVsTruthR, hRecoStateYErrorVsTruthPhi, hRecoStateYErrorVsTruthZ,
      hRecoStateZErrorVsTruthR, hRecoStateZErrorVsTruthPhi, hRecoStateZErrorVsTruthZ};
  for (std::size_t index = 0; index < statePositionErrorMaps.size(); ++index)
  {
    statePositionErrorVsTruthPositionCanvas->cd(index + 1);
    drawColz(statePositionErrorMaps[index]);
  }
  statePositionErrorVsTruthPositionCanvas->Write();

  auto* trackPtCanvas = new TCanvas(
      "c_track_pt", "Truth and reconstructed track pT", 850, 700);
  trackPtCanvas->SetLeftMargin(0.14);
  trackPtCanvas->SetBottomMargin(0.14);
  hTruthTrackPt->SetStats(false);
  hRecoTrackPt->SetStats(false);
  hTruthTrackPt->SetLineColor(kRed + 1);
  hRecoTrackPt->SetLineColor(kBlue + 1);
  hTruthTrackPt->SetLineWidth(2);
  hRecoTrackPt->SetLineWidth(2);
  hTruthTrackPt->SetMaximum(
      1.15 * std::max(hTruthTrackPt->GetMaximum(), hRecoTrackPt->GetMaximum()));
  hTruthTrackPt->Draw("HIST");
  hRecoTrackPt->Draw("HIST SAME");
  auto* trackPtLegend = new TLegend(0.62, 0.72, 0.86, 0.86);
  trackPtLegend->SetBorderSize(0);
  trackPtLegend->SetFillStyle(0);
  trackPtLegend->AddEntry(hTruthTrackPt, "Truth track", "l");
  trackPtLegend->AddEntry(hRecoTrackPt, "Reco track", "l");
  trackPtLegend->Draw();
  trackPtCanvas->Write();

  auto* deltaTrackPtCanvas = new TCanvas(
      "c_delta_track_pt", "Track transverse-momentum residual", 850, 700);
  drawOneDimensional(hDeltaTrackPt);
  deltaTrackPtCanvas->Write();

  auto* residualVsTruthPtCanvas = new TCanvas(
      "c_residuals_vs_truth_pt", "State residuals versus truth track pT", 1400, 1100);
  residualVsTruthPtCanvas->Divide(2, 2);
  residualVsTruthPtCanvas->cd(1);
  drawColz(hDxVsTruthPt);
  residualVsTruthPtCanvas->cd(2);
  drawColz(hDyVsTruthPt);
  residualVsTruthPtCanvas->cd(3);
  drawColz(hDzVsTruthPt);
  residualVsTruthPtCanvas->cd(4);
  drawColz(hDphiVsTruthPt);
  residualVsTruthPtCanvas->Write();

  auto* residualTruthPtBinsCanvas = new TCanvas(
      "c_residuals_truth_pt_bins", "State residuals in truth track pT intervals", 1800, 1350);
  residualTruthPtBinsCanvas->Divide(4, 3);
  const std::array<const char*, 3> truthPtBinLabels = {
      "0.2 #leq p_{T,truth} < 0.6 GeV/c",
      "0.6 #leq p_{T,truth} #leq 1.0 GeV/c",
      "p_{T,truth} > 1.0 GeV/c"};
  const auto drawPtBinnedResidual = [&drawOneDimensional](TH1* histogram, const char* label)
  {
    drawOneDimensional(histogram);
    auto* binLabel = new TLatex(0.20, 0.84, label);
    binLabel->SetNDC();
    binLabel->SetTextSize(0.045);
    binLabel->Draw();
  };
  for (std::size_t ptBin = 0; ptBin < hDxTruthPtBins.size(); ++ptBin)
  {
    residualTruthPtBinsCanvas->cd(4 * ptBin + 1);
    drawPtBinnedResidual(hDxTruthPtBins[ptBin], truthPtBinLabels[ptBin]);
    residualTruthPtBinsCanvas->cd(4 * ptBin + 2);
    drawPtBinnedResidual(hDyTruthPtBins[ptBin], truthPtBinLabels[ptBin]);
    residualTruthPtBinsCanvas->cd(4 * ptBin + 3);
    drawPtBinnedResidual(hDzTruthPtBins[ptBin], truthPtBinLabels[ptBin]);
    residualTruthPtBinsCanvas->cd(4 * ptBin + 4);
    drawPtBinnedResidual(hDphiTruthPtBins[ptBin], truthPtBinLabels[ptBin]);
  }
  residualTruthPtBinsCanvas->Write();

  auto* layerCanvas = new TCanvas(
      "c_position_residuals_vs_layer", "Position residuals versus layer", 1800, 550);
  layerCanvas->Divide(3, 1);
  layerCanvas->cd(1);
  drawColz(hDxVsLayer);
  layerCanvas->cd(2);
  drawColz(hDyVsLayer);
  layerCanvas->cd(3);
  drawColz(hDzVsLayer);
  layerCanvas->Write();

  auto* truthPositionCanvas = new TCanvas(
      "c_residuals_vs_truth_position", "Residuals versus truth position", 1800, 1000);
  truthPositionCanvas->Divide(3, 2);
  truthPositionCanvas->cd(1);
  drawColz(hDphiVsTruthR);
  truthPositionCanvas->cd(2);
  drawColz(hDphiVsTruthZ);
  truthPositionCanvas->cd(3);
  drawColz(hDphiVsTruthPhi);
  truthPositionCanvas->cd(4);
  drawColz(hDzVsTruthR);
  truthPositionCanvas->cd(5);
  drawColz(hDzVsTruthZ);
  truthPositionCanvas->cd(6);
  drawColz(hDzVsTruthPhi);
  truthPositionCanvas->Write();

  auto* meanResidualCanvas = new TCanvas(
      "c_mean_residuals_vs_truth_position",
      "Mean residuals versus truth position", 1800, 1000);
  meanResidualCanvas->Divide(3, 2);
  meanResidualCanvas->cd(1);
  pDphiVsTruthR->Draw("HIST");
  meanResidualCanvas->cd(2);
  pDphiVsTruthZ->Draw("HIST");
  meanResidualCanvas->cd(3);
  pDphiVsTruthPhi->Draw("HIST");
  meanResidualCanvas->cd(4);
  pDzVsTruthR->Draw("HIST");
  meanResidualCanvas->cd(5);
  pDzVsTruthZ->Draw("HIST");
  meanResidualCanvas->cd(6);
  pDzVsTruthPhi->Draw("HIST");
  meanResidualCanvas->Write();

  auto* residualMapCanvas = new TCanvas(
      "c_residual_maps_truth_r_phi", "Mean residual maps in truth coordinate space", 1800, 1000);
  residualMapCanvas->Divide(3, 2);
  residualMapCanvas->cd(1);
  drawColz(pDphiTruthRPhi);
  residualMapCanvas->cd(2);
  drawColz(pDphiTruthZPhi);
  residualMapCanvas->cd(3);
  drawColz(pDphiTruthZR);
  residualMapCanvas->cd(4);
  drawColz(pDzTruthRPhi);
  residualMapCanvas->cd(5);
  drawColz(pDzTruthZPhi);
  residualMapCanvas->cd(6);
  drawColz(pDzTruthZR);
  residualMapCanvas->Write();

  // Export canvases as standalone figures. The input filename is included in
  // each output name so Acts and GenFit plots do not overwrite one another.
  gSystem->mkdir(figureDirectory, true);
  TString inputTag = gSystem->BaseName(inputFileName);
  if (inputTag.EndsWith(".root"))
  {
    inputTag.Remove(inputTag.Length() - 5);
  }

  const TString figurePrefix = TString::Format("%s/%s", figureDirectory, inputTag.Data());
  canvas->SaveAs(figurePrefix + "_position_residuals.png");
  canvas->SaveAs(figurePrefix + "_position_residuals.pdf");
  drCanvas->SaveAs(figurePrefix + "_delta_r.png");
  drCanvas->SaveAs(figurePrefix + "_delta_r.pdf");
  truthClusterRadiusCanvas->SaveAs(
      figurePrefix + "_truth_state_r_minus_reco_cluster_r.png");
  truthClusterRadiusCanvas->SaveAs(
      figurePrefix + "_truth_state_r_minus_reco_cluster_r.pdf");
  statePositionCanvas->SaveAs(figurePrefix + "_state_positions.png");
  statePositionCanvas->SaveAs(figurePrefix + "_state_positions.pdf");
  statePositionErrorCanvas->SaveAs(figurePrefix + "_reco_state_position_errors.png");
  statePositionErrorCanvas->SaveAs(figurePrefix + "_reco_state_position_errors.pdf");
  statePositionErrorVsTruthPositionCanvas->SaveAs(
      figurePrefix + "_reco_state_position_errors_vs_truth_position.png");
  statePositionErrorVsTruthPositionCanvas->SaveAs(
      figurePrefix + "_reco_state_position_errors_vs_truth_position.pdf");
  trackPtCanvas->SaveAs(figurePrefix + "_track_pt.png");
  trackPtCanvas->SaveAs(figurePrefix + "_track_pt.pdf");
  deltaTrackPtCanvas->SaveAs(figurePrefix + "_delta_track_pt.png");
  deltaTrackPtCanvas->SaveAs(figurePrefix + "_delta_track_pt.pdf");
  residualVsTruthPtCanvas->SaveAs(figurePrefix + "_residuals_vs_truth_pt.png");
  residualVsTruthPtCanvas->SaveAs(figurePrefix + "_residuals_vs_truth_pt.pdf");
  residualTruthPtBinsCanvas->SaveAs(figurePrefix + "_residuals_truth_pt_bins.png");
  residualTruthPtBinsCanvas->SaveAs(figurePrefix + "_residuals_truth_pt_bins.pdf");
  layerCanvas->SaveAs(figurePrefix + "_position_residuals_vs_layer.png");
  layerCanvas->SaveAs(figurePrefix + "_position_residuals_vs_layer.pdf");
  truthPositionCanvas->SaveAs(figurePrefix + "_residuals_vs_truth_position.png");
  truthPositionCanvas->SaveAs(figurePrefix + "_residuals_vs_truth_position.pdf");
  meanResidualCanvas->SaveAs(figurePrefix + "_mean_residuals_vs_truth_position.png");
  meanResidualCanvas->SaveAs(figurePrefix + "_mean_residuals_vs_truth_position.pdf");
  residualMapCanvas->SaveAs(figurePrefix + "_residual_maps_truth_r_phi.png");
  residualMapCanvas->SaveAs(figurePrefix + "_residual_maps_truth_r_phi.pdf");

  outputFile.Close();
  std::cout << "Wrote " << outputFileName << std::endl;
}
