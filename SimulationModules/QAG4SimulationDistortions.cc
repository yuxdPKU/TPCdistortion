#include "QAG4SimulationDistortions.h"

#include <micromegas/MicromegasDefs.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Acts/Definitions/Algebra.hpp>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <utility>  // for pair
#include <vector>

namespace
{

  // square
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  // radius
  template <class T>
  T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  template <class T>
  constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }

    return phi;
  }

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }
}  // namespace

//____________________________________________________________________________..
QAG4SimulationDistortions::QAG4SimulationDistortions(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::Init(PHCompositeNode* /*unused*/)
{
  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_states = 0;
  m_accepted_states = 0;

  // histogram manager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_beta = new TH2F(std::string(get_histo_prefix() + "betadz").c_str(), ";tan#beta; #Deltaz [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hm->registerHisto(h_beta);

  h_alpha = new TH2F(std::string(get_histo_prefix() + "alphardphi").c_str(), ";tan#alpha; r#Delta#phi [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hm->registerHisto(h_alpha);

  h_rphiResid = new TH2F(std::string(get_histo_prefix() + "rphiResid").c_str(), ";r [cm]; #Deltar#phi [cm]", 60, 20, 80, 500, -2, 2);
  hm->registerHisto(h_rphiResid);

  h_zResid = new TH2F(std::string(get_histo_prefix() + "zResid").c_str(), ";z [cm]; #Deltaz [cm]", 200, -100, 100, 1000, -2, 2);
  hm->registerHisto(h_zResid);

  h_etaResid = new TH2F(std::string(get_histo_prefix() + "etaResid").c_str(), ";#eta;#Delta#eta", 20, -1, 1, 500, -0.2, 0.2);
  hm->registerHisto(h_etaResid);

  h_etaResidLayer = new TH2F(std::string(get_histo_prefix() + "etaResidLayer").c_str(), ";r [cm]; #Delta#eta", 60, 20, 80, 500, -0.2, 0.2);
  hm->registerHisto(h_etaResidLayer);

  h_zResidLayer = new TH2F(std::string(get_histo_prefix() + "zResidLayer").c_str(), ";r [cm]; #Deltaz [cm]", 60, 20, 80, 1000, -2, 2);
  hm->registerHisto(h_zResidLayer);

  h_deltarphi_layer = new TH2F(std::string(get_histo_prefix() + "deltarphi_layer").c_str(), ";layer; r.#Delta#phi_{cluster-track} (cm)", 57, 0, 57, 500, -2, 2);
  hm->registerHisto(h_deltarphi_layer);

  h_deltaz_layer = new TH2F(std::string(get_histo_prefix() + "deltaz_layer").c_str(), ";layer; #Deltaz_{cluster-track} (cm)", 57, 0, 57, 100, -2, 2);
  hm->registerHisto(h_deltaz_layer);

  h_statez_pulls = new TH2F(std::string(get_histo_prefix() + "statez_pulls").c_str(), ";layer; #Deltaz_{cluster-track}/#sigma_{z}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h_statez_pulls);

  h_staterphi_pulls = new TH2F(std::string(get_histo_prefix() + "staterphi_pulls").c_str(), ";layer; #Deltar#phi_{cluster-track}/#sigma_{rphi}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h_staterphi_pulls);

  h_clusz_pulls = new TH2F(std::string(get_histo_prefix() + "clusz_pulls").c_str(), ";layer; #Deltaz_{cluster-track}/#sigma_{z}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h_clusz_pulls);

  h_clusrphi_pulls = new TH2F(std::string(get_histo_prefix() + "clusrphi_pulls").c_str(), ";layer; #Deltar#phi_{cluster-track}/#sigma_{r#phi}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h_clusrphi_pulls);

  h_nstates_vs_nclus = new TH2F(std::string(get_histo_prefix() + "nstates_vs_nclus").c_str(), ";Number of states on track; Number of clusters on track", 57, 0, 57, 57, 0, 57);
  hm->registerHisto(h_nstates_vs_nclus);

  h_tpot_deltarphi = new TH1F(std::string(get_histo_prefix() + "tpot_deltarphi").c_str(), ";TPOT  r.#Delta#phi_{cluster-track} (cm); Number of TPOT clusters", 100, -1, 1);
  hm->registerHisto(h_tpot_deltarphi);

  h_tpot_deltaz = new TH1F(std::string(get_histo_prefix() + "tpot_deltaz").c_str(), ";TPOT  #Deltaz_{cluster-track} (cm); Number of TPOT clusters", 100, -2, 2);
  hm->registerHisto(h_tpot_deltaz);

  t_tree = new TTree(std::string(get_histo_prefix() + "residTree").c_str(), "tpc residual info");
  t_tree->Branch("tanAlpha", &m_tanAlpha, "tanAlpha/F");
  t_tree->Branch("tanBeta", &m_tanBeta, "tanBeta/F");
  t_tree->Branch("drphi", &m_drphi, "drphi/F");
  t_tree->Branch("dz", &m_dz, "dz/F");
  t_tree->Branch("clusR", &m_clusR, "clusR/F");
  t_tree->Branch("clusPhi", &m_clusPhi, "clusPhi/F");
  t_tree->Branch("clusZ", &m_clusZ, "clusZ/F");
  t_tree->Branch("clusEta", &m_clusEta, "clusEta/F");
  t_tree->Branch("stateR", &m_stateR, "stateR/F");
  t_tree->Branch("statePhi", &m_statePhi, "statePhi/F");
  t_tree->Branch("stateZ", &m_stateZ, "stateZ/F");
  t_tree->Branch("stateEta", &m_stateEta, "stateEta/F");
  t_tree->Branch("stateRPhiErr", &m_stateRPhiErr, "stateRPhiErr/F");
  t_tree->Branch("stateZErr", &m_stateZErr, "stateZErr/F");
  t_tree->Branch("clusRPhiErr", &m_clusRPhiErr, "clusRPhiErr/F");
  t_tree->Branch("clusZErr", &m_clusZErr, "clusZErr/F");
  t_tree->Branch("cluskey", &m_cluskey, "cluskey/l");
  t_tree->Branch("event", &m_event, "event/I");

  t_tree->Branch("layer", &m_layer, "layer/I");
  t_tree->Branch("statePt", &m_statePt, "statePt/F");
  t_tree->Branch("statePz", &m_statePz, "statePz/F");

  t_tree->Branch("trackPt", &m_trackPt, "trackPt/F");
  t_tree->Branch("charge", &m_charge, "charge/I");
  t_tree->Branch("crossing", &m_crossing, "crossing/I");
  t_tree->Branch("trackdEdx", &m_trackdEdx, "trackdEdx/F");

  if (m_detailedInfo)
  {
  t_tree->Branch("cluskey_mvtx", &m_cluskey_mvtx);
  t_tree->Branch("layer_mvtx", &m_layer_mvtx);
  t_tree->Branch("sclusgx_mvtx", &m_sclusgx_mvtx);
  t_tree->Branch("sclusgy_mvtx", &m_sclusgy_mvtx);
  t_tree->Branch("sclusgz_mvtx", &m_sclusgz_mvtx);
  t_tree->Branch("sclusgr_mvtx", &m_sclusgr_mvtx);
  t_tree->Branch("sclusphi_mvtx", &m_sclusphi_mvtx);
  t_tree->Branch("scluseta_mvtx", &m_scluseta_mvtx);
  t_tree->Branch("stategx_mvtx", &m_stategx_mvtx);
  t_tree->Branch("stategy_mvtx", &m_stategy_mvtx);
  t_tree->Branch("stategz_mvtx", &m_stategz_mvtx);
  t_tree->Branch("stategr_mvtx", &m_stategr_mvtx);
  t_tree->Branch("statephi_mvtx", &m_statephi_mvtx);
  t_tree->Branch("stateeta_mvtx", &m_stateeta_mvtx);

  t_tree->Branch("cluskey_intt", &m_cluskey_intt);
  t_tree->Branch("layer_intt", &m_layer_intt);
  t_tree->Branch("sclusgx_intt", &m_sclusgx_intt);
  t_tree->Branch("sclusgy_intt", &m_sclusgy_intt);
  t_tree->Branch("sclusgz_intt", &m_sclusgz_intt);
  t_tree->Branch("sclusgr_intt", &m_sclusgr_intt);
  t_tree->Branch("sclusphi_intt", &m_sclusphi_intt);
  t_tree->Branch("scluseta_intt", &m_scluseta_intt);
  t_tree->Branch("stategx_intt", &m_stategx_intt);
  t_tree->Branch("stategy_intt", &m_stategy_intt);
  t_tree->Branch("stategz_intt", &m_stategz_intt);
  t_tree->Branch("stategr_intt", &m_stategr_intt);
  t_tree->Branch("statephi_intt", &m_statephi_intt);
  t_tree->Branch("stateeta_intt", &m_stateeta_intt);

  t_tree->Branch("cluskey_tpot", &m_cluskey_tpot);
  t_tree->Branch("layer_tpot", &m_layer_tpot);
  t_tree->Branch("segtype_tpot", &m_segtype_tpot);
  t_tree->Branch("tileid_tpot", &m_tileid_tpot);
  t_tree->Branch("sclusgx_tpot", &m_sclusgx_tpot);
  t_tree->Branch("sclusgy_tpot", &m_sclusgy_tpot);
  t_tree->Branch("sclusgz_tpot", &m_sclusgz_tpot);
  t_tree->Branch("sclusgr_tpot", &m_sclusgr_tpot);
  t_tree->Branch("sclusphi_tpot", &m_sclusphi_tpot);
  t_tree->Branch("scluseta_tpot", &m_scluseta_tpot);
  t_tree->Branch("stategx_tpot", &m_stategx_tpot);
  t_tree->Branch("stategy_tpot", &m_stategy_tpot);
  t_tree->Branch("stategz_tpot", &m_stategz_tpot);
  t_tree->Branch("stategr_tpot", &m_stategr_tpot);
  t_tree->Branch("statephi_tpot", &m_statephi_tpot);
  t_tree->Branch("stateeta_tpot", &m_stateeta_tpot);
  }

  hm->registerHisto(t_tree);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::InitRun(PHCompositeNode* topNode)
{
  // track map
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // cluster map
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // find tpc geometry
  m_tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  // load geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  // load distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);
  if (m_disable_module_edge_corr)
  {
    m_globalPositionWrapper.set_enable_module_edge_corr(false);
  }
  if (m_disable_static_corr)
  {
    m_globalPositionWrapper.set_enable_static_corr(false);
  }
  if (m_disable_average_corr)
  {
    m_globalPositionWrapper.set_enable_average_corr(false);
  }
  if (m_disable_fluctuation_corr)
  {
    m_globalPositionWrapper.set_enable_fluctuation_corr(false);
  }

  if (!m_trackMap || !m_clusterContainer || !m_tGeometry)
  {
    std::cout << PHWHERE << "Necessary distortion container not on node tree. Bailing."
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::process_event(PHCompositeNode* /*unused*/)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  std::cout << "QAG4SimulationDistortions::process_event - tracks: " << m_trackMap->size() << std::endl;
  for (const auto& [key, track] : *m_trackMap)
  {
    // total track counter
    ++m_total_tracks;

    // get track crossing and check
    const auto crossing = track->get_crossing();
    if (crossing == std::numeric_limits<short>::max())
    {
      std::cout << "QAG4SimulationDistortions::process_event - invalid crossing. Track skipped." << std::endl;
      continue;
    }

    // check track quality
    if (!checkTrack(track))
    {
      if (Verbosity() > 1)
      {
        std::cout << "failed track selection" << std::endl;
      }
      continue;
    }
    if (Verbosity() > 1)
    {
      std::cout << "pass track selection" << std::endl;
    }

    // get seeeds
    auto* tpcSeed = track->get_tpc_seed();
    auto* siliconSeed = track->get_silicon_seed();

    /// Should have never been added to the map...
    if (!tpcSeed || !siliconSeed)
    {
      continue;
    }

    if (Verbosity() > 0)
    {
      std::cout << " tpc seed " << tpcSeed->get_tpc_seed_index()
                << " with si seed " << siliconSeed->get_silicon_seed_index()
                << " crossing " << siliconSeed->get_crossing()
                << std::endl;
    }

    if (m_detailedInfo) {get_MvtxInttTpot_info(track);}

    // accepted track counter
    ++m_accepted_tracks;

    for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
    {
      ++m_total_states;

      auto& state = iter->second;
      const auto ckey = state->get_cluskey();
      const auto trkrId = TrkrDefs::getTrkrId(ckey);

      if (trkrId != TrkrDefs::tpcId)
      {
        continue;
      }

      ++m_accepted_states;

      auto* cluster = m_clusterContainer->findCluster(ckey);

      const auto clusGlobPosition = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, crossing);

      const float clusR = get_r(clusGlobPosition(0), clusGlobPosition(1));
      const float clusPhi = std::atan2(clusGlobPosition(1), clusGlobPosition(0));
      const float clusZ = clusGlobPosition(2);

      // cluster errors
      const float clusRPhiErr = cluster->getRPhiError();
      const float clusZErr = cluster->getZError();

      const Acts::Vector3 stateGlobPosition = Acts::Vector3(state->get_x(),
                                                            state->get_y(),
                                                            state->get_z());
      const Acts::Vector3 stateGlobMom = Acts::Vector3(state->get_px(),
                                                       state->get_py(),
                                                       state->get_pz());

      const float stateRPhiErr = state->get_rphi_error();
      const float stateZErr = state->get_z_error();

      const float stateR = get_r(stateGlobPosition(0), stateGlobPosition(1));

      const auto dr = clusR - stateR;
      const auto trackDrDt = (stateGlobPosition(0) * stateGlobMom(0) + stateGlobPosition(1) * stateGlobMom(1)) / stateR;
      const auto trackDxDr = stateGlobMom(0) / trackDrDt;
      const auto trackDyDr = stateGlobMom(1) / trackDrDt;
      const auto trackDzDr = stateGlobMom(2) / trackDrDt;

      const auto trackX = stateGlobPosition(0) + dr * trackDxDr;
      const auto trackY = stateGlobPosition(1) + dr * trackDyDr;
      const auto trackZ = stateGlobPosition(2) + dr * trackDzDr;
      const float statePhi = std::atan2(trackY, trackX);
      const float stateZ = trackZ;

      // Calculate residuals
      const float drphi = clusR * deltaPhi(clusPhi - statePhi);
      const float dz = clusZ - stateZ;

      const auto trackPPhi = -stateGlobMom(0) * std::sin(statePhi) + stateGlobMom(1) * std::cos(statePhi);
      const auto trackPR = stateGlobMom(0) * std::cos(statePhi) + stateGlobMom(1) * std::sin(statePhi);
      const auto trackPZ = stateGlobMom(2);

      const auto trackAlpha = -trackPPhi / trackPR;
      const auto trackBeta = -trackPZ / trackPR;
      const auto trackEta = std::atanh(stateGlobMom(2) / stateGlobMom.norm());
      const auto clusEta = std::atanh(clusZ / clusGlobPosition.norm());

      h_alpha->Fill(trackAlpha, drphi);
      h_beta->Fill(trackBeta, dz);
      h_rphiResid->Fill(clusR, drphi);
      h_zResid->Fill(stateZ, dz);
      h_etaResid->Fill(trackEta, clusEta - trackEta);
      h_zResidLayer->Fill(clusR, dz);
      h_etaResidLayer->Fill(clusR, clusEta - trackEta);

      const auto layer = TrkrDefs::getLayer(ckey);
      h_deltarphi_layer->Fill(layer, drphi);
      h_deltaz_layer->Fill(layer, dz);

      h_statez_pulls->Fill(layer, dz / stateZErr);
      h_staterphi_pulls->Fill(layer, drphi / stateRPhiErr);
      h_clusz_pulls->Fill(layer, dz / clusZErr);
      h_clusrphi_pulls->Fill(layer, drphi / clusRPhiErr);

      m_tanAlpha = trackAlpha;
      m_tanBeta = trackBeta;
      m_drphi = drphi;
      m_dz = dz;
      m_clusR = clusR;
      m_clusPhi = clusPhi;
      m_clusZ = clusZ;
      m_statePhi = statePhi;
      m_stateZ = stateZ;
      m_stateR = stateR;
      m_stateRPhiErr = stateRPhiErr;
      m_stateZErr = stateZErr;
      m_clusRPhiErr = clusRPhiErr;
      m_clusZErr = clusZErr;
      m_cluskey = ckey;

      m_clusEta = clusEta;
      m_stateEta = trackEta;
      m_layer = layer;
      m_statePt = sqrt(stateGlobMom(0) * stateGlobMom(0) + stateGlobMom(1) * stateGlobMom(1));
      m_statePz = stateGlobMom(2);
      m_trackPt = track->get_pt();
      m_charge = track->get_charge();
      m_crossing = track->get_crossing();

      float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
      // These are randomly chosen layer thicknesses for the TPC, to get the
      // correct region thicknesses in an easy to pass way to the helper fxn
      layerThicknesses[0] = m_tpcGeom->GetLayerCellGeom(7)->get_thickness();
      layerThicknesses[1] = m_tpcGeom->GetLayerCellGeom(8)->get_thickness();
      layerThicknesses[2] = m_tpcGeom->GetLayerCellGeom(27)->get_thickness();
      layerThicknesses[3] = m_tpcGeom->GetLayerCellGeom(50)->get_thickness();

      m_trackdEdx = TrackAnalysisUtils::calc_dedx(tpcSeed, m_clusterContainer, m_tGeometry, layerThicknesses);

      t_tree->Fill();
    }
    int nstates = track->size_states();
    int nclus = (track->get_silicon_seed()->size_cluster_keys()) + (track->get_tpc_seed()->size_cluster_keys());
    h_nstates_vs_nclus->Fill(nstates, nclus);
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int QAG4SimulationDistortions::End(PHCompositeNode* /*topNode*/)
{
  // print counters
  std::cout
      << "QAG4SimulationDistortions::End -"
      << " track statistics total: " << m_total_tracks
      << " accepted: " << m_accepted_tracks
      << " fraction: " << 100. * m_accepted_tracks / m_total_tracks << "%"
      << std::endl;

  std::cout
      << "QAG4SimulationDistortions::End -"
      << " state statistics total: " << m_total_states
      << " accepted: " << m_accepted_states << " fraction: "
      << 100. * m_accepted_states / m_total_states << "%"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________________
bool QAG4SimulationDistortions::checkTrack(SvtxTrack* track)
{
  if (track->get_pt() < m_minPT)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt, tpc and micromegas hits
  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 3)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(cluster_keys) < 35)
  {
    return false;
  }

  const auto state_keys(get_state_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(state_keys) < 3)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(state_keys) < 2)
  {
    return false;
  }

  if (m_useMicromegas && checkTPOTResidual(track) == false)
  {
    return false;
  }

  return true;
}

bool QAG4SimulationDistortions::checkTPOTResidual(SvtxTrack* track)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  bool flag = true;

  int nTPOTcluster = 0;
  int nTPOTstate = 0;
  int TPOTtileID = -1;
  for (const auto& cluskey : get_cluster_keys(track))
  {
    // make sure cluster is from TPOT
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if (detId != TrkrDefs::micromegasId)
    {
      continue;
    }
    TPOTtileID = MicromegasDefs::getTileId(cluskey);
    nTPOTcluster++;

    auto* cluster = m_clusterContainer->findCluster(cluskey);

    SvtxTrackState* state = nullptr;

    // the track states from the Acts fit are fitted to fully corrected clusters, and are on the surface
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      if (stateckey == cluskey)
      {
        state = tstate;
        break;
      }
    }

    const auto layer = TrkrDefs::getLayer(cluskey);

    if (!state)
    {
      if (Verbosity() > 1)
      {
        std::cout << "   no state for cluster " << cluskey << "  in layer " << layer << std::endl;
      }
      continue;
    }
    nTPOTstate++;

    const auto crossing = track->get_crossing();
    assert(crossing != std::numeric_limits<short>::max());

    // calculate residuals with respect to cluster
    // Get all the relevant information for residual calculation
    const auto globClusPos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluskey, cluster, crossing);
    const double clusR = get_r(globClusPos(0), globClusPos(1));
    const double clusPhi = std::atan2(globClusPos(1), globClusPos(0));
    const double clusZ = globClusPos(2);

    const double globStateX = state->get_x();
    const double globStateY = state->get_y();
    const double globStateZ = state->get_z();
    const double globStatePx = state->get_px();
    const double globStatePy = state->get_py();
    const double globStatePz = state->get_pz();

    const double trackR = std::sqrt(square(globStateX) + square(globStateY));

    const double dr = clusR - trackR;
    const double trackDrDt = (globStateX * globStatePx + globStateY * globStatePy) / trackR;
    const double trackDxDr = globStatePx / trackDrDt;
    const double trackDyDr = globStatePy / trackDrDt;
    const double trackDzDr = globStatePz / trackDrDt;

    const double trackX = globStateX + dr * trackDxDr;
    const double trackY = globStateY + dr * trackDyDr;
    const double trackZ = globStateZ + dr * trackDzDr;
    const double trackPhi = std::atan2(trackY, trackX);

    // Calculate residuals
    // need to be calculated in local coordinates in the future
    const double drphi = clusR * deltaPhi(clusPhi - trackPhi);
    if (std::isnan(drphi))
    {
      continue;
    }

    const double dz = clusZ - trackZ;
    if (std::isnan(dz))
    {
      continue;
    }

    h_tpot_deltarphi->Fill(drphi);
    h_tpot_deltaz->Fill(dz);

    if (Verbosity() > 3)
    {
      std::cout << "QAG4SimulationDistortions::checkTPOTResidual -"
                << " drphi: " << drphi
                << " dz: " << dz
                << std::endl;
    }

    // check rphi residual for layer 55
    if (layer == 55 && std::fabs(drphi) > 0.1)
    {
      flag = false;
      break;
    }

    // check z residual for layer 56
    if (layer == 56 && std::fabs(dz) > 1)
    {
      flag = false;
      break;
    }
  }

  if (flag)
  {
    // SCOZ has a half dead tile
    // only require one TPOT cluster/state from SCOP
    if (TPOTtileID == 0)
    {
      if (nTPOTcluster < 1 || nTPOTstate < 1)
      {
        flag = false;
      }
    }
    else if (TPOTtileID > 0)
    {
      if (nTPOTcluster < 2 || nTPOTstate < 2)
      {
        flag = false;
      }
    }
    else if (TPOTtileID < 0)
    {
      flag = false;
    }
  }

  return flag;
}

std::vector<TrkrDefs::cluskey> QAG4SimulationDistortions::get_cluster_keys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}

std::vector<TrkrDefs::cluskey> QAG4SimulationDistortions::get_state_keys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState* tstate = state_iter->second;
    auto stateckey = tstate->get_cluskey();
    out.push_back(stateckey);
  }
  return out;
}

void QAG4SimulationDistortions::clearVector()
{
  m_cluskey_mvtx.clear();
  m_layer_mvtx.clear();
  m_sclusgx_mvtx.clear();
  m_sclusgy_mvtx.clear();
  m_sclusgz_mvtx.clear();
  m_sclusgr_mvtx.clear();
  m_sclusphi_mvtx.clear();
  m_scluseta_mvtx.clear();
  m_stategx_mvtx.clear();
  m_stategy_mvtx.clear();
  m_stategz_mvtx.clear();
  m_stategr_mvtx.clear();
  m_statephi_mvtx.clear();
  m_stateeta_mvtx.clear();

  m_cluskey_intt.clear();
  m_layer_intt.clear();
  m_sclusgx_intt.clear();
  m_sclusgy_intt.clear();
  m_sclusgz_intt.clear();
  m_sclusgr_intt.clear();
  m_sclusphi_intt.clear();
  m_scluseta_intt.clear();
  m_stategx_intt.clear();
  m_stategy_intt.clear();
  m_stategz_intt.clear();
  m_stategr_intt.clear();
  m_statephi_intt.clear();
  m_stateeta_intt.clear();

  m_cluskey_tpot.clear();
  m_layer_tpot.clear();
  m_segtype_tpot.clear();
  m_tileid_tpot.clear();
  m_sclusgx_tpot.clear();
  m_sclusgy_tpot.clear();
  m_sclusgz_tpot.clear();
  m_sclusgr_tpot.clear();
  m_sclusphi_tpot.clear();
  m_scluseta_tpot.clear();
  m_stategx_tpot.clear();
  m_stategy_tpot.clear();
  m_stategz_tpot.clear();
  m_stategr_tpot.clear();
  m_statephi_tpot.clear();
  m_stateeta_tpot.clear();
}

void QAG4SimulationDistortions::get_MvtxInttTpot_info(SvtxTrack* track)
{
  clearVector();
  for (const auto& ckey : get_cluster_keys(track))
  {
    auto cluster = m_clusterContainer->findCluster(ckey);

    unsigned int layer = TrkrDefs::getLayer(ckey);
    Acts::Vector3 glob = m_tGeometry->getGlobalPosition(ckey, cluster);
    float m_sclusgx = glob.x();
    float m_sclusgy = glob.y();
    float m_sclusgz = glob.z();
    float m_sclusgr = get_r(m_sclusgx, m_sclusgy);
    float m_sclusphi = atan2(glob.y(), glob.x());
    float m_scluseta = acos(glob.z() / std::sqrt(square(glob.x()) + square(glob.y()) + square(glob.z())));

    SvtxTrackState* state = nullptr;
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      if (stateckey == ckey)
      {
        state = tstate;
        break;
      }
    }

    float m_stategx = NAN;
    float m_stategy = NAN;
    float m_stategz = NAN;
    float m_stategr = NAN;
    float m_statephi = NAN;
    float m_stateeta = NAN;
    if (state)
    {
      Acts::Vector3 stateglob(state->get_x(), state->get_y(), state->get_z());
      m_stategx = stateglob.x();
      m_stategy = stateglob.y();
      m_stategz = stateglob.z();
      m_stategr = get_r(m_stategx, m_stategy);
      m_statephi = atan2(stateglob.y(), stateglob.x());
      m_stateeta = acos(stateglob.z() / std::sqrt(square(stateglob.x()) + square(stateglob.y()) + square(stateglob.z())));
    }

    if (TrkrDefs::getTrkrId(ckey)==TrkrDefs::mvtxId)
    {
      m_cluskey_mvtx.push_back(ckey);
      m_layer_mvtx.push_back(layer);
      m_sclusgx_mvtx.push_back(m_sclusgx);
      m_sclusgy_mvtx.push_back(m_sclusgy);
      m_sclusgz_mvtx.push_back(m_sclusgz);
      m_sclusgr_mvtx.push_back(m_sclusgr);
      m_sclusphi_mvtx.push_back(m_sclusphi);
      m_scluseta_mvtx.push_back(m_scluseta);
      m_stategx_mvtx.push_back(m_stategx);
      m_stategy_mvtx.push_back(m_stategy);
      m_stategz_mvtx.push_back(m_stategz);
      m_stategr_mvtx.push_back(m_stategr);
      m_statephi_mvtx.push_back(m_statephi);
      m_stateeta_mvtx.push_back(m_stateeta);
    }
    else if (TrkrDefs::getTrkrId(ckey)==TrkrDefs::inttId)
    {
      m_cluskey_intt.push_back(ckey);
      m_layer_intt.push_back(layer);
      m_sclusgx_intt.push_back(m_sclusgx);
      m_sclusgy_intt.push_back(m_sclusgy);
      m_sclusgz_intt.push_back(m_sclusgz);
      m_sclusgr_intt.push_back(m_sclusgr);
      m_sclusphi_intt.push_back(m_sclusphi);
      m_scluseta_intt.push_back(m_scluseta);
      m_stategx_intt.push_back(m_stategx);
      m_stategy_intt.push_back(m_stategy);
      m_stategz_intt.push_back(m_stategz);
      m_stategr_intt.push_back(m_stategr);
      m_statephi_intt.push_back(m_statephi);
      m_stateeta_intt.push_back(m_stateeta);
    }
    else if (TrkrDefs::getTrkrId(ckey)==TrkrDefs::micromegasId)
    {
      int m_segtype = (int) MicromegasDefs::getSegmentationType(ckey);
      int m_tileid = MicromegasDefs::getTileId(ckey);
      m_cluskey_tpot.push_back(ckey);
      m_layer_tpot.push_back(layer);
      m_segtype_tpot.push_back(m_segtype);
      m_tileid_tpot.push_back(m_tileid);
      m_sclusgx_tpot.push_back(m_sclusgx);
      m_sclusgy_tpot.push_back(m_sclusgy);
      m_sclusgz_tpot.push_back(m_sclusgz);
      m_sclusgr_tpot.push_back(m_sclusgr);
      m_sclusphi_tpot.push_back(m_sclusphi);
      m_scluseta_tpot.push_back(m_scluseta);
      m_stategx_tpot.push_back(m_stategx);
      m_stategy_tpot.push_back(m_stategy);
      m_stategz_tpot.push_back(m_stategz);
      m_stategr_tpot.push_back(m_stategr);
      m_statephi_tpot.push_back(m_statephi);
      m_stateeta_tpot.push_back(m_stateeta);
    }
  }
}
