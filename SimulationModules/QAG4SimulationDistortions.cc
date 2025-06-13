
#include "QAG4SimulationDistortions.h"

#include <fun4all/SubsysReco.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <micromegas/MicromegasDefs.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TTree.h>

#include <Acts/Definitions/Algebra.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>  // for pair
#include <vector>

namespace
{

  // square
  template <class T>
  inline constexpr T square(const T& x)
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
  inline constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    else if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }
    else
    {
      return phi;
    }
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
QAG4SimulationDistortions::~QAG4SimulationDistortions() = default;

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

  TH1* h(nullptr);

  h = new TH2F(TString(get_histo_prefix()) + "betadz", ";tan#beta; #Deltaz [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);

  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "alphardphi", ";tan#alpha; r#Delta#phi [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "rphiResid", ";r [cm]; #Deltar#phi [cm]", 60, 20, 80, 500, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "zResid", ";z [cm]; #Deltaz [cm]", 200, -100, 100, 1000, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "etaResid", ";#eta;#Delta#eta", 20, -1, 1, 500, -0.2, 0.2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "etaResidLayer", ";r [cm]; #Delta#eta", 60, 20, 80, 500, -0.2, 0.2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "zResidLayer", ";r [cm]; #Deltaz [cm]", 60, 20, 80, 1000, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "deltarphi_layer", ";layer; r.#Delta#phi_{cluster-track} (cm)", 57, 0, 57, 500, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "deltaz_layer", ";layer; #Deltaz_{cluster-track} (cm)", 57, 0, 57, 100, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "statez_pulls", ";layer; #Deltaz_{cluster-track}/#sigma_{z}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "staterphi_pulls", ";layer; #Deltar#phi_{cluster-track}/#sigma_{rphi}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "clusz_pulls", ";layer; #Deltaz_{cluster-track}/#sigma_{z}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "clusrphi_pulls", ";layer; #Deltar#phi_{cluster-track}/#sigma_{r#phi}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "nstates_vs_nclus", ";Number of states on track; Number of clusters on track", 57, 0, 57, 57, 0, 57);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "tpot_deltarphi", ";TPOT  r.#Delta#phi_{cluster-track} (cm); Number of TPOT clusters", 100, -1, 1);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "tpot_deltaz", ";TPOT  #Deltaz_{cluster-track} (cm); Number of TPOT clusters", 100, -2, 2);
  hm->registerHisto(h);

  TTree* t(nullptr);

  t = new TTree(TString(get_histo_prefix()) + "residTree", "tpc residual info");
  t->Branch("tanAlpha", &m_tanAlpha, "tanAlpha/F");
  t->Branch("tanBeta", &m_tanBeta, "tanBeta/F");
  t->Branch("drphi", &m_drphi, "drphi/F");
  t->Branch("dz", &m_dz, "dz/F");
  t->Branch("clusR", &m_clusR, "clusR/F");
  t->Branch("clusPhi", &m_clusPhi, "clusPhi/F");
  t->Branch("clusZ", &m_clusZ, "clusZ/F");
  t->Branch("statePhi", &m_statePhi, "statePhi/F");
  t->Branch("stateZ", &m_stateZ, "stateZ/F");
  t->Branch("stateR", &m_stateR, "stateR/F");
  t->Branch("stateRPhiErr", &m_stateRPhiErr, "stateRPhiErr/F");
  t->Branch("stateZErr", &m_stateZErr, "stateZErr/F");
  t->Branch("clusRPhiErr", &m_clusRPhiErr, "clusRPhiErr/F");
  t->Branch("clusZErr", &m_clusZErr, "clusZErr/F");
  t->Branch("cluskey", &m_cluskey, "cluskey/l");
  t->Branch("event", &m_event, "event/I");

  t->Branch("clusEta", &m_clusEta, "clusEta/F");
  t->Branch("stateEta", &m_stateEta, "stateEta/F");
  t->Branch("layer", &m_layer, "layer/I");
  t->Branch("statePt", &m_statePt, "statePt/F");
  t->Branch("statePz", &m_statePz, "statePz/F");

  t->Branch("tanAlpha_mover", &m_tanAlpha_mover, "tanAlpha_mover/F");
  t->Branch("tanBeta_mover", &m_tanBeta_mover, "tanBeta_mover/F");
  t->Branch("drphi_mover", &m_drphi_mover, "drphi_mover/F");
  t->Branch("dz_mover", &m_dz_mover, "dz_mover/F");
  t->Branch("clusR_mover", &m_clusR_mover, "clusR_mover/F");
  t->Branch("clusPhi_mover", &m_clusPhi_mover, "clusPhi_mover/F");
  t->Branch("clusZ_mover", &m_clusZ_mover, "clusZ_mover/F");
  t->Branch("statePhi_mover", &m_statePhi_mover, "statePhi_mover/F");
  t->Branch("stateZ_mover", &m_stateZ_mover, "stateZ_mover/F");
  t->Branch("stateR_mover", &m_stateR_mover, "stateR_mover/F");

  t->Branch("clusEta_mover", &m_clusEta_mover, "clusEta_mover/F");
  t->Branch("stateEta_mover", &m_stateEta_mover, "stateEta_mover/F");

  t->Branch("trackPt", &m_trackPt, "trackPt/F");
  t->Branch("charge", &m_charge, "charge/I");
  t->Branch("crossing", &m_crossing, "crossing/I");
  t->Branch("trackdEdx", &m_trackdEdx, "trackdEdx/F");
  t->Branch("track_nmvtx", &m_track_nmvtx, "track_nmvtx/I");
  t->Branch("track_nmvtxstate", &m_track_nmvtxstate, "track_nmvtxstate/I");
  t->Branch("track_nintt", &m_track_nintt, "track_nintt/I");
  t->Branch("track_ninttstate", &m_track_ninttstate, "track_ninttstate/I");
  t->Branch("track_ntpc", &m_track_ntpc, "track_ntpc/I");
  t->Branch("track_ntpcstate", &m_track_ntpcstate, "track_ntpcstate/I");
  t->Branch("track_ntpot", &m_track_ntpot, "track_ntpot/I");
  t->Branch("track_ntpotstate", &m_track_ntpotstate, "track_ntpotstate/I");

  t->Branch("cluskey_mvtx", &m_cluskey_mvtx);
  t->Branch("layer_mvtx", &m_layer_mvtx);
  t->Branch("sclusgx_mvtx", &m_sclusgx_mvtx);
  t->Branch("sclusgy_mvtx", &m_sclusgy_mvtx);
  t->Branch("sclusgz_mvtx", &m_sclusgz_mvtx);
  t->Branch("sclusgr_mvtx", &m_sclusgr_mvtx);
  t->Branch("sclusphi_mvtx", &m_sclusphi_mvtx);
  t->Branch("scluseta_mvtx", &m_scluseta_mvtx);
  t->Branch("stategx_mvtx", &m_stategx_mvtx);
  t->Branch("stategy_mvtx", &m_stategy_mvtx);
  t->Branch("stategz_mvtx", &m_stategz_mvtx);
  t->Branch("stategr_mvtx", &m_stategr_mvtx);
  t->Branch("statephi_mvtx", &m_statephi_mvtx);
  t->Branch("stateeta_mvtx", &m_stateeta_mvtx);

  t->Branch("cluskey_intt", &m_cluskey_intt);
  t->Branch("layer_intt", &m_layer_intt);
  t->Branch("sclusgx_intt", &m_sclusgx_intt);
  t->Branch("sclusgy_intt", &m_sclusgy_intt);
  t->Branch("sclusgz_intt", &m_sclusgz_intt);
  t->Branch("sclusgr_intt", &m_sclusgr_intt);
  t->Branch("sclusphi_intt", &m_sclusphi_intt);
  t->Branch("scluseta_intt", &m_scluseta_intt);
  t->Branch("stategx_intt", &m_stategx_intt);
  t->Branch("stategy_intt", &m_stategy_intt);
  t->Branch("stategz_intt", &m_stategz_intt);
  t->Branch("stategr_intt", &m_stategr_intt);
  t->Branch("statephi_intt", &m_statephi_intt);
  t->Branch("stateeta_intt", &m_stateeta_intt);

  t->Branch("cluskey_tpot", &m_cluskey_tpot);
  t->Branch("layer_tpot", &m_layer_tpot);
  t->Branch("segtype_tpot", &m_segtype_tpot);
  t->Branch("tileid_tpot", &m_tileid_tpot);
  t->Branch("sclusgx_tpot", &m_sclusgx_tpot);
  t->Branch("sclusgy_tpot", &m_sclusgy_tpot);
  t->Branch("sclusgz_tpot", &m_sclusgz_tpot);
  t->Branch("sclusgr_tpot", &m_sclusgr_tpot);
  t->Branch("sclusphi_tpot", &m_sclusphi_tpot);
  t->Branch("scluseta_tpot", &m_scluseta_tpot);
  t->Branch("stategx_tpot", &m_stategx_tpot);
  t->Branch("stategy_tpot", &m_stategy_tpot);
  t->Branch("stategz_tpot", &m_stategz_tpot);
  t->Branch("stategr_tpot", &m_stategr_tpot);
  t->Branch("statephi_tpot", &m_statephi_tpot);
  t->Branch("stateeta_tpot", &m_stateeta_tpot);

  hm->registerHisto(t);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::InitRun(PHCompositeNode* topNode)
{
  // track map
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // cluster map
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // tpc geometry
  m_tpcGeom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  // load geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  // load distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);
  if (m_disable_module_edge_corr) { m_globalPositionWrapper.set_enable_module_edge_corr(false); }
  if (m_disable_static_corr) { m_globalPositionWrapper.set_enable_static_corr(false); }
  if (m_disable_average_corr) { m_globalPositionWrapper.set_enable_average_corr(false); }
  if (m_disable_fluctuation_corr) { m_globalPositionWrapper.set_enable_fluctuation_corr(false); }

  // clusterMover needs the correct radii of the TPC layers
  m_clusterMover.initialize_geometry(m_tpcGeom);
  m_clusterMover.set_verbosity(0);

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

  auto h_beta = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "betadz"));
  assert(h_beta);

  auto h_alpha = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "alphardphi"));
  assert(h_alpha);

  auto h_rphiResid = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "rphiResid"));
  assert(h_rphiResid);

  auto h_zResid = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "zResid"));
  assert(h_zResid);

  auto h_etaResid = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "etaResid"));
  assert(h_etaResid);

  auto h_etaResidLayer = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "etaResidLayer"));
  assert(h_etaResidLayer);

  auto h_zResidLayer = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "zResidLayer"));
  assert(h_zResidLayer);

  auto h_deltarphi_layer = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "deltarphi_layer"));
  assert(h_deltarphi_layer);

  auto h_deltaz_layer = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "deltaz_layer"));
  assert(h_deltaz_layer);

  auto h_statez_pulls = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "statez_pulls"));
  assert(h_statez_pulls);

  auto h_staterphi_pulls = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "staterphi_pulls"));
  assert(h_staterphi_pulls);

  auto h_clusz_pulls = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "clusz_pulls"));
  assert(h_clusz_pulls);

  auto h_clusrphi_pulls = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "clusrphi_pulls"));
  assert(h_clusrphi_pulls);

  auto h_nstates_vs_nclus = dynamic_cast<TH2*>(hm->getHisto(get_histo_prefix() + "nstates_vs_nclus"));
  assert(h_nstates_vs_nclus);

  auto t_tree = dynamic_cast<TTree*>(hm->getHisto(get_histo_prefix() + "residTree"));
  assert(t_tree);

  std::cout << "QAG4SimulationDistortions::process_event - tracks: " << m_trackMap->size() << std::endl;
  for (const auto& [key, track] : *m_trackMap)
  {

    // total track counter
    ++m_total_tracks;

    // get track crossing and check
    const auto crossing = track->get_crossing();
    if(crossing == SHRT_MAX)
    {
      std::cout << "QAG4SimulationDistortions::process_event - invalid crossing. Track skipped." << std::endl;
      continue;
    }

    // check track quality
    if (!checkTrack(track))
    {
      if (Verbosity() > 1) { std::cout << "failed track selection" << std::endl;}
      continue;
    }
    if (Verbosity() > 1) { std::cout << "pass track selection" << std::endl;}

    // get seeeds
    auto tpcSeed = track->get_tpc_seed();
    auto siliconSeed = track->get_silicon_seed();

    /*
    std::cout<<"track crossing "<<crossing<<" , tpc_pt = "<<fabs(1. / tpcSeed->get_qOverR()) * (0.3 / 100.) * 1.4<<" , tpc_phi = "<<tpcSeed->get_phi()<<" , si_phi = "<<siliconSeed->get_phi()<<" , dphi = "<<(tpcSeed->get_phi() - siliconSeed->get_phi())<<std::endl;
    std::cout<<"tpc_eta = "<<tpcSeed->get_eta()<<" , si_eta = "<<siliconSeed->get_eta()<<" , deta = "<<(tpcSeed->get_eta() - siliconSeed->get_eta())<<std::endl;
    std::cout<<"tpc_x = "<<(TrackSeedHelper::get_xyz(tpcSeed)).x()<<" , si_x = "<<(TrackSeedHelper::get_xyz(siliconSeed)).x()<<" , dx = "<<((TrackSeedHelper::get_xyz(tpcSeed)).x() - (TrackSeedHelper::get_xyz(siliconSeed)).x())<<std::endl;
    std::cout<<"tpc_y = "<<(TrackSeedHelper::get_xyz(tpcSeed)).y()<<" , si_y = "<<(TrackSeedHelper::get_xyz(siliconSeed)).y()<<" , dy = "<<((TrackSeedHelper::get_xyz(tpcSeed)).y() - (TrackSeedHelper::get_xyz(siliconSeed)).y())<<std::endl;
    std::cout<<"tpc_z = "<<(TrackSeedHelper::get_xyz(tpcSeed)).z()<<" , si_z = "<<(TrackSeedHelper::get_xyz(siliconSeed)).z()<<" , dz = "<<((TrackSeedHelper::get_xyz(tpcSeed)).z() - (TrackSeedHelper::get_xyz(siliconSeed)).z())<<std::endl;
    */

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

    get_MvtxInttTpot_info(track);

    // accepted track counter
    ++m_accepted_tracks;

    // get the fully corrected cluster global positions
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
    for (const auto& ckey : get_cluster_keys(track))
    {
      auto cluster = m_clusterContainer->findCluster(ckey);

      // Fully correct the cluster positions for the crossing and all distortions
      Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, crossing );

      // add the global positions to a vector to give to the cluster mover
      global_raw.emplace_back(std::make_pair(ckey, global));
    }
    // move the corrected cluster positions back to the original readout surface
    auto global_moved = m_clusterMover.processTrack(global_raw);

    for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
    {
      ++m_total_states;

      auto& state = iter->second;
      const auto ckey = state->get_cluskey();
      const auto trkrId = TrkrDefs::getTrkrId(ckey);

      if( trkrId != TrkrDefs::tpcId )
      { continue; }

      ++m_accepted_states;

      auto cluster = m_clusterContainer->findCluster(ckey);

      const auto clusGlobPosition = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, crossing);

      Acts::Vector3 clusGlobPosition_moved(0, 0, 0);
      for (const auto& pair : global_moved)
      {
        auto thiskey = pair.first;
        clusGlobPosition_moved = pair.second;
        if (thiskey == ckey)
        {
          break;
        }
      }

      const float clusR = get_r(clusGlobPosition(0), clusGlobPosition(1));
      const float clusPhi = std::atan2(clusGlobPosition(1), clusGlobPosition(0));
      const float clusZ = clusGlobPosition(2);

      const float clusR_moved = get_r(clusGlobPosition_moved(0), clusGlobPosition_moved(1));
      const float clusPhi_moved = std::atan2(clusGlobPosition_moved(1), clusGlobPosition_moved(0));
      const float clusZ_moved = clusGlobPosition_moved(2);

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

      const auto stateX_unmoved = stateGlobPosition(0);
      const auto stateY_unmoved = stateGlobPosition(1);
      const auto stateZ_unmoved = stateGlobPosition(2);
      const float statePhi_unmoved = std::atan2(stateY_unmoved,stateX_unmoved);

      // Calculate residuals
      const float drphi = clusR * deltaPhi(clusPhi - statePhi);
      const float dz = clusZ - stateZ;

      const float drphi_mover = clusR_moved * deltaPhi(clusPhi_moved - statePhi_unmoved);
      const float dz_mover = clusZ_moved - stateZ_unmoved;

      const auto trackPPhi = -stateGlobMom(0) * std::sin(statePhi) + stateGlobMom(1) * std::cos(statePhi);
      const auto trackPR = stateGlobMom(0) * std::cos(statePhi) + stateGlobMom(1) * std::sin(statePhi);
      const auto trackPZ = stateGlobMom(2);

      const auto trackPPhi_mover = -stateGlobMom(0) * std::sin(statePhi_unmoved) + stateGlobMom(1) * std::cos(statePhi_unmoved);
      const auto trackPR_mover = stateGlobMom(0) * std::cos(statePhi_unmoved) + stateGlobMom(1) * std::sin(statePhi_unmoved);
      const auto trackPZ_mover = stateGlobMom(2);

      const auto trackAlpha = -trackPPhi / trackPR;
      const auto trackBeta = -trackPZ / trackPR;
      const auto trackEta = std::atanh(stateGlobMom(2) / stateGlobMom.norm());
      const auto clusEta = std::atanh(clusZ / clusGlobPosition.norm());

      const auto trackAlpha_mover = -trackPPhi_mover / trackPR_mover;
      const auto trackBeta_mover = -trackPZ_mover / trackPR_mover;
      const auto trackEta_mover = std::atanh(stateGlobMom(2) / stateGlobMom.norm());
      const auto clusEta_mover = std::atanh(clusZ_moved / clusGlobPosition_moved.norm());

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
      m_statePt = sqrt(stateGlobMom(0)*stateGlobMom(0)+stateGlobMom(1)*stateGlobMom(1));
      m_statePz = stateGlobMom(2);
      m_trackPt = track->get_pt();
      m_charge = track->get_charge();
      m_crossing = track->get_crossing();
      m_trackdEdx = calc_dedx(tpcSeed, m_clusterContainer, m_tpcGeom);

      m_tanAlpha_mover = trackAlpha_mover;
      m_tanBeta_mover = trackBeta_mover;
      m_drphi_mover = drphi_mover;
      m_dz_mover = dz_mover;
      m_clusR_mover = clusR_moved;
      m_clusPhi_mover = clusPhi_moved;
      m_clusZ_mover = clusZ_moved;
      m_statePhi_mover = statePhi_unmoved;
      m_stateZ_mover = stateZ_unmoved;
      m_stateR_mover = stateR;

      m_clusEta_mover = clusEta_mover;
      m_stateEta_mover = trackEta_mover;

      const auto cluster_keys(get_cluster_keys(track));
      m_track_nmvtx = count_clusters<TrkrDefs::mvtxId>(cluster_keys);
      m_track_nintt = count_clusters<TrkrDefs::inttId>(cluster_keys);
      m_track_ntpc = count_clusters<TrkrDefs::tpcId>(cluster_keys);
      m_track_ntpot = count_clusters<TrkrDefs::micromegasId>(cluster_keys);
      const auto state_keys(get_state_keys(track));
      m_track_nmvtxstate = count_clusters<TrkrDefs::mvtxId>(state_keys);
      m_track_ninttstate = count_clusters<TrkrDefs::inttId>(state_keys);
      m_track_ntpcstate = count_clusters<TrkrDefs::tpcId>(state_keys);
      m_track_ntpotstate = count_clusters<TrkrDefs::micromegasId>(state_keys);

      t_tree->Fill();
    }
    int nstates = track->size_states();
    int nclus = (track->get_silicon_seed()->size_cluster_keys()) + (track->get_tpc_seed()->size_cluster_keys());
    h_nstates_vs_nclus->Fill(nstates,nclus);
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

  if (track->get_pt() < 0.5)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt tpc and micromegas hits
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
//  if (count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2)
//  {
//    return false;
//  }

  const auto state_keys(get_state_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(state_keys) < 3)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(state_keys) < 2)
  {
    return false;
  }
//  if (count_clusters<TrkrDefs::micromegasId>(state_keys) < 2)
//  {
//    return false;
//  }

  if (checkTPOTResidual(track)==false)
  {
    return false;
  }

  return true;
}

bool QAG4SimulationDistortions::checkTPOTResidual(SvtxTrack* track)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h_tpot_deltarphi = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "tpot_deltarphi"));
  assert(h_tpot_deltarphi);

  auto h_tpot_deltaz = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "tpot_deltaz"));
  assert(h_tpot_deltaz);

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

    const auto cluster = m_clusterContainer->findCluster(cluskey);

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
    assert(crossing != SHRT_MAX);

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
    if (layer==55 && std::fabs(drphi)>0.1)
    {
      flag = false;
      break;
    }

    // check z residual for layer 56
    if (layer==56 && std::fabs(dz)>1)
    {
      flag = false;
      break;
    }

  }

  if (flag)
  {
    // SCOZ has a half dead tile
    // only require one TPOT cluster/state from SCOP
    if (TPOTtileID==0)
    {
      if (nTPOTcluster<1 || nTPOTstate<1)
      {
        flag = false;
      }
    }
    else if (TPOTtileID>0)
    {
      if (nTPOTcluster<2 || nTPOTstate<2)
      {
        flag = false;
      }
    }
    else if (TPOTtileID<0)
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

float QAG4SimulationDistortions::calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcGeom)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                     tpcseed->end_cluster_keys());

  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    auto detid = TrkrDefs::getTrkrId(cluster_key);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    TrkrCluster* cluster = clustermap->findCluster(cluster_key);
    float adc = cluster->getAdc();
    PHG4TpcCylinderGeom* GeoLayer_local = tpcGeom->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();
    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
    float beta = atan(tpcseed->get_slope());
    float alphacorr = cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  if (dedxlist.size() < 1)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;
}
