
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
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>

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
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
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
  t->Branch("trackPt", &m_trackPt, "trackPt/F");
  t->Branch("charge", &m_charge, "charge/I");
  t->Branch("crossing", &m_crossing, "crossing/I");

  t->Branch("segtype_tpot", &m_segtype_tpot);
  t->Branch("tileid_tpot", &m_tileid_tpot);
  t->Branch("sclusgx_tpot", &m_sclusgx_tpot);
  t->Branch("sclusgy_tpot", &m_sclusgy_tpot);
  t->Branch("sclusgz_tpot", &m_sclusgz_tpot);
  t->Branch("sclusgr_tpot", &m_sclusgr_tpot);
  t->Branch("sclusphi_tpot", &m_sclusphi_tpot);
  t->Branch("scluseta_tpot", &m_scluseta_tpot);

  hm->registerHisto(t);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::InitRun(PHCompositeNode* topNode)
{
  // track map
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxSiliconMMTrackMap");

  // cluster map
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // load geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  // load distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);


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
      continue;
    }

    // get seeeds
    auto tpcSeed = track->get_tpc_seed();
    auto siliconSeed = track->get_silicon_seed();

    /// Should have never been added to the map...
    if (!tpcSeed || !siliconSeed)
    {
      continue;
    }

    get_Tpot_info(track);

    for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
    {
      auto state = iter->second;

      /// If the state name wasn't set to the ckey, it wasn't analyzed
      /// in PHTpcResiduals (i.e. it isn't in the tpc)
      if ((state->get_name()).find("UNKNOWN") != std::string::npos)
      {
        continue;
      }

      TrkrDefs::cluskey const ckey = std::stoll(state->get_name());

      auto cluster = m_clusterContainer->findCluster(ckey);

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
      m_statePt = sqrt(stateGlobMom(0)*stateGlobMom(0)+stateGlobMom(1)*stateGlobMom(1));
      m_statePz = stateGlobMom(2);
      m_trackPt = track->get_pt();
      m_charge = track->get_charge();
      m_crossing = track->get_crossing();

      t_tree->Fill();
    }
    int nstates = track->size_states();
    int nclus = (track->get_silicon_seed()->size_cluster_keys()) + (track->get_tpc_seed()->size_cluster_keys());
    h_nstates_vs_nclus->Fill(nstates,nclus);
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool QAG4SimulationDistortions::checkTrack(SvtxTrack* track)
{

  if (track->get_pt() < 0.5)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt tpc and micromegas hits
  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(cluster_keys) < 20)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2)
  {
    return false;
  }

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

  for (const auto& cluskey : get_cluster_keys(track))
  {

    // make sure cluster is from TPOT
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if (detId != TrkrDefs::micromegasId)
    {
      continue;
    }

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
    }

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

    // check rphi and z error
    if (std::fabs(drphi)>0.1)
    {
      flag = false;
      break;
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

void QAG4SimulationDistortions::clearVector()
{
  m_segtype_tpot.clear();
  m_tileid_tpot.clear();
  m_sclusgx_tpot.clear();
  m_sclusgy_tpot.clear();
  m_sclusgz_tpot.clear();
  m_sclusgr_tpot.clear();
  m_sclusphi_tpot.clear();
  m_scluseta_tpot.clear();
}

void QAG4SimulationDistortions::get_Tpot_info(SvtxTrack* track)
{
  clearVector();
  for (const auto& ckey : get_cluster_keys(track))
  {
    if (TrkrDefs::getTrkrId(ckey)!=TrkrDefs::micromegasId) continue;
    auto cluster = m_clusterContainer->findCluster(ckey);

    int m_segtype = (int) MicromegasDefs::getSegmentationType(ckey);
    int m_tileid = MicromegasDefs::getTileId(ckey);
    Acts::Vector3 glob = m_tGeometry->getGlobalPosition(ckey, cluster);
    float m_sclusgx = glob.x();
    float m_sclusgy = glob.y();
    float m_sclusgz = glob.z();
    float m_sclusgr = get_r(m_sclusgx, m_sclusgy);
    float m_sclusphi = atan2(glob.y(), glob.x());
    float m_scluseta = acos(glob.z() / std::sqrt(square(glob.x()) + square(glob.y()) + square(glob.z())));
    m_segtype_tpot.push_back(m_segtype);
    m_tileid_tpot.push_back(m_tileid);
    m_sclusgx_tpot.push_back(m_sclusgx);
    m_sclusgy_tpot.push_back(m_sclusgy);
    m_sclusgz_tpot.push_back(m_sclusgz);
    m_sclusgr_tpot.push_back(m_sclusgr);
    m_sclusphi_tpot.push_back(m_sclusphi);
    m_scluseta_tpot.push_back(m_scluseta);
  }
}
