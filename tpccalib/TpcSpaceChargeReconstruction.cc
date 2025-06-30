/**
 * \file TpcSpaceChargeReconstruction.cc
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeReconstruction.h"
#include "TpcSpaceChargeMatrixContainerv2.h"
#include "TpcSpaceChargeMatrixContainer1D.h"
#include "TpcSpaceChargeMatrixContainer2D.h"
#include "TpcSpaceChargeReconstructionHelper.h"

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <micromegas/MicromegasDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <boost/format.hpp>

#include <cassert>
#include <memory>

namespace
{

  /// square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  /// calculate delta_phi between -pi and pi
  template <class T>
  inline constexpr T delta_phi(const T& phi)
  {
    if (phi >= M_PI)
    {
      return phi - 2 * M_PI;
    }
    else if (phi < -M_PI)
    {
      return phi + 2 * M_PI;
    }
    else
    {
      return phi;
    }
  }

  /// radius
  template <class T>
  T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  /// get sector median angle associated to a given index
  /** this assumes that sector 0 is centered on phi=0, then numbered along increasing phi */
  inline constexpr double get_sector_phi(int isec)
  {
    return isec * M_PI / 6;
  }

  // specify bins for which one will save histograms
  static const std::vector<float> phi_rec = {get_sector_phi(9)};
  static const std::vector<float> z_rec = {5.};

  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2. * M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;

  static constexpr int m_minClusCount = 10;

  static constexpr float m_layermin = 7;
  static constexpr float m_layermax = 55;

  /// get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
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

  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack* track)
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

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }

  [[maybe_unused]] std::ostream& operator<<(std::ostream& out, const Acts::Vector3& vector)
  {
    out << "(" << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

}  // namespace

//_____________________________________________________________________
TpcSpaceChargeReconstruction::TpcSpaceChargeReconstruction(const std::string& name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_matrix_container(new TpcSpaceChargeMatrixContainerv2)
  , m_matrix_container_1D_layer_negz(new TpcSpaceChargeMatrixContainer1D)
  , m_matrix_container_1D_layer_posz(new TpcSpaceChargeMatrixContainer1D)
  , m_matrix_container_1D_radius_negz(new TpcSpaceChargeMatrixContainer1D)
  , m_matrix_container_1D_radius_posz(new TpcSpaceChargeMatrixContainer1D)
  , m_matrix_container_2D_radius_z(new TpcSpaceChargeMatrixContainer2D)
{
  InitializeParameters();
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::set_grid_dimensions(int phibins, int rbins, int zbins)
{
  m_matrix_container->set_grid_dimensions(phibins, rbins, zbins);
  m_matrix_container_1D_radius_negz->set_grid_dimensions(rbins);
  m_matrix_container_1D_radius_posz->set_grid_dimensions(rbins);
  m_matrix_container_2D_radius_z->set_grid_dimensions(phibins, rbins, zbins);
}

void TpcSpaceChargeReconstruction::set_grid_dimensions(const int layerBins)
{
  m_matrix_container_1D_layer_negz->set_grid_dimensions(layerBins);
  m_matrix_container_1D_layer_posz->set_grid_dimensions(layerBins);
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::Init(PHCompositeNode* /*topNode*/)
{
  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_clusters = 0;
  m_accepted_clusters = 0;

  // histogram evaluation
  if (m_savehistograms)
  {
    create_histograms();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::InitRun(PHCompositeNode* /*unused*/)
{
  // load parameters
  UpdateParametersWithMacro();
  m_max_talpha = get_double_param("spacecharge_max_talpha");
  m_max_drphi = get_double_param("spacecharge_max_drphi");
  m_max_tbeta = get_double_param("spacecharge_max_tbeta");
  m_max_dz = get_double_param("spacecharge_max_dz");

  // print
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_outputfile: " << m_outputfile << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_use_micromegas: " << std::boolalpha << m_use_micromegas << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_talpha: " << m_max_talpha << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_drphi: " << m_max_drphi << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_tbeta: " << m_max_tbeta << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_max_dz: " << m_max_dz << std::endl;
  std::cout << "TpcSpaceChargeReconstruction::InitRun - m_min_pt: " << m_min_pt << " GeV/c" << std::endl;

  // also identify the matrix container
  m_matrix_container->identify();
  m_matrix_container_1D_radius_negz->identify();
  m_matrix_container_1D_radius_posz->identify();
  m_matrix_container_2D_radius_z->identify();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::End(PHCompositeNode* /*topNode*/)
{
  // save matrix container in output file
  if (m_matrix_container && m_matrix_container_1D_layer_negz && m_matrix_container_1D_layer_posz && m_matrix_container_1D_radius_negz && m_matrix_container_1D_radius_posz && m_matrix_container_2D_radius_z)
  {
    std::unique_ptr<TFile> outputfile(TFile::Open(m_outputfile.c_str(), "RECREATE"));
    outputfile->cd();
    m_matrix_container->Write("TpcSpaceChargeMatrixContainer");
    m_matrix_container_1D_layer_negz->Write("TpcSpaceChargeMatrixContainer_1D_layer_negz");
    m_matrix_container_1D_layer_posz->Write("TpcSpaceChargeMatrixContainer_1D_layer_posz");
    m_matrix_container_1D_radius_negz->Write("TpcSpaceChargeMatrixContainer_1D_radius_negz");
    m_matrix_container_1D_radius_posz->Write("TpcSpaceChargeMatrixContainer_1D_radius_posz");
    m_matrix_container_2D_radius_z->Write("TpcSpaceChargeMatrixContainer_2D_radius_z");
  }

  // save histograms
  if (m_savehistograms && m_histogramfile)
  {
    m_histogramfile->cd();
    for (const auto& [cell, h] : m_h_drphi)
    {
      if (h)
      {
        h->Write();
      }
    }
    for (const auto& [cell, h] : m_h_dz)
    {
      if (h)
      {
        h->Write();
      }
    }
    for (const auto& [cell, h] : m_h_drphi_alpha)
    {
      if (h)
      {
        h->Write();
      }
    }
    for (const auto& [cell, h] : m_h_dz_beta)
    {
      if (h)
      {
        h->Write();
      }
    }
    m_histogramfile->Close();
  }

  // print counters
  std::cout
      << "TpcSpaceChargeReconstruction::End -"
      << " track statistics total: " << m_total_tracks
      << " accepted: " << m_accepted_tracks
      << " fraction: " << 100. * m_accepted_tracks / m_total_tracks << "%"
      << std::endl;

  std::cout
      << "TpcSpaceChargeReconstruction::End -"
      << " cluster statistics total: " << m_total_clusters
      << " accepted: " << m_accepted_clusters << " fraction: "
      << 100. * m_accepted_clusters / m_total_clusters << "%"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void TpcSpaceChargeReconstruction::SetDefaultParameters()
{
  // residual cuts
  set_default_double_param("spacecharge_max_talpha", 0.6);
  set_default_double_param("spacecharge_max_drphi", 0.5);
  set_default_double_param("spacecharge_max_tbeta", 1.5);
  set_default_double_param("spacecharge_max_dz", 0.5);
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::load_nodes(PHCompositeNode* topNode)
{

  // tracks
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, _svtx_track_map_name);
  assert(m_track_map);

  // clusters
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if(!m_cluster_map)
  { m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); }

  assert(m_cluster_map);

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);
  if (m_disable_module_edge_corr) { m_globalPositionWrapper.set_enable_module_edge_corr(false); }
  if (m_disable_static_corr) { m_globalPositionWrapper.set_enable_static_corr(false); }
  if (m_disable_average_corr) { m_globalPositionWrapper.set_enable_average_corr(false); }
  if (m_disable_fluctuation_corr) { m_globalPositionWrapper.set_enable_fluctuation_corr(false); }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::create_histograms()
{
  std::cout << "TpcSpaceChargeReconstruction::create_histograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset(new TFile(m_histogramfilename.c_str(), "RECREATE"));
  m_histogramfile->cd();

  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions(phibins, rbins, zbins);

  // get bins corresponding to selected angles
  std::set<int> phibin_rec;
  std::transform(phi_rec.begin(), phi_rec.end(), std::inserter(phibin_rec, phibin_rec.end()), [&](const float& phi)
                 { return phibins * (phi - m_phimin) / (m_phimax - m_phimin); });

  std::set<int> zbin_rec;
  std::transform(z_rec.begin(), z_rec.end(), std::inserter(zbin_rec, zbin_rec.end()), [&](const float& z)
                 { return zbins * (z - m_zmin) / (m_zmax - m_zmin); });

  // keep track of all cell ids that match selected histograms
  for (int iphi = 0; iphi < phibins; ++iphi)
  {
    for (int ir = 0; ir < rbins; ++ir)
    {
      for (int iz = 0; iz < zbins; ++iz)
      {
        if (phibin_rec.find(iphi) == phibin_rec.end() || zbin_rec.find(iz) == zbin_rec.end())
        {
          continue;
        }
        const auto icell = m_matrix_container->get_cell_index(iphi, ir, iz);

        {
          // rphi residuals
	  std::string hname = (boost::format("residual_drphi_p%i_r%i_z%i") %iphi %ir %iz).str();
          auto h = new TH1F(hname.c_str(), hname.c_str(), 100, -m_max_drphi, +m_max_drphi);
          h->GetXaxis()->SetTitle("r.#Delta#phi_{cluster-track} (cm)");
          m_h_drphi.insert(std::make_pair(icell, h));
        }

        {
          // 2D histograms
	  std::string hname = (boost::format("residual_2d_drphi_p%i_r%i_z%i") %iphi %ir %iz).str();
          auto h = new TH2F(hname.c_str(), hname.c_str(), 100, -m_max_talpha, m_max_talpha, 100, -m_max_drphi, +m_max_drphi);
          h->GetXaxis()->SetTitle("tan#alpha");
          h->GetYaxis()->SetTitle("r.#Delta#phi_{cluster-track} (cm)");
          m_h_drphi_alpha.insert(std::make_pair(icell, h));
        }

        {
          // z residuals
	  std::string hname = (boost::format("residual_dz_p%i_r%i_z%i") %iphi %ir %iz).str();
          auto h = new TH1F(hname.c_str(), hname.c_str(), 100, -m_max_dz, +m_max_dz);
          h->GetXaxis()->SetTitle("#Deltaz_{cluster-track} (cm)");
          m_h_dz.insert(std::make_pair(icell, h));
        }

        {
          // 2D histograms
          static constexpr double max_tbeta = 0.5;
	  std::string hname = (boost::format("residual_2d_dz_p%i_r%i_z%i") %iphi %ir %iz).str();
          auto h = new TH2F(hname.c_str(), hname.c_str(), 100, -max_tbeta, max_tbeta, 100, -m_max_dz, +m_max_dz);
          h->GetXaxis()->SetTitle("tan#beta");
          h->GetYaxis()->SetTitle("#Deltaz_{cluster-track} (cm)");
          m_h_dz_beta.insert(std::make_pair(icell, h));
        }
      }
    }
  }
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_tracks()
{
  if (!(m_track_map && m_cluster_map))
  {
    return;
  }

  for (const auto& [trackKey, track] : *m_track_map)
  {
    ++m_total_tracks;
    if (accept_track(track))
    {
      ++m_accepted_tracks;
      process_track(track);
    }
  }
}

//_____________________________________________________________________
bool TpcSpaceChargeReconstruction::accept_track(SvtxTrack* track) const
{
  if (m_requireCrossing && track->get_crossing()!=0)
  {
    return false;
  }

  if (track->get_pt() < m_min_pt)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt and micromegas hits
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
//  if (m_use_micromegas && count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2)
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
//  if (m_use_micromegas && count_clusters<TrkrDefs::micromegasId>(state_keys) < 2)
//  {
//    return false;
//  }

  if (m_use_micromegas && checkTPOTResidual(track)==false)
  {
    return false;
  }

  // all tests passed
  return true;
}

//___________________________________________________________________________________
bool TpcSpaceChargeReconstruction::checkTPOTResidual(SvtxTrack* track) const
{
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

    const auto cluster = m_cluster_map->findCluster(cluskey);

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
    const double drp = clusR * delta_phi(clusPhi - trackPhi);
    if (std::isnan(drp))
    {
      continue;
    }

    const double dz = clusZ - trackZ;
    if (std::isnan(dz))
    {
      continue;
    }

    if (Verbosity() > 3)
    {
      std::cout << "TpcSpaceChargeReconstruction::checkTPOTResidual -"
                << " drp: " << drp
                << " dz: " << dz
                << std::endl;
    }

    // check rphi residual for layer 55
    if (layer==55 && std::fabs(drp)>0.1)
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

//___________________________________________________________________________________
bool TpcSpaceChargeReconstruction::checkTrackCM(SvtxTrack* track) const
{
  if (Verbosity() > 2)
  {
    std::cout << "TpcSpaceChargeReconstruction::checkTrackCM - pcaz: " << track->get_z() << std::endl;
  }

  if (m_requireCM && (fabs(track->get_z()) > m_pcazcut || fabs(track->get_eta()) > m_etacut))
  {
    return false;
  }

  return true;
}

//_____________________________________________________________________
void TpcSpaceChargeReconstruction::process_track(SvtxTrack* track)
{
  // printout all track state
  if (Verbosity())
  {
    for (auto&& iter = track->begin_states(); iter != track->end_states(); ++iter)
    {
      const auto& [pathlength, state] = *iter;
      const auto r = std::sqrt(square(state->get_x()) + square(state->get_y()));
      const auto phi = std::atan2(state->get_y(), state->get_x());
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
                << " pathlength: " << pathlength
                << " radius: " << r
                << " phi: " << phi
                << " z: " << state->get_z()
                << std::endl;
    }
  }

  // store crossing
  // it is needed for geting cluster's global position
  const auto crossing = track->get_crossing();
  assert(crossing != SHRT_MAX);

  bool nearCM = checkTrackCM(track);

  // running track state
  // loop over clusters
  for (const auto& cluster_key : get_cluster_keys(track))
  {
    auto cluster = m_cluster_map->findCluster(cluster_key);
    if (!cluster)
    {
      std::cout << PHWHERE << " unable to find cluster for key " << cluster_key << std::endl;
      continue;
    }

    ++m_total_clusters;

    // make sure
    const auto detId = TrkrDefs::getTrkrId(cluster_key);
    if (detId != TrkrDefs::tpcId)
    {
      continue;
    }

    // cluster r, phi and z
    const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
    const auto cluster_r = get_r(global_position.x(), global_position.y());
    const auto cluster_phi = std::atan2(global_position.y(), global_position.x());
    const auto cluster_z = global_position.z();

    // cluster errors
    double cluster_rphi_error = cluster->getRPhiError();
    double cluster_z_error = cluster->getZError();

    // find track state that match cluster
    bool found = false;
    auto state_iter = track->begin_states();
    for(; state_iter != track->end_states(); ++state_iter)
    {
      const auto& [pathlengh, state] = *state_iter;
      if( state->get_cluskey() == cluster_key )
      {
        found = true;
        break;
      }
    }

    if( !found )
    {
      if( Verbosity() )
      { std::cout << "TpcSpaceChargeReconstruction::process_track - could not find track state for layer: " <<  (int) TrkrDefs::getLayer(cluster_key) << std::endl; }
      continue;
    }

    // get relevant track state
    const auto state = state_iter->second;

    // extrapolate track parameters to the cluster r
    const auto track_r = get_r(state->get_x(), state->get_y());
    const auto dr = cluster_r - track_r;
    const auto track_drdt = (state->get_x() * state->get_px() + state->get_y() * state->get_py()) / track_r;
    const auto track_dxdr = state->get_px() / track_drdt;
    const auto track_dydr = state->get_py() / track_drdt;
    const auto track_dzdr = state->get_pz() / track_drdt;

    // store state position
    const auto track_x = state->get_x() + dr * track_dxdr;
    const auto track_y = state->get_y() + dr * track_dydr;
    const auto track_z = state->get_z() + dr * track_dzdr;
    const auto track_phi = std::atan2(track_y, track_x);

    // get track angles
    const auto cosphi(std::cos(track_phi));
    const auto sinphi(std::sin(track_phi));
    const auto track_pphi = -state->get_px() * sinphi + state->get_py() * cosphi;
    const auto track_pr = state->get_px() * cosphi + state->get_py() * sinphi;
    const auto track_pz = state->get_pz();
    const auto talpha = -track_pphi / track_pr;
    const auto tbeta = -track_pz / track_pr;

    // sanity check
    if (std::isnan(talpha))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - talpha is nan" << std::endl;
      continue;
    }

    if (std::isnan(tbeta))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - tbeta is nan" << std::endl;
      continue;
    }

    // track errors
    const auto track_rphi_error = state->get_rphi_error();
    const auto track_z_error = state->get_z_error();

    // residuals
    const auto drp = cluster_r * delta_phi(cluster_phi - track_phi);
    const auto dz = cluster_z - track_z;

    // sanity checks
    if (std::isnan(drp))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - drp is nan" << std::endl;
      continue;
    }

    if (std::isnan(dz))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - dz is nan" << std::endl;
      continue;
    }

    // residual errors squared
    const auto erp = square(track_rphi_error) + square(cluster_rphi_error);
    const auto ez = square(track_z_error) + square(cluster_z_error);

    // sanity check
    // TODO: check whether this happens and fix upstream
    if (std::isnan(erp))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - erp is nan" << std::endl;
      continue;
    }

    if (std::isnan(ez))
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - ez is nan" << std::endl;
      continue;
    }

    // get cell
    const auto i = get_cell_index(global_position);
    const auto i_1D_layer = get_cell_index_layer((int) TrkrDefs::getLayer(cluster_key));
    const auto i_1D_radius = get_cell_index_radius(global_position);
    const auto i_2D_rz = get_cell_index_rz(global_position);
    if (Verbosity() > 3)
    {
      std::cout << "Bin index found is " << i << std::endl;
      std::cout << "Bin index (1D layer) found is " << i_1D_layer << std::endl;
      std::cout << "Bin index (1D radius) found is " << i_1D_radius << std::endl;
      std::cout << "Bin index (2D radius & z) found is " << i_2D_rz << std::endl;
    }
    if (i < 0 || i_1D_layer < 0 || i_1D_radius < 0 || i_2D_rz < 0)
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - invalid cell index" << std::endl;
      continue;
    }

    if (m_savehistograms)
    {
      {
        const auto iter = m_h_drphi.find(i);
        if (iter != m_h_drphi.end())
        {
          iter->second->Fill(drp);
        }
      }
      {
        const auto iter = m_h_drphi_alpha.find(i);
        if (iter != m_h_drphi_alpha.end())
        {
          iter->second->Fill(talpha, drp);
        }
      }
      {
        const auto iter = m_h_dz.find(i);
        if (iter != m_h_dz.end())
        {
          iter->second->Fill(dz);
        }
      }
      {
        const auto iter = m_h_dz_beta.find(i);
        if (iter != m_h_dz_beta.end())
        {
          iter->second->Fill(tbeta, dz);
        }
      }
    }

    // check against limits
    if (std::abs(talpha) > m_max_talpha)
    {
      continue;
    }
    if (std::abs(tbeta) > m_max_tbeta)
    {
      continue;
    }

    // check against limits
    if (std::abs(drp) > m_max_drphi)
    {
      continue;
    }
    if (std::abs(dz) > m_max_dz)
    {
      continue;
    }

    // check rphi and z error
    if (std::sqrt(erp) < m_minRPhiErr)
    {
      continue;
    }

    if (std::sqrt(ez) < m_minZErr)
    {
      continue;
    }

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeReconstruction::process_track - layer: " << (int) TrkrDefs::getLayer(cluster_key) << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
                << " cluster: (" << cluster_r << ", " << cluster_r * cluster_phi << ", " << cluster_z << ")"
                << " (" << cluster_rphi_error << ", " << cluster_z_error << ")"
                << std::endl;
      std::cout << "TpcSpaceChargeReconstruction::process_track -"
                << " track: (" << track_r << ", " << cluster_r * track_phi << ", " << track_z << ")"
                << " (" << talpha << ", " << tbeta << ")"
                << " (" << track_rphi_error << ", " << track_z_error << ")"
                << std::endl;
      std::cout << std::endl;
    }

    if (nearCM)
    {
      // Fill distortion matrices
      if (cluster_z<0)
      {
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 0, 0, square(cluster_r) / erp);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 0, 1, 0);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 0, 2, cluster_r*talpha / erp);

        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 1, 0, 0);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 1, 1, 1. / ez);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 1, 2, tbeta / ez);

        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 2, 0, cluster_r * talpha / erp);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 2, 1, tbeta / ez);
        m_matrix_container_1D_layer_negz->add_to_lhs(i_1D_layer, 2, 2, square(talpha) / erp + square(tbeta) / ez);

        m_matrix_container_1D_layer_negz->add_to_rhs(i_1D_layer, 0, cluster_r*drp / erp);
        m_matrix_container_1D_layer_negz->add_to_rhs(i_1D_layer, 1, dz / ez);
        m_matrix_container_1D_layer_negz->add_to_rhs(i_1D_layer, 2, talpha * drp / erp + tbeta * dz / ez);

        // update entries in cell
        m_matrix_container_1D_layer_negz->add_to_entries(i_1D_layer);


        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 0, 0, square(cluster_r) / erp);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 0, 1, 0);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 0, 2, cluster_r*talpha / erp);

        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 1, 0, 0);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 1, 1, 1. / ez);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 1, 2, tbeta / ez);

        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 2, 0, cluster_r * talpha / erp);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 2, 1, tbeta / ez);
        m_matrix_container_1D_radius_negz->add_to_lhs(i_1D_radius, 2, 2, square(talpha) / erp + square(tbeta) / ez);

        m_matrix_container_1D_radius_negz->add_to_rhs(i_1D_radius, 0, cluster_r*drp / erp);
        m_matrix_container_1D_radius_negz->add_to_rhs(i_1D_radius, 1, dz / ez);
        m_matrix_container_1D_radius_negz->add_to_rhs(i_1D_radius, 2, talpha * drp / erp + tbeta * dz / ez);

        // update entries in cell
        m_matrix_container_1D_radius_negz->add_to_entries(i_1D_radius);
      }
      else if (cluster_z>=0)
      {
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 0, 0, square(cluster_r) / erp);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 0, 1, 0);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 0, 2, cluster_r*talpha / erp);

        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 1, 0, 0);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 1, 1, 1. / ez);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 1, 2, tbeta / ez);

        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 2, 0, cluster_r * talpha / erp);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 2, 1, tbeta / ez);
        m_matrix_container_1D_layer_posz->add_to_lhs(i_1D_layer, 2, 2, square(talpha) / erp + square(tbeta) / ez);

        m_matrix_container_1D_layer_posz->add_to_rhs(i_1D_layer, 0, cluster_r*drp / erp);
        m_matrix_container_1D_layer_posz->add_to_rhs(i_1D_layer, 1, dz / ez);
        m_matrix_container_1D_layer_posz->add_to_rhs(i_1D_layer, 2, talpha * drp / erp + tbeta * dz / ez);

        // update entries in cell
        m_matrix_container_1D_layer_posz->add_to_entries(i_1D_layer);


        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 0, 0, square(cluster_r) / erp);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 0, 1, 0);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 0, 2, cluster_r*talpha / erp);

        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 1, 0, 0);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 1, 1, 1. / ez);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 1, 2, tbeta / ez);

        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 2, 0, cluster_r * talpha / erp);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 2, 1, tbeta / ez);
        m_matrix_container_1D_radius_posz->add_to_lhs(i_1D_radius, 2, 2, square(talpha) / erp + square(tbeta) / ez);

        m_matrix_container_1D_radius_posz->add_to_rhs(i_1D_radius, 0, cluster_r*drp / erp);
        m_matrix_container_1D_radius_posz->add_to_rhs(i_1D_radius, 1, dz / ez);
        m_matrix_container_1D_radius_posz->add_to_rhs(i_1D_radius, 2, talpha * drp / erp + tbeta * dz / ez);

        // update entries in cell
        m_matrix_container_1D_radius_posz->add_to_entries(i_1D_radius);
      }
    }

    // Fill distortion matrices
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 0, 0, square(cluster_r) / erp);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 0, 1, 0);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 0, 2, cluster_r*talpha / erp);

    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 1, 0, 0);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 1, 1, 1. / ez);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 1, 2, tbeta / ez);

    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 2, 0, cluster_r * talpha / erp);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 2, 1, tbeta / ez);
    m_matrix_container_2D_radius_z->add_to_lhs(i_2D_rz, 2, 2, square(talpha) / erp + square(tbeta) / ez);

    m_matrix_container_2D_radius_z->add_to_rhs(i_2D_rz, 0, cluster_r*drp / erp);
    m_matrix_container_2D_radius_z->add_to_rhs(i_2D_rz, 1, dz / ez);
    m_matrix_container_2D_radius_z->add_to_rhs(i_2D_rz, 2, talpha * drp / erp + tbeta * dz / ez);

    // update entries in cell
    m_matrix_container_2D_radius_z->add_to_entries(i_2D_rz);


    // update matrices
    // see https://indico.bnl.gov/event/7440/contributions/43328/attachments/31334/49446/talk.pdf for details
    m_matrix_container->add_to_lhs(i, 0, 0, square(cluster_r) / erp);
    m_matrix_container->add_to_lhs(i, 0, 1, 0);
    m_matrix_container->add_to_lhs(i, 0, 2, cluster_r*talpha / erp);

    m_matrix_container->add_to_lhs(i, 1, 0, 0);
    m_matrix_container->add_to_lhs(i, 1, 1, 1. / ez);
    m_matrix_container->add_to_lhs(i, 1, 2, tbeta / ez);

    m_matrix_container->add_to_lhs(i, 2, 0, cluster_r * talpha / erp);
    m_matrix_container->add_to_lhs(i, 2, 1, tbeta / ez);
    m_matrix_container->add_to_lhs(i, 2, 2, square(talpha) / erp + square(tbeta) / ez);

    m_matrix_container->add_to_rhs(i, 0, cluster_r*drp / erp);
    m_matrix_container->add_to_rhs(i, 1, dz / ez);
    m_matrix_container->add_to_rhs(i, 2, talpha * drp / erp + tbeta * dz / ez);

    // also update rphi reduced matrices
    m_matrix_container->add_to_lhs_rphi(i, 0, 0, square(cluster_r) / erp);
    m_matrix_container->add_to_lhs_rphi(i, 0, 1, cluster_r*talpha / erp);
    m_matrix_container->add_to_lhs_rphi(i, 1, 0, cluster_r * talpha / erp);
    m_matrix_container->add_to_lhs_rphi(i, 1, 1, square(talpha) / erp);

    m_matrix_container->add_to_rhs_rphi(i, 0, cluster_r*drp / erp);
    m_matrix_container->add_to_rhs_rphi(i, 1, talpha * drp / erp);

    // also update z reduced matrices
    m_matrix_container->add_to_lhs_z(i, 0, 0, 1. / ez);
    m_matrix_container->add_to_lhs_z(i, 0, 1, tbeta / ez);
    m_matrix_container->add_to_lhs_z(i, 1, 0, tbeta / ez);
    m_matrix_container->add_to_lhs_z(i, 1, 1, square(tbeta) / ez);

    m_matrix_container->add_to_rhs_z(i, 0, dz / ez);
    m_matrix_container->add_to_rhs_z(i, 1, tbeta * dz / ez);

    // update entries in cell
    m_matrix_container->add_to_entries(i);

    // increment number of accepted clusters
    ++m_accepted_clusters;
  }
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell_index(const Acts::Vector3& global_position)
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions(phibins, rbins, zbins);

  // phi
  // bound check
  float phi = std::atan2(global_position.y(), global_position.x());
  while (phi < m_phimin)
  {
    phi += 2. * M_PI;
  }
  while (phi >= m_phimax)
  {
    phi -= 2. * M_PI;
  }
  int iphi = phibins * (phi - m_phimin) / (m_phimax - m_phimin);

  // radius
  const float r = get_r(global_position.x(), global_position.y());
  if (r < m_rmin || r >= m_rmax)
  {
    return -1;
  }
  int ir = rbins * (r - m_rmin) / (m_rmax - m_rmin);

  // z
  const float z = global_position.z();
  if (z < m_zmin || z >= m_zmax)
  {
    return -1;
  }
  int iz = zbins * (z - m_zmin) / (m_zmax - m_zmin);

  return m_matrix_container->get_cell_index(iphi, ir, iz);
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell_index_layer(const int layer)
{
  // get grid dimensions from matrix container
  int layerbins = 0;
  m_matrix_container_1D_layer_negz->get_grid_dimensions(layerbins);

  // layer
  if (layer < m_layermin || layer >= m_layermax)
  {
    return -1;
  }
  const int ilayer = layerbins * (layer - m_layermin) / (m_layermax - m_layermin);

  // get index from matrix container
  return m_matrix_container_1D_layer_negz->get_cell_index(ilayer);
}

//_____________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell_index_radius(const Acts::Vector3& loc)
{
  // get grid dimensions from matrix container
  int rbins = 0;
  m_matrix_container_1D_radius_negz->get_grid_dimensions(rbins);

  // r
  const auto r = get_r(loc(0), loc(1));
  if (r < m_rmin || r >= m_rmax)
  {
    return -1;
  }
  const int ir = rbins * (r - m_rmin) / (m_rmax - m_rmin);


  // get index from matrix container
  return m_matrix_container_1D_radius_negz->get_cell_index(ir);
}

//_______________________________________________________________________________
int TpcSpaceChargeReconstruction::get_cell_index_rz(const Acts::Vector3& loc)
{
  // get grid dimensions from matrix container
  int pbins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container_2D_radius_z->get_grid_dimensions(pbins, rbins, zbins);

  // r
  const auto r = get_r(loc(0), loc(1));
  if (r < m_rmin || r >= m_rmax)
  {
    return -1;
  }
  const int ir = rbins * (r - m_rmin) / (m_rmax - m_rmin);

  // z
  const auto z = loc(2);
  if (z < m_zmin || z >= m_zmax)
  {
    return -1;
  }
  const int iz = zbins * (z - m_zmin) / (m_zmax - m_zmin);

  // get index from matrix container
  return m_matrix_container_2D_radius_z->get_cell_index(ir, iz);
}
