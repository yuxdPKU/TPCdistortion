// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONDISTORTIONS_H
#define QAG4SIMULATIONDISTORTIONS_H

#include <fun4all/SubsysReco.h>
#include <tpc/TpcClusterMover.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeed.h>

#include <math.h>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class TrkrClusterContainer;
class SvtxTrack;
class ActsGeometry;
class PHG4TpcCylinderGeomContainer;

class QAG4SimulationDistortions : public SubsysReco
{
 public:
  QAG4SimulationDistortions(const std::string& name = "QAG4SimulationDistortions");

  ~QAG4SimulationDistortions() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  //! track map name
  void set_trackmap_name( const std::string& value )
  { m_trackmapname = value; }

  void disableModuleEdgeCorr() { m_disable_module_edge_corr = true; }
  void disableStaticCorr() { m_disable_static_corr = true; }
  void disableAverageCorr() { m_disable_average_corr = true; }
  void disableFluctuationCorr() { m_disable_fluctuation_corr = true; }

 private:

  //! track map name
  std::string m_trackmapname = "SvtxSiliconMMTrackMap";

  std::string get_histo_prefix()
  {
    return std::string("h_") + Name() + std::string("_");
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack* track);
  bool checkTrack(SvtxTrack* track);
  bool checkTPOTResidual(SvtxTrack* track);
  float calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcGeom);
  void clearVector();
  void get_MvtxInttTpot_info(SvtxTrack* track);
  SvtxTrackMap* m_trackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TpcCylinderGeomContainer *m_tpcGeom = nullptr;

  TpcClusterMover m_clusterMover;

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  int m_event = 0;
  float m_tanAlpha = NAN;
  float m_tanBeta = NAN;
  float m_drphi = NAN;
  float m_dz = NAN;
  float m_clusR = NAN;
  float m_clusPhi = NAN;
  float m_clusZ = NAN;
  float m_statePhi = NAN;
  float m_stateZ = NAN;
  float m_stateR = NAN;
  float m_stateRPhiErr = NAN;
  float m_stateZErr = NAN;
  float m_clusRPhiErr = NAN;
  float m_clusZErr = NAN;

  float m_clusEta = NAN;
  float m_stateEta = NAN;

  float m_tanAlpha_mover = NAN;
  float m_tanBeta_mover = NAN;
  float m_drphi_mover = NAN;
  float m_dz_mover = NAN;
  float m_clusR_mover = NAN;
  float m_clusPhi_mover = NAN;
  float m_clusZ_mover = NAN;
  float m_statePhi_mover = NAN;
  float m_stateZ_mover = NAN;
  float m_stateR_mover = NAN;

  float m_clusEta_mover = NAN;
  float m_stateEta_mover = NAN;

  int m_layer = -1;
  float m_statePt = NAN;
  float m_statePz = NAN;
  float m_trackPt = NAN;
  float m_trackdEdx = NAN;
  int m_track_nmvtx = 0;
  int m_track_nmvtxstate = 0;
  int m_track_nintt = 0;
  int m_track_ninttstate = 0;
  int m_track_ntpc = 0;
  int m_track_ntpcstate = 0;
  int m_track_ntpot = 0;
  int m_track_ntpotstate = 0;
  int m_charge = -10;
  int m_crossing = -10;

  std::vector<TrkrDefs::cluskey> m_cluskey_mvtx;
  std::vector<int> m_layer_mvtx;
  std::vector<float> m_sclusgx_mvtx;
  std::vector<float> m_sclusgy_mvtx;
  std::vector<float> m_sclusgz_mvtx;
  std::vector<float> m_sclusgr_mvtx;
  std::vector<float> m_sclusphi_mvtx;
  std::vector<float> m_scluseta_mvtx;
  std::vector<float> m_stategx_mvtx;
  std::vector<float> m_stategy_mvtx;
  std::vector<float> m_stategz_mvtx;
  std::vector<float> m_stategr_mvtx;
  std::vector<float> m_statephi_mvtx;
  std::vector<float> m_stateeta_mvtx;

  std::vector<TrkrDefs::cluskey> m_cluskey_intt;
  std::vector<int> m_layer_intt;
  std::vector<float> m_sclusgx_intt;
  std::vector<float> m_sclusgy_intt;
  std::vector<float> m_sclusgz_intt;
  std::vector<float> m_sclusgr_intt;
  std::vector<float> m_sclusphi_intt;
  std::vector<float> m_scluseta_intt;
  std::vector<float> m_stategx_intt;
  std::vector<float> m_stategy_intt;
  std::vector<float> m_stategz_intt;
  std::vector<float> m_stategr_intt;
  std::vector<float> m_statephi_intt;
  std::vector<float> m_stateeta_intt;

  std::vector<TrkrDefs::cluskey> m_cluskey_tpot;
  std::vector<int> m_layer_tpot;
  std::vector<int> m_segtype_tpot;
  std::vector<int> m_tileid_tpot;
  std::vector<float> m_sclusgx_tpot;
  std::vector<float> m_sclusgy_tpot;
  std::vector<float> m_sclusgz_tpot;
  std::vector<float> m_sclusgr_tpot;
  std::vector<float> m_sclusphi_tpot;
  std::vector<float> m_scluseta_tpot;
  std::vector<float> m_stategx_tpot;
  std::vector<float> m_stategy_tpot;
  std::vector<float> m_stategz_tpot;
  std::vector<float> m_stategr_tpot;
  std::vector<float> m_statephi_tpot;
  std::vector<float> m_stateeta_tpot;

  TrkrDefs::cluskey m_cluskey = TrkrDefs::CLUSKEYMAX;

  /// disable distortion correction
  bool m_disable_module_edge_corr = false;
  bool m_disable_static_corr = false;
  bool m_disable_average_corr = false;
  bool m_disable_fluctuation_corr = false;

  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_states = 0;
  int m_accepted_states = 0;
  //@}


};

#endif  // QAG4SIMULATIONDISTORTIONS_H
