// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONDISTORTIONS_H
#define QAG4SIMULATIONDISTORTIONS_H

#include <fun4all/SubsysReco.h>
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

 private:
  std::string get_histo_prefix()
  {
    return std::string("h_") + Name() + std::string("_");
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  bool checkTrack(SvtxTrack* track);
  bool checkTPOTResidual(SvtxTrack* track);
  float calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcGeom);
  void clearVector();
  void get_Tpot_info(SvtxTrack* track);
  SvtxTrackMap* m_trackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TpcCylinderGeomContainer *m_tpcGeom = nullptr;

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
  int m_layer = -1;
  float m_statePt = NAN;
  float m_statePz = NAN;
  float m_trackPt = NAN;
  float m_trackdEdx = NAN;
  int m_charge = -10;
  int m_crossing = -10;

  std::vector<int> m_segtype_tpot;
  std::vector<int> m_tileid_tpot;
  std::vector<float> m_sclusgx_tpot;
  std::vector<float> m_sclusgy_tpot;
  std::vector<float> m_sclusgz_tpot;
  std::vector<float> m_sclusgr_tpot;
  std::vector<float> m_sclusphi_tpot;
  std::vector<float> m_scluseta_tpot;

  TrkrDefs::cluskey m_cluskey = TrkrDefs::CLUSKEYMAX;
};

#endif  // QAG4SIMULATIONDISTORTIONS_H
