#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
/**
 * \file TpcSpaceChargeReconstruction.h
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <tpc/TpcGlobalPositionWrapper.h>

/// Acts includes to create all necessary definitions
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <vector>

// forward declaration
class SvtxTrack;
class SvtxTrackMap;
class TpcSpaceChargeMatrixContainer;
class TrkrCluster;
class TrkrClusterContainer;
class ActGeometry;

class TFile;
class TH1;
class TH2;

/**
 * \class TpcSpaceChargeReconstruction
 * \brief performs space charge distortion reconstruction using tracks
 * \detail To reconstruct the distortions dr0, drphi0 and dz0 in a given volume element, the following chisquare is minimized
 chisquare = sum_cluster (drphi - (drphi0 + dr0 tan alpha))**2/error**2 + sum_cluster ( dz - (dz0 + dr0 tan beta))**2/error**2
 with
 - drphi and dz the residuals (track - cluster) measured for a given cluster
 - alpha and beta the track angles at the cluster in the rphi,r plane and the z,r plane, respectively
 The chisquare being quadratic in dr0, drphi0 and dz0, it can be minimized analytically.
 This results in a linear equation lhs[i].[corrections] = rhs[i], and thus [corrections] = lhs[i]**(-1).rhs[i]
 The lhs and rhs matrices are filled in TpcSpaceChargeReconstruction::process_track
 The actual inversion is performed in TpcSpaceChargeMatrixInversion::calculate_distortions
 */

class TpcSpaceChargeReconstruction : public SubsysReco, public PHParameterInterface
{
 public:
  /// constructor
  TpcSpaceChargeReconstruction(const std::string& = "TPCSPACECHARGERECONSTRUCTION");

  ///@name configuration
  //@{

  /// set whether to use only tracks with micromegas or not
  void set_use_micromegas(bool value)
  {
    m_use_micromegas = value;
  }

  void setMinRPhiErr(float minRPhiErr)
  {
    m_minRPhiErr = minRPhiErr;
  }

  void setMinZErr(float minZErr)
  {
    m_minZErr = minZErr;
  }

  /// track min pT
  void set_min_pt(double value)
  {
    m_min_pt = value;
  }

  /// track pcaz cut
  void setPCAzcut(double value)
  {
    m_pcazcut = value;
  }

  /// track eta cut
  void setEtacut(double value)
  {
    m_etacut = value;
  }

  /// track crossing
  void requireCrossing(bool flag = true)
  {
    m_requireCrossing = flag;
  }

  /// track near CM
  void requireCM(bool flag = true)
  {
    m_requireCM = flag;
  }

  /// set grid dimensions
  /**
  \param phibins the number of bins in the azimuth direction
  \param zbins the number of bins along z
  */
  void set_grid_dimensions(int phibins, int rbins, int zbins);
  void set_grid_dimensions(const int layerBins);

  /// set to true to store evaluation histograms and ntuples
  void set_save_histograms(bool value) { m_savehistograms = value; }

  /// output file name for evaluation histograms
  void set_histogram_outputfile(const std::string& outputfile) { m_histogramfilename = outputfile; }

  /// output file
  /**
   * this is the file where space charge matrix container is stored
   */
  void set_outputfile(const std::string& filename)
  {
    m_outputfile = filename;
  }

  //@}

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// parameters
  void SetDefaultParameters() override;

  /// disable distortion correction
  void disableModuleEdgeCorr() { m_disable_module_edge_corr = true; }
  void disableStaticCorr() { m_disable_static_corr = true; }
  void disableAverageCorr() { m_disable_average_corr = true; }
  void disableFluctuationCorr() { m_disable_fluctuation_corr = true; }

  void set_svtx_track_map_name(std::string name) {_svtx_track_map_name = name;}

 private:
  /// load nodes
  int load_nodes(PHCompositeNode*);

  /// create evaluation histograms
  void create_histograms();

  /// process tracks
  void process_tracks();

  /// returns true if track fulfills basic requirement for distortion calculations
  bool accept_track(SvtxTrack*) const;
  bool checkTrackCM(SvtxTrack *track) const;
  bool checkTPOTResidual(SvtxTrack* track) const;

  /// process track
  void process_track(SvtxTrack*);

  /// get relevant cell for a given cluster
  int get_cell_index(const Acts::Vector3&);
  int get_cell_index_layer(const int layer);
  int get_cell_index_radius(const Acts::Vector3 &loc);
  int get_cell_index_rz(const Acts::Vector3 &loc);

  /// output file
  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  /// true if only tracks with micromegas must be used
  bool m_use_micromegas = true;

  /// minimum pT required for track to be considered in residuals calculation (GeV/c)
  double m_min_pt = 0.5;

  /// track pca z cut, only apply if m_requireCM is enabled
  double m_pcazcut = 10; // +/- 10 cm

  /// track eta cut, only apply if m_requireCM is enabled
  double m_etacut = 0.25; // +/- 0.25

  /// require track crossing zero
  bool m_requireCrossing = false;

  /// require track near CM, only for 1D distortion map (layer or radius)
  bool m_requireCM = false;

  ///@name selection parameters
  //@{
  // residual cuts in r, phi plane
  float m_max_talpha = 0.6;
  float m_max_drphi = 0.5;

  // residual cuts in r, z plane
  float m_max_tbeta = 1.5;
  float m_max_dz = 0.5;
  //@}

  float m_zmin = 0;
  float m_zmax = 0;
  
  // cuts on rphi, z residuals errors
  float m_minRPhiErr = 0.005;  // 0.005cm -- 50um
  float m_minZErr = 0.01;  // 0.01cm -- 100um

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_1D_layer_negz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_1D_layer_posz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_1D_radius_negz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_1D_radius_posz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_2D_radius_z;

  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_clusters = 0;
  int m_accepted_clusters = 0;
  //@}

  ///@name evaluation histograms
  //@{
  /// Output root histograms
  bool m_savehistograms = false;

  /// histogram output file name
  std::string m_histogramfilename = "TpcSpaceChargeReconstruction.root";
  std::unique_ptr<TFile> m_histogramfile;

  using TH1_map_t = std::map<int, TH1*>;
  using TH2_map_t = std::map<int, TH2*>;

  TH1_map_t m_h_drphi;
  TH1_map_t m_h_dz;
  TH2_map_t m_h_drphi_alpha;
  TH2_map_t m_h_dz_beta;
  //@}

  //! tracks
  SvtxTrackMap* m_track_map = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

 ActsGeometry *m_tGeometry = nullptr;

  /// disable distortion correction
  bool m_disable_module_edge_corr = false;
  bool m_disable_static_corr = false;
  bool m_disable_average_corr = false;
  bool m_disable_fluctuation_corr = false;

  //
  std::string _svtx_track_map_name = "SvtxTrackMap";
};

#endif
