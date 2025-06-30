#ifndef TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
#define TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
/**
 * \file TpcSpaceChargeMatrixInversion.h
 * \brief aggregate space charge reconstruction matrices from several jobs and perform the inversion
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcSpaceChargeMatrixContainer.h"

#include <fun4all/Fun4AllBase.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <memory>

class TH3;

/**
 * \class TpcSpaceChargeMatrixInversion
 * \brief performs space charge distortion reconstruction using tracks
 */

class TpcSpaceChargeMatrixInversion : public Fun4AllBase
{
 public:
  /// constructor
  TpcSpaceChargeMatrixInversion(const std::string& = "TPCSPACECHARGEMATRIXINVERSION");

  ///@name modifiers
  //@{

  /// load central membrane distortion correction
  void load_cm_distortion_corrections(const std::string& /*filename*/);

  /// load average distortion correction
  void load_average_distortion_corrections(const std::string& /*filename*/);

  /// add space charge correction matrix to current. Returns true on success
  bool add(const TpcSpaceChargeMatrixContainer&);

  /// add space charge correction matrix, loaded from file, to current. Returns true on success
  bool add_from_file(const std::string& /*filename*/, const std::string& /*objectname*/ = "TpcSpaceChargeMatrixContainer");

  enum class InversionMode
  {
    FullInversion,        // use 3D matrices (phi,z,r)
    ReducedInversion_phi, // use 3D reduced matrices, use phi residuals to get delta_r
    ReducedInversion_z    // use 3D reduced matrices, use z residuals to get delta_r
  };

  /// calculate distortions by inverting stored matrices, and save relevant histograms
  void calculate_distortion_corrections(const InversionMode = InversionMode::FullInversion);

  /// extrapolate distortions
  void extrapolate_distortion_corrections();

  /// save distortions
  void save_distortion_corrections(const std::string& /*filename*/ = "DistortionCorrections.root");

  // set min cluster count cut
  void set_min_cluster_count(int min) {min_cluster_count = min;}

  //@}

 private:
  TH3* transform_dphi_rdphi(const TH3* source, const TString& name);

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

  /// output distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_average;

  /// central membrane distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_cm;

  std::array<TH3*, 2> hDRPint = {{nullptr, nullptr}};
  int min_cluster_count = 5;
};

#endif
