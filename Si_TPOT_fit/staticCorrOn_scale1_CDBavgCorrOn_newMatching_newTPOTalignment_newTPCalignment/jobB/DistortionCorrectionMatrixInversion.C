#include <tpccalib/TpcSpaceChargeMatrixInversion.h>
#include <tpccalib/TpcSpaceChargeMatrixInversion1D.h>
#include <tpccalib/TpcSpaceChargeMatrixInversion2D.h>
#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <cstdio>
#include <sstream>

R__LOAD_LIBRARY(libtpccalib.so)

//_______________________________________________
// get list of files matching selection
std::vector<std::string> list_files( const std::string& selection )
{
  std::vector<std::string> out;

  std::cout << "list_files - selection: " << selection << std::endl;
  if( selection.empty() ) return out;

  const std::string command = std::string("ls -1 ") + selection;
  auto tmp = popen( command.c_str(), "r" );
  char line[512];
  while( fgets( line, 512, tmp ) )
  {

    std::istringstream istr( line );
    std::string filename;
    istr >> filename;

    if( filename.empty() ) continue;
    if( access( filename.c_str(), R_OK ) ) continue;

    out.push_back( filename );
  }
  pclose( tmp );
  return out;
}

//_______________________________________________
void DistortionCorrectionMatrixInversion(int run=53285)
{
  TpcSpaceChargeReconstructionHelper::set_phi_range_central( {-1.73246,-1.43608} );
  TpcSpaceChargeReconstructionHelper::set_phi_range_east( {-2.26272,-1.96089} );
  TpcSpaceChargeReconstructionHelper::set_phi_range_west( {-1.21241,-0.909953} );

  TpcSpaceChargeReconstructionHelper::set_theta_range_central( {{-0.92,-0.613404}, {-0.568656,-0.0325273}, {0.0276558,0.566484}, {0.609928,0.917104}} );
  //TpcSpaceChargeReconstructionHelper::set_theta_range_central( {{-0.788212,-0.613404}, {-0.568656,-0.0325273}, {0.0276558,0.566484}, {0.609928,0.917104}} );
  TpcSpaceChargeReconstructionHelper::set_theta_range_east( {{-0.641231,-0.140061}, {0.133491,0.638103}} );
  TpcSpaceChargeReconstructionHelper::set_theta_range_west( {{-0.643933,-0.14138}, {0.13084,0.63809}} );

  TpcSpaceChargeReconstructionHelper::Verbosity( 2 );

  // input files
  /*
   * this is the list of distortion correction matrices files coming out from Job A
   * that needs to be inverted, to get track-based, beam-induced distortions inside
   * TPOT acceptance
   */
  //const TString tag = "_flat_genfit_truth_notpc_distorted-new";

//53877 - 400khz
//53876 - 430khz
//53756 - 380khz
//53744 - 300khz
//53630 - 550khz
//53534 - 250khz
//53285 - 70khz

  //int run = 53285;
  //int run = 53630;
  //int run = 53756;
  //int run = 53877;
  //int run = 53534;
  //int run = 53744;
  //int run = 53876;

  //const TString inputFile = Form( "DST/CONDOR%s/TpcSpaceChargeMatrices%s_*.root", tag.Data(), tag.Data() );
  //const TString inputFile = Form( "/sphenix/u/xyu3/hftg01/DST_FOR_DISTORTION/Reconstructed/%d/clusters_seeds_%d-*.root_PhTpcResiduals.root", run, run );
  const TString inputFile = Form( "../Reconstructed/%d/clusters_seeds_%d-*.root_PhTpcResiduals.root", run, run );

  // Central membrane distortion corrections
  /*
   * this is the 2D distortion corrections measured at the central membrane using diffuse lasers
   * it is used to extrapolate the distortions measured in the TPOT acceptance to the rest of the TPC acceptance
   * see: https://indico.bnl.gov/event/22887/contributions/90413/attachments/54020/92443/distortion_extrapolation_hp.pdf
   */
  //const std::string inputfile_cm = "distortion_maps/average_minus_static_distortion_cm.root";
  const std::string inputfile_cm = "/sphenix/tg/tg01/jets/bkimelman/BenProduction/Jan23_2025/adc_100/Laminations_52622_beamOn_-1ZSChris_2.00MIP_00011.root";

  // output file
  const TString outputFile = Form( "Rootfiles/Distortions_full_mm_%d.root", run );
  const TString outputFile_phi = Form( "Rootfiles/Distortions_full_mm_%d_phi.root", run );
  const TString outputFile_z = Form( "Rootfiles/Distortions_full_mm_%d_z.root", run );
  const TString outputFile1D = Form( "Rootfiles/Distortions_1D_mm_%d", run );
  const TString outputFile2D = Form( "Rootfiles/Distortions_2D_mm_%d_rz.root", run );

  std::cout << "DistortionCorrectionMatrixInversion - inputFile: " << inputFile << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - inputfile_cm: " << inputfile_cm << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - outputFile: " << outputFile << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - outputFile_phi: " << outputFile_phi << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - outputFile_z: " << outputFile_z << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - outputFile1D: " << outputFile1D << std::endl;
  std::cout << "DistortionCorrectionMatrixInversion - outputFile2D: " << outputFile2D << std::endl;

  auto filenames = list_files( inputFile.Data() );
  std::cout << "SpaceChargeMatrixInversion - loaded " << filenames.size() << " files" << std::endl;

  /*
  std::cout << "DistortionCorrectionMatrixInversion - 1D begin." << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion1D spaceChargeMatrixInversion1D;
  spaceChargeMatrixInversion1D.Verbosity(1);

  // load input files
  for( const auto& file:filenames )
  { spaceChargeMatrixInversion1D.add_from_file( file , "TpcSpaceChargeMatrixContainer_1D"); }

  // calculate the distortions in TPOT acceptance
  spaceChargeMatrixInversion1D.calculate_distortion_corrections();

  // write to output
  spaceChargeMatrixInversion1D.save_distortion_corrections( outputFile1D.Data() );

  std::cout << "DistortionCorrectionMatrixInversion - 1D done." << std::endl;
  */

  std::cout << "DistortionCorrectionMatrixInversion - 2D begin." << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion2D spaceChargeMatrixInversion2D;
  spaceChargeMatrixInversion2D.Verbosity(0);
  spaceChargeMatrixInversion2D.set_min_cluster_count(50);

  // load input files
  for( const auto& file:filenames )
  { spaceChargeMatrixInversion2D.add_from_file( file , "TpcSpaceChargeMatrixContainer_2D_radius_z"); }

  // calculate the distortions in TPOT acceptance
  spaceChargeMatrixInversion2D.calculate_distortion_corrections();

  // load central membrane corrections
  //spaceChargeMatrixInversion.load_cm_distortion_corrections( inputfile_cm );
  spaceChargeMatrixInversion2D.extrapolate_distortion_corrections();

  // write to output
  spaceChargeMatrixInversion2D.save_distortion_corrections( outputFile2D.Data() );

  std::cout << "DistortionCorrectionMatrixInversion - 3D begin." << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion spaceChargeMatrixInversion;
  spaceChargeMatrixInversion.Verbosity(0);
  spaceChargeMatrixInversion.set_min_cluster_count(15);

  // load input files
  for( const auto& file:filenames )
  { spaceChargeMatrixInversion.add_from_file( file , "TpcSpaceChargeMatrixContainer"); }

  // calculate the distortions in TPOT acceptance
  spaceChargeMatrixInversion.calculate_distortion_corrections();

  // load central membrane corrections
  //spaceChargeMatrixInversion.load_cm_distortion_corrections( inputfile_cm );
  //spaceChargeMatrixInversion.extrapolate_distortion_corrections();

  // write to output
  spaceChargeMatrixInversion.save_distortion_corrections( outputFile.Data() );

  std::cout << "DistortionCorrectionMatrixInversion - all done." << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion spaceChargeMatrixInversion_phi;
  spaceChargeMatrixInversion_phi.Verbosity(0);
  spaceChargeMatrixInversion_phi.set_min_cluster_count(15);

  // load input files
  for( const auto& file:filenames )
  { spaceChargeMatrixInversion_phi.add_from_file( file , "TpcSpaceChargeMatrixContainer"); }

  // calculate the distortions in TPOT acceptance
  spaceChargeMatrixInversion_phi.calculate_distortion_corrections();

  // load central membrane corrections
  //spaceChargeMatrixInversion_phi.load_cm_distortion_corrections( inputfile_cm );
  //spaceChargeMatrixInversion_phi.extrapolate_distortion_corrections();

  // write to output
  spaceChargeMatrixInversion_phi.save_distortion_corrections( outputFile_phi.Data() );

  std::cout << "DistortionCorrectionMatrixInversion -- phi - all done." << std::endl;

  // perform matrix inversion
  TpcSpaceChargeMatrixInversion spaceChargeMatrixInversion_r;
  spaceChargeMatrixInversion_r.Verbosity(0);
  spaceChargeMatrixInversion_r.set_min_cluster_count(15);

  // load input files
  for( const auto& file:filenames )
  { spaceChargeMatrixInversion_r.add_from_file( file , "TpcSpaceChargeMatrixContainer"); }

  // calculate the distortions in TPOT acceptance
  spaceChargeMatrixInversion_r.calculate_distortion_corrections();

  // load central membrane corrections
  //spaceChargeMatrixInversion_r.load_cm_distortion_corrections( inputfile_cm );
  //spaceChargeMatrixInversion_r.extrapolate_distortion_corrections();

  // write to output
  spaceChargeMatrixInversion_r.save_distortion_corrections( outputFile_r.Data() );

  std::cout << "DistortionCorrectionMatrixInversion -- r - all done." << std::endl;


}
