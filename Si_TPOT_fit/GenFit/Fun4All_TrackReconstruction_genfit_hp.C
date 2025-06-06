#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>
#include <qautils/QAHistManagerDef.h>

#include <simqa_modules/QAG4SimulationDistortions.h>

#include "Trkr_RecoInit.C"
#include "Trkr_TpcReadoutInit.C"

#include "Trkr_Reco.C"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libsimqa_modules.so)

//____________________________________________________________________
int Fun4All_TrackReconstruction_genfit_hp(
  const int nEvents = 1000,
  const int nSkipEvents = 0,
  const char* tag = "ana475_2024p018_v001",
  const char* tag_clusters = "ana466_2024p012_v001",
  const int runnumber = 53877,
  const int segment = 2,
  const char* outputFile =  "TRACKS-00053877-0002.root",
  const char* residualsFile = "TpcResiduals-00053285-0002.root",
  const char* qaFile = "DistortionQA-00053285-0002.root"
  )
{
  // print inputs
  std::cout << "Fun4All_TrackReconstruction - nEvents: " << nEvents << std::endl;
  std::cout << "Fun4All_TrackReconstruction - nSkipEvents: " << nSkipEvents << std::endl;
  std::cout << "Fun4All_TrackReconstruction - tag: " << tag << std::endl;
  std::cout << "Fun4All_TrackReconstruction - tag_clusters: " << tag_clusters << std::endl;
  std::cout << "Fun4All_TrackReconstruction - runnumber: " << runnumber << std::endl;
  std::cout << "Fun4All_TrackReconstruction - segment: " << segment << std::endl;
  std::cout << "Fun4All_TrackReconstruction - outputFile: " << outputFile << std::endl;
  std::cout << "Fun4All_TrackReconstruction - residualsFile: " << residualsFile << std::endl;

  // get run range
  const int range_begin = int( runnumber/100 ) * 100;
  const int range_end = int( runnumber/100 + 1 ) * 100;

  // generate track input file
  const auto inputFile = Form( "/sphenix/lustre01/sphnxpro/production/run2pp/physics/%s/DST_TRKR_SEED/run_%08i_%08i/dst/DST_TRKR_SEED_run2pp_%s-%08i-%05i.root",
    tag, range_begin, range_end,
    tag, runnumber, segment );
  std::cout << "Fun4All_TrackReconstruction - inputFile: " << inputFile << std::endl;

  // generate cluster input file
  const auto inputFile_clusters = Form( "/sphenix/lustre01/sphnxpro/production/run2pp/physics/%s/DST_TRKR_CLUSTER/run_%08i_%08i/dst/DST_TRKR_CLUSTER_run2pp_%s-%08i-%05i.root",
    tag_clusters, range_begin, range_end,
    tag_clusters, runnumber, segment );
  std::cout << "Fun4All_TrackReconstruction - inputFile_clusters: " << inputFile_clusters << std::endl;

  TRACKING::pp_mode = true;
  G4TRACKING::SC_CALIBMODE = true;
  G4TRACKING::convert_seeds_to_svtxtracks = false;

  // condition database
  Enable::CDB = true;

  // reco const
  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  // tpc readout initialization
  TpcReadoutInit( runnumber );

  // printout
  std::cout<< "Fun4All_TrackReconstruction - samples: " << TRACKING::reco_tpc_maxtime_sample << std::endl;
  std::cout<< "Fun4All_TrackReconstruction - pre: " << TRACKING::reco_tpc_time_presample << std::endl;
  std::cout<< "Fun4All_TrackReconstruction - vdrift: " << G4TPC::tpc_drift_velocity_reco << std::endl;

  // server
  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  PHRandomSeed::Verbosity(1);

  // tracking geometry
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  std::cout << "Fun4All_TrackReconstruction_hp - geofile: " << geofile << std::endl;

  auto ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  // distortion correction
  if(true)
  {
    // module edge corrections
    G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

    // static distortions
    G4TPC::ENABLE_STATIC_CORRECTIONS = true;

    // average distortions
    G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
  }

  // tpc zero suppression
  TRACKING::tpc_zero_supp = true;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  // input managers
  {
    auto in = new Fun4AllDstInputManager("DST_clusters_in");
    in->fileopen(inputFile_clusters);
    se->registerInputManager(in);
  }

  {
    auto in = new Fun4AllDstInputManager("DST_seed_in");
    in->fileopen(inputFile);
    se->registerInputManager(in);
  }

  std::cout << "Fun4All_TrackReconstruction_new_hp - done with input managers" << std::endl;

  {
    // matching to silicons
    auto silicon_match = new PHSiliconTpcTrackMatching;
    silicon_match->Verbosity(0);

    // narrow matching windows
    silicon_match->set_x_search_window(0.36);
    silicon_match->set_y_search_window(0.36);
    silicon_match->set_z_search_window(2.5);
    silicon_match->set_phi_search_window(0.014);
    silicon_match->set_eta_search_window(0.0091);
    silicon_match->set_test_windows_printout(false);

    silicon_match->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(silicon_match);
  }

  if( true )
  {
    // matching with micromegas
    auto mm_match = new PHMicromegasTpcTrackMatching;
    mm_match->Verbosity(0);
    mm_match->set_pp_mode(TRACKING::pp_mode);

    mm_match->set_rphi_search_window_lyr1(3.0);
    mm_match->set_rphi_search_window_lyr2(15.0);

    mm_match->set_z_search_window_lyr1(30.0);
    mm_match->set_z_search_window_lyr2(3.0);

    mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
    mm_match->set_test_windows_printout(false);  // used for tuning search windows only
    se->registerSubsystem(mm_match);
  }

  {
    // track fit
    se->registerSubsystem(new PHTpcDeltaZCorrection);

    // perform final track fit with GENFIT
    auto genfitFit = new PHGenFitTrkFitter;
    genfitFit->set_fit_silicon_mms(G4TRACKING::SC_CALIBMODE);
    genfitFit->set_use_micromegas(G4TRACKING::SC_USE_MICROMEGAS);
    se->registerSubsystem(genfitFit);
  }

  if (G4TRACKING::SC_CALIBMODE)
  {
    // Genfit based Tpc space charge Reconstruction
    auto tpcSpaceChargeReconstruction = new TpcSpaceChargeReconstruction;
    tpcSpaceChargeReconstruction->set_use_micromegas(G4TRACKING::SC_USE_MICROMEGAS);
    tpcSpaceChargeReconstruction->set_outputfile(residualsFile);

    // reconstructed distortion grid size (phi, r, z)
    tpcSpaceChargeReconstruction->set_grid_dimensions(36, 48, 80);
    se->registerSubsystem(tpcSpaceChargeReconstruction);
  }

  if (G4TRACKING::SC_CALIBMODE)
  {
    auto simulationDistortionQA = new QAG4SimulationDistortions();
    simulationDistortionQA->set_trackmap_name("SvtxTrackMap");
    se->registerSubsystem(simulationDistortionQA);
  }

  // output manager
  if( false )
  {
    auto out = new Fun4AllDstOutputManager("DSTOUT", outputFile);
    se->registerOutputManager(out);
  }

  // skip events if any specified
  if( nSkipEvents > 0 )
  { se->skip( nSkipEvents ); }

  // run
  se->run(nEvents);

  // terminate
  se->End();
  se->PrintTimer();

  // QA
  QAHistManagerDef::saveQARootFile(qaFile);

  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}

// This function is only used to test if we can load this as root6 macro
// without running into unresolved libraries and include files
void RunFFALoadTest() {}
