/*
 * This macro shows a working example of running TrackSeeding over the cluster DST
 * This has track residuals as default output but has KFParticle set up with a togglable flag
 * with the default set up for K Short reconstruction
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/TpcSeedsQA.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>

//#include <distortionanalysis/DistortionAnalysis.h>
#include <trackreco/PHTrackPruner.h>

#include <stdio.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libDistortionAnalysis.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
void Fun4All_TrackAnalysis(
    const int nEvents = 10,
    const std::string clusterfilename = "DST_TRKR_CLUSTER_run2pp_ana494_2024p021_v001-00053877-00000.root",
    const std::string clusterdir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana494_2024p021_v001/DST_TRKR_CLUSTER/run_00053800_00053900/dst/",
    const std::string seedfilename = "DST_TRKR_SEED_run2pp_ana494_2024p021_v001-00053877-00000.root",
    const std::string seeddir = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana494_2024p021_v001/DST_TRKR_SEED/run_00053800_00053900/dst/",
    const std::string outdir = "root/",
    const std::string outfilename = "clusters_seeds",
    const int index = 0,
    const int stepsize = 10,
    const bool convertSeeds = false)
{
  std::string inputclusterFile = clusterdir + clusterfilename;
  std::string inputseedFile = seeddir + seedfilename;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
  std::pair<int, int>
      runseg = Fun4AllUtils::GetRunSegment(clusterfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  TpcReadoutInit(runnumber);
  // these lines show how to override the drift velocity and time offset values set in TpcReadoutInit
  // G4TPC::tpc_drift_velocity_reco = 0.0073844; // cm/ns
  // TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
  // G4TPC::tpc_tzero_reco = -5*50;  // ns
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  G4TRACKING::SC_CALIBMODE = true;
  G4TRACKING::SC_USE_MICROMEGAS = false;
  TRACKING::pp_mode = true;
  
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  
  string outDir = outdir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outfilename + "_" + runnumber + "-" + segment + ".root";
  std::cout<<"outfile "<<outfile<<std::endl;
  std::string theOutfile = outfile.Data();

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
   // to use a custom file instead of the database file:
  G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  std::cout<<"Average distortion map used: "<<G4TPC::average_correction_filename<<std::endl;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsinseed = new Fun4AllDstInputManager("SeedInputManager");
  hitsinseed->fileopen(inputseedFile);
  se->registerInputManager(hitsinseed);

  auto hitsinclus = new Fun4AllDstInputManager("ClusterInputManager");
  hitsinclus->fileopen(inputclusterFile);
  se->registerInputManager(hitsinclus);

  // Always apply preliminary distortion corrections to TPC clusters before silicon matching
  // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
  auto prelim_distcorr = new PrelimDistortionCorrection;
  prelim_distcorr->set_pp_mode(true);
  prelim_distcorr->Verbosity(0);
  se->registerSubsystem(prelim_distcorr);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  if(G4TPC::ENABLE_AVERAGE_CORRECTIONS)
  {
    // for general tracking
    // Eta/Phi window is determined by 3 sigma window
    // X/Y/Z window is determined by 4 sigma window
    silicon_match->window_deta.set_posQoverpT_maxabs({-0.014,0.0331,0.48});
    silicon_match->window_deta.set_negQoverpT_maxabs({-0.006,0.0235,0.52});
    silicon_match->set_deltaeta_min(0.03);
    silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.15,0,0});
    silicon_match->window_dx.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dy.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dz.set_posQoverpT_maxabs({1.138,0.3919,0.84});
    silicon_match->window_dz.set_negQoverpT_maxabs({0.719,0.6485,0.65});
    silicon_match->set_crossing_deltaz_max(30);
    silicon_match->set_crossing_deltaz_min(2);

    // for distortion correction using SI-TPOT fit and track pT>0.5
    if (G4TRACKING::SC_CALIBMODE)
    {
      silicon_match->window_deta.set_posQoverpT_maxabs({0.016,0.0060,1.13});
      silicon_match->window_deta.set_negQoverpT_maxabs({0.022,0.0022,1.44});
      silicon_match->set_deltaeta_min(0.03);
      silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.09,0,0});
      silicon_match->window_dx.set_QoverpT_maxabs({2.0,0,0});
      silicon_match->window_dy.set_QoverpT_maxabs({1.5,0,0});
      silicon_match->window_dz.set_posQoverpT_maxabs({1.213,0.0211,2.09});
      silicon_match->window_dz.set_negQoverpT_maxabs({1.307,0.0001,4.52});
      silicon_match->set_crossing_deltaz_min(1.2);
    }
  }
  silicon_match->print_windows(true);
  se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_pp_mode(TRACKING::pp_mode);
  //mm_match->set_rphi_search_window_lyr1(3.);
  mm_match->set_rphi_search_window_lyr1(1.5);//test value
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);


  std::string tpcresidstring;
  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("SvtxTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
//    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
//    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->fitSiliconMMs(false);//by default
    actsFit->setUseMicromegas(true);//by default
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

  }

  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  
  //new cuts
  finder->setDcaCut(0.05);
  finder->setTrackPtCut(0.1);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(300);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.10);
  
  se->registerSubsystem(finder);

  // Propagate track positions to the vertex position
  auto vtxProp = new PHActsVertexPropagator;
  vtxProp->Verbosity(0);
  vtxProp->fieldMap(G4MAGNET::magfield_tracking);
  se->registerSubsystem(vtxProp);
 
  //prune acts full tracks, create new SvtxTrackMap
  auto trackpruner = new PHTrackPruner;
  trackpruner->Verbosity(0);
  trackpruner->set_pruned_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
  trackpruner->set_track_pt_low_cut(0.5);
  trackpruner->set_track_quality_high_cut(100);
  trackpruner->set_nmvtx_clus_low_cut(3);
  trackpruner->set_nintt_clus_low_cut(2);
  trackpruner->set_ntpc_clus_low_cut(35);
  trackpruner->set_ntpot_clus_low_cut(0);
  trackpruner->set_nmvtx_states_low_cut(3);
  trackpruner->set_nintt_states_low_cut(2);
  trackpruner->set_ntpc_states_low_cut(35);
  trackpruner->set_ntpot_states_low_cut(0);
  se->registerSubsystem(trackpruner);

  // perform final track fit with ACTS
  // Si-TPOT fit
  auto actsFit_SiTpotFit = new PHActsTrkFitter;
  actsFit_SiTpotFit->Verbosity(0);
  actsFit_SiTpotFit->commissioning(G4TRACKING::use_alignment);
  // in calibration mode, fit only Silicons and Micromegas hits
  actsFit_SiTpotFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
  actsFit_SiTpotFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
  actsFit_SiTpotFit->set_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
  actsFit_SiTpotFit->set_pp_mode(TRACKING::pp_mode);
  actsFit_SiTpotFit->set_use_clustermover(true);  // default is true for now
  actsFit_SiTpotFit->useActsEvaluator(false);
  actsFit_SiTpotFit->useOutlierFinder(false);
  actsFit_SiTpotFit->setFieldMap(G4MAGNET::magfield_tracking);
  se->registerSubsystem(actsFit_SiTpotFit);

  if (G4TRACKING::SC_CALIBMODE)
  {
    /*
    * in calibration mode, calculate residuals between TPC and fitted tracks,
    * store in dedicated structure for distortion correction
    */
    auto residuals = new PHTpcResiduals;
    const TString tpc_residoutfile = theOutfile + "_PhTpcResiduals.root";
    tpcresidstring = tpc_residoutfile.Data();
    residuals->setOutputfile(tpc_residoutfile.Data());
    residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    residuals->disableAverageCorr();

    // matches Tony's analysis
    residuals->setMinPt( 0.5 );
    residuals->requireCrossing(false);
    residuals->requireCM(true);
    residuals->setPCAzcut(10);
    residuals->setEtacut(0.25);

    residuals->setMaxTrackAlpha(0.6);
    residuals->setMaxTrackBeta(1.5);
    residuals->setMaxTrackResidualDrphi(2);
    residuals->setMaxTrackResidualDz(5);

    residuals->setMinRPhiErr(0.005);
    residuals->setMinZErr(0.01);

    // reconstructed distortion grid size (layer)
    residuals->setGridDimensions(48);

    // reconstructed distortion grid size (phi, r, z)
    residuals->setGridDimensions(36, 16, 80);
    se->registerSubsystem(residuals);
  }

  TString residoutfile = theOutfile + "_resid.root";
  std::string residstring(residoutfile.Data());

  /*
  //auto resid = new TrackResiduals("TrackResiduals");
  auto resid = new DistortionAnalysis("DistortionAnalysis");
  resid->outfileName(residstring);
  resid->alignment(false);

  // adjust track map name
  if(G4TRACKING::SC_CALIBMODE && !G4TRACKING::convert_seeds_to_svtxtracks)
  {
    resid->trackmapName("SvtxSiliconMMTrackMap");
    if( G4TRACKING::SC_USE_MICROMEGAS )
    { resid->set_doMicromegasOnly(true); }
  }

  //resid->clusterTree();
  //resid->hitTree();
  resid->noEventTree();
  resid->noVertexTree();
  resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);

  resid->Verbosity(0);
  se->registerSubsystem(resid);
  */

  TString dstfile = theOutfile + "_dst.root";
  std::string dststring(dstfile.Data());
  /*
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("out", dststring);
  out->AddNode("Sync");
  out->AddNode("EventHeader");
  out->AddNode("PrunedSvtxTrackSeedContainer");
  out->AddNode("SvtxSiliconMMTrackMap");
  se->registerOutputManager(out);
  */

  //auto ntuplizer = new TrkrNtuplizer("TrkrNtuplizer");
  //se->registerSubsystem(ntuplizer);

  Enable::QA = true;

  // Fun4AllOutputManager *out = new Fun4AllDstOutputManager("out", "/sphenix/tg/tg01/hf/jdosbo/tracking_development/Run24/Beam/41626/hitsets.root");
  // se->registerOutputManager(out);
  if (Enable::QA)
  {
    Distortions_QA();
  }
  se->skip(stepsize*index);
  se->run(nEvents);
  se->End();
  se->PrintTimer();
  CDBInterface::instance()->Print();

  std::string qaOutputFileName;
  if (Enable::QA)
  {
    TString qaname = theOutfile + "_qa.root";
    qaOutputFileName = qaname.Data();
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }

  ifstream file_trackresid(residstring.c_str(), ios::binary | ios::ate);
  if (file_trackresid.good() && (file_trackresid.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + residstring + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_tpcresid(tpcresidstring.c_str(), ios::binary | ios::ate);
  if (file_tpcresid.good() && (file_tpcresid.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + tpcresidstring + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_qa(qaOutputFileName.c_str(), ios::binary | ios::ate);
  if (file_qa.good() && (file_qa.tellg() > 100))
  {
    string outputDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + qaOutputFileName + " " + outputDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_dst(dststring.c_str(), ios::binary | ios::ate);
  if (file_dst.good() && (file_dst.tellg() > 100))
  {
    string outputDstDirMove = outdir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputDstDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + dststring + " " + outputDstDirMove;
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
