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

#include <distortionanalysis/DistortionAnalysis.h>
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
void Fun4All_ZFAllTrackers(
    const int nEvents = 10,
    const std::string clusterfilename = "DST_TRKR_CLUSTER_run2pp_ana517_2024p024_v001-00053877-00000.root",
    const std::string outdir = "root/",
    const std::string outfilename = "clusters_seeds",
    const int index = 0,
    const int stepsize = 10,
    const bool convertSeeds = false)
{
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
  // refer to https://indico.bnl.gov/event/29480/contributions/115897/attachments/65795/113043/Distortions_2025_10_27_hp.pdf
  G4TPC::tpc_drift_velocity_reco = 0.00745455; // cm/ns
  TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
  G4TPC::tpc_tzero_reco = -203.8238;  // ns
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
  G4TRACKING::SC_CALIBMODE = false;
  G4TRACKING::SC_USE_MICROMEGAS = true;
  TRACKING::pp_mode = true;
  
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  
  string outDir = outdir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outfilename + "_" + runnumber + "-" + segment + "-" + index + ".root";
  std::cout<<"outfile "<<outfile<<std::endl;
  std::string theOutfile = outfile.Data();

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4MAGNET::magfield = "0.01";
  G4MAGNET::magfield_tracking = G4MAGNET::magfield;
  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto hitsinclus = new Fun4AllDstInputManager("ClusterInputManager");
  hitsinclus->fileopen(clusterfilename);
  se->registerInputManager(hitsinclus);

  TRACKING::tpc_zero_supp = true;
  /*
  Mvtx_HitUnpacking();
  Intt_HitUnpacking();
  Tpc_HitUnpacking();
  Micromegas_HitUnpacking();

  Mvtx_Clustering();
  Intt_Clustering();

  auto tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(0);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  se->registerSubsystem(tpcclusterizer);

  Micromegas_Clustering();
  */

  // this now does Si and TPC seeding only
  Tracking_Reco_TrackSeed_ZeroField();

  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  // set search windows matching Silicon to TPC seeds
  // Selected for tracks with ntpc>34,|z_Si-z_TPC|<30,crossing==0
  // see https://indico.bnl.gov/event/26202/attachments/59703/102575/2025_01_31_ZeroField.pdf
  // Will probably be narrowed from improved alignment soon.
  silicon_match->window_dx.set_QoverpT_maxabs({2.6,0,0});
  silicon_match->window_dy.set_QoverpT_maxabs({2.3,0,0});
  silicon_match->window_dz.set_QoverpT_range({-2.9,0,0},{4.2,0,0});
  silicon_match->window_deta.set_QoverpT_maxabs({0.06,0,0});
  silicon_match->window_dphi.set_QoverpT_maxabs({0.11,0,0});
  silicon_match->set_test_windows_printout(false);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  silicon_match->zeroField(true);
  se->registerSubsystem(silicon_match);

  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_pp_mode(TRACKING::pp_mode);

  mm_match->set_rphi_search_window_lyr1(3.);
  //mm_match->set_rphi_search_window_lyr1(1.5);//test value
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);            // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(true);  // used for tuning search windows only
  mm_match->zeroField(true);
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
    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);
  }



  PHSimpleVertexFinder *finder = new PHSimpleVertexFinder;
  finder->zeroField(true);
  finder->Verbosity(0);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(-99999.);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000000000);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);

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
    residuals->setMinPt( -9999 );
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
  resid->zeroField();
  resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);
  resid->Verbosity(0);
  se->registerSubsystem(resid);

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

  Enable::QA = false;

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
