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

#include <trackreco/PHTrackPruner.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>

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
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)
R__LOAD_LIBRARY(libtrack_reco.so)

namespace RecoFitMode
{
  std::string normalize(std::string value)
  {
    std::transform(
        value.begin(), value.end(), value.begin(),
        [](const unsigned char character)
        { return static_cast<char>(std::tolower(character)); });
    value.erase(
        std::remove_if(
            value.begin(), value.end(),
            [](const unsigned char character)
            { return character == '_' || character == '-' || std::isspace(character); }),
        value.end());
    return value;
  }

  bool resolve_use_acts(const std::string& requested_mode, bool& valid)
  {
    valid = true;
    const std::string mode = normalize(requested_mode);
    if (mode == "acts" || mode == "actsfit")
    {
      return true;
    }
    if (mode == "genfit" || mode == "gen")
    {
      return false;
    }
    valid = false;
    return false;
  }
}

void Fun4All_TrackAnalysis(
    const int nEvents = 10,
    const std::string seedfilename = "DST_TRKR_SEED_run3pp_ana573_2026p003_v001-00079516-00000.root",
    const std::string outdir = "root/",
    const std::string outfilename = "polyseed",
    const std::string fitMode = "acts",
    const bool usemms = false,
    const int index = 0,
    const int stepsize = 10)
{
  bool validFitMode = false;
  const bool useActsFit = RecoFitMode::resolve_use_acts(fitMode, validFitMode);
  if (!validFitMode)
  {
    std::cout << "Fun4All_TrackAnalysis - invalid fitMode: " << fitMode
              << " (expected acts, actsfit, genfit, or gen)" << std::endl;
    return;
  }

  std::cout << "Fun4All_TrackAnalysis - fit mode: "
            << (useActsFit ? "actsfit" : "genfit") << std::endl;

  std::pair<int, int>
      runseg = Fun4AllUtils::GetRunSegment(seedfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;

  auto *rc = recoConsts::instance();
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
  G4TRACKING::SC_USE_MICROMEGAS = usemms;
  TRACKING::streaming_mode = true;
  
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

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = false;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
   // to use a custom file instead of the database file:
  G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  std::cout<<"Average distortion map used: "<<G4TPC::average_correction_filename<<std::endl;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();

  auto *hitsinseed = new Fun4AllDstInputManager("SeedInputManager");
  hitsinseed->fileopen(seedfilename);
  se->registerInputManager(hitsinseed);

  std::string tpcresidstring;

  Reject_Laser_Events();

  const std::string& clusterMapName = "TRKR_CLUSTER_SEED";

  Tracking_Reco_TrackMatching_run2pp(clusterMapName);
  //Tracking_Reco_TrackFit_run2pp("", clusterMapName);
  //Tracking_Reco_Vertex_run2pp(clusterMapName);

  auto deltazcorr = new PHTpcDeltaZCorrection;
  deltazcorr->Verbosity(0);
  deltazcorr->setTrkrClusterContainerName(clusterMapName);
  se->registerSubsystem(deltazcorr);

  if (useActsFit)
  {
    std::cout << "Using ACTS fit" << std::endl;

    // The first pass must use the full detector so that the track pruner
    // can apply its TPC and Micromegas cluster/state requirements.
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    actsFit->setTrkrClusterContainerName(clusterMapName);
    actsFit->fitSiliconMMs(false);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::streaming_mode);
    actsFit->set_use_clustermover(true);
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::streaming_mode);
    se->registerSubsystem(cleaner);

    auto trackpruner = new PHTrackPruner;
    trackpruner->Verbosity(0);
    trackpruner->set_cluster_map_name(clusterMapName);
    trackpruner->set_pruned_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
    trackpruner->set_track_pt_low_cut(0.5);
    trackpruner->set_track_quality_high_cut(100);
    trackpruner->set_nmvtx_clus_low_cut(3);
    trackpruner->set_nintt_clus_low_cut(2);
    trackpruner->set_ntpc_clus_low_cut(35);
    if (G4TRACKING::SC_USE_MICROMEGAS) { trackpruner->set_ntpot_clus_low_cut(1); }
    else { trackpruner->set_ntpot_clus_low_cut(0); }
    trackpruner->set_nmvtx_states_low_cut(3);
    trackpruner->set_nintt_states_low_cut(2);
    trackpruner->set_ntpc_states_low_cut(35);
    if (G4TRACKING::SC_USE_MICROMEGAS) { trackpruner->set_ntpot_states_low_cut(1); }
    else { trackpruner->set_ntpot_states_low_cut(0); }
    se->registerSubsystem(trackpruner);

    auto actsFit_SiTpotFit = new PHActsTrkFitter;
    actsFit_SiTpotFit->Verbosity(0);
    actsFit_SiTpotFit->commissioning(G4TRACKING::use_alignment);
    actsFit_SiTpotFit->setTrkrClusterContainerName(clusterMapName);
    actsFit_SiTpotFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit_SiTpotFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit_SiTpotFit->set_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
    actsFit_SiTpotFit->set_pp_mode(TRACKING::streaming_mode);
    actsFit_SiTpotFit->set_use_clustermover(true);
    actsFit_SiTpotFit->useActsEvaluator(false);
    actsFit_SiTpotFit->useOutlierFinder(false);
    actsFit_SiTpotFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit_SiTpotFit);
  }
  else
  {
    std::cout << "Using GENFIT" << std::endl;

    // Full-detector fit used to populate SvtxTrackMap for pruning.
    auto genfitFit = new PHGenFitTrkFitter;
    // need setter for TRKR_CLUSTER_SEED
    genfitFit->set_fit_silicon_mms(false);
    genfitFit->set_use_micromegas(G4TRACKING::SC_USE_MICROMEGAS);
    se->registerSubsystem(genfitFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::streaming_mode);
    se->registerSubsystem(cleaner);

    auto trackpruner = new PHTrackPruner;
    trackpruner->Verbosity(0);
    trackpruner->set_cluster_map_name(clusterMapName);
    trackpruner->set_pruned_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
    trackpruner->set_track_pt_low_cut(0.5);
    trackpruner->set_track_quality_high_cut(100);
    trackpruner->set_nmvtx_clus_low_cut(3);
    trackpruner->set_nintt_clus_low_cut(2);
    trackpruner->set_ntpc_clus_low_cut(35);
    if (G4TRACKING::SC_USE_MICROMEGAS) { trackpruner->set_ntpot_clus_low_cut(1); }
    else { trackpruner->set_ntpot_clus_low_cut(0); }
    trackpruner->set_nmvtx_states_low_cut(3);
    trackpruner->set_nintt_states_low_cut(2);
    trackpruner->set_ntpc_states_low_cut(35);
    if (G4TRACKING::SC_USE_MICROMEGAS) { trackpruner->set_ntpot_states_low_cut(1); }
    else { trackpruner->set_ntpot_states_low_cut(0); }
    se->registerSubsystem(trackpruner);

    auto genfitFit_SiTpotFit = new PHGenFitTrkFitter;
    // need setter for TRKR_CLUSTER_SEED
    genfitFit_SiTpotFit->set_fit_silicon_mms(G4TRACKING::SC_CALIBMODE);
    genfitFit_SiTpotFit->set_use_micromegas(G4TRACKING::SC_USE_MICROMEGAS);
    genfitFit_SiTpotFit->set_svtx_track_map_name("SvtxSiliconMMTrackMap");
    genfitFit_SiTpotFit->set_svtx_seed_map_name("PrunedSvtxTrackSeedContainer");
    se->registerSubsystem(genfitFit_SiTpotFit);
  }

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
    residuals->setClusterMapName("TRKR_CLUSTER_SEED");
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
    residuals->Verbosity(0);
    se->registerSubsystem(residuals);
  }

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

  Enable::QA = true;

  if (Enable::QA)
  {
    Distortions_QA(G4TRACKING::SC_USE_MICROMEGAS);
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
