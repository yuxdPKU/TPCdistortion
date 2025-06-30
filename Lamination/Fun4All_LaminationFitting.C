/*
 * This macro shows running TPC lamination clustering
 */
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <QA.C>
#include <fun4all/Fun4AllUtils.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <ffamodules/HeadReco.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>

#include <phool/recoConsts.h>

#include <stdio.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libtrackingqa.so)
bool isGood(const string &infile);

void Fun4All_LaminationFitting(
    const int nEvents = 0,
    const int runnumber = 53018,
    const std::string &filelist = "lamination_dst_list_53018",
    const std::string &outfilename = "LAMINATION_Fit",
    const std::string &outdir = "./",
    bool ppmode = true
)
{
  auto se = Fun4AllServer::instance();
  se->Verbosity(0);
  auto rc = recoConsts::instance();

  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  //! flags to set
  Enable::QA = false;
  TRACKING::tpc_zero_supp = true;
  TRACKING::pp_mode = ppmode;

  TString outfile = outdir + outfilename + "_" + runnumber;
  std::cout<<"outfile "<<outfile<<std::endl;
  std::string theOutfile = outfile.Data();

  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  auto laminationinput = new Fun4AllDstInputManager("LaminationInput");
  laminationinput->AddListFile(filelist);
  se->registerInputManager(laminationinput);

  Enable::CDB = true;
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  std::cout << "CDB tracking geometry file "<<geofile << std::endl;
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;

  if(runnumber == 53098)
  {
    G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  }

  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS=false;

  G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
  //G4TPC::average_correction_filename = std::string(Form("/sphenix/tg/tg01/jets/bkimelman/BenProduction/Feb25_2025/Laminations_run2pp_ana466_2024p012_v001-%08d.root",runnumber));
  //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS=false;
  //G4TPC::average_correction_interpolate = false;

  //G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;

  //TRACKING::reco_tpc_maxtime_sample = 1023;

  G4TPC::laser_adc_threshold = 100;

  TpcReadoutInit( runnumber );
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING = true;
  TString laminationoutfile = theOutfile + ".root";
  std::string laminationstring(laminationoutfile.Data());
  G4TPC::LaminationOutputName = laminationstring;

  TString laminationqafile = theOutfile + "_qa.pdf";
  std::string laminationqastring(laminationqafile.Data());
  G4TPC::LaminationQAName = laminationqastring;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();
  
  TPC_LaminationFitting();

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  CDBInterface::instance()->Print();

  if (Enable::QA)
  {
    TString qaname = theOutfile + "_qa.root";
    std::string qaOutputFileName(qaname.Data());
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
