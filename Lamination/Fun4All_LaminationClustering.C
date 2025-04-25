/*
 * This macro shows running TPC lamination clustering
 */
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_RecoInit.C>
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

void Fun4All_LaminationClustering(
    const int nEvents = 0,
    const int runnumber = 54966,
    const int segment = 0,
    const std::string &filelist = "filelist",
    const std::string &outfilename = "DST_LAMINATION_CLUSTER",
    const std::string &outdir = "./",
    const int index = 0,
    const int stepsize = 10
)
{
  auto se = Fun4AllServer::instance();
  se->Verbosity(1);
  auto rc = recoConsts::instance();

  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");

  //! flags to set
  Enable::QA = false;
  TRACKING::tpc_zero_supp = true;
  TRACKING::pp_mode = false;

  string outDir = outdir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outDir;
  system(makeDirectory.c_str());
  TString outfile = outDir + outfilename + "_" + runnumber + "-" + segment + "-" + index;
  std::cout<<"outfile "<<outfile<<std::endl;
  std::string theOutfile = outfile.Data();

  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  std::ifstream ifs(filelist);
  std::string filepath;
  int i = 0;
  while(std::getline(ifs,filepath))
    {
      std::cout << "Adding DST with filepath: " << filepath << std::endl; 
      if(i==0)
	{
	   std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(filepath);
	   int runNumber = runseg.first;
	   int segment = runseg.second;
	   rc->set_IntFlag("RUNNUMBER", runNumber);
	   rc->set_uint64Flag("TIMESTAMP", runNumber);
	}
      std::string inputname = "InputManager" + std::to_string(i);
      auto hitsin = new Fun4AllDstInputManager(inputname);
      hitsin->fileopen(filepath);
      se->registerInputManager(hitsin);
      i++;
    }

  Enable::CDB = true;
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  std::cout << "CDB tracking geometry file "<<geofile << std::endl;
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING = true;

  G4MAGNET::magfield_rescale = 1;
  TrackingInit();
  
  Tpc_HitUnpacking();
  
  Tpc_LaserEventIdentifying();

  TPC_LaminationClustering();

  /*
  TPC_LaserClustering();

  auto tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(0);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  tpcclusterizer->set_reject_event(G4TPC::REJECT_LASER_EVENTS);
  se->registerSubsystem(tpcclusterizer);
  */

  TString dstfile = theOutfile + "_dst.root";
  std::string dststring(dstfile.Data());
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("out", dststring);
  out->AddNode("Sync");
  out->AddNode("EventHeader");
  out->AddNode("LAMINATION_CLUSTER");
  out->AddNode("LaserEventInfo");
  se->registerOutputManager(out);

  se->skip(stepsize*index);
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
