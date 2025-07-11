#ifndef MACRO_TRKRQA_C
#define MACRO_TRKRQA_C

#include <GlobalVariables.C>
#include <G4_TrkrVariables.C>
#include <Trkr_TruthTables.C>
#include <QA.C>
#include <fun4all/Fun4AllServer.h>
#include <simqa_modules/QAG4SimulationMvtx.h>
#include <simqa_modules/QAG4SimulationIntt.h>
#include <simqa_modules/QAG4SimulationTpc.h>
#include <simqa_modules/QAG4SimulationMicromegas.h>
#include <simqa_modules/QAG4SimulationTracking.h>
#include <simqa_modules/QAG4SimulationUpsilon.h>
#include <simqa_modules/QAG4SimulationVertex.h>
#include <simqa_modules/QAG4SimulationDistortions.h>
#include <trackingqa/VertexQA.h>

R__LOAD_LIBRARY(libsimqa_modules.so)
R__LOAD_LIBRARY(libtrackingqa.so)

void Mvtx_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::MVTX_VERBOSITY);

  Fun4AllServer* se = Fun4AllServer::instance();
  QAG4SimulationMvtx* qa = new QAG4SimulationMvtx;
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);
}

void Intt_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::INTT_VERBOSITY);

  Fun4AllServer* se = Fun4AllServer::instance();
  QAG4SimulationIntt* qa = new QAG4SimulationIntt;
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);
}

void TPC_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::TPC_VERBOSITY);

  Fun4AllServer* se = Fun4AllServer::instance();
  QAG4SimulationTpc * qa =  new QAG4SimulationTpc;
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);
}

void Micromegas_QA()
{
  auto se = Fun4AllServer::instance();
  auto qa_mm = new QAG4SimulationMicromegas;
  qa_mm->Verbosity(Enable::QA_VERBOSITY);
  se->registerSubsystem(qa_mm);
}


void Tracking_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::TRACKING_VERBOSITY);

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

std::cout<<"Hey this is my Tracking_QA"<<std::endl;
  VertexQA* qa_track = new VertexQA();
  qa_track->Verbosity(verbosity);
  se->registerSubsystem(qa_track);

  build_truthreco_tables();

  QAG4SimulationTracking* qa = new QAG4SimulationTracking();
  //  qa->addEmbeddingID(2);
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);

  QAG4SimulationVertex* qa2 = new QAG4SimulationVertex();
  // qa2->addEmbeddingID(2);
  qa2->Verbosity(verbosity);
  se->registerSubsystem(qa2);

  //  Acts Kalman Filter vertex finder
  //=================================
  QAG4SimulationVertex* qav = new QAG4SimulationVertex();
  // qav->addEmbeddingID(2);
  qav->Verbosity(verbosity);
  qav->setVertexMapName("SvtxVertexMapActs");
  se->registerSubsystem(qav);

  if (Input::UPSILON)
  {
    QAG4SimulationUpsilon* qa = new QAG4SimulationUpsilon();

    for (int id : Input::UPSILON_EmbedIds)
    {
      qa->addEmbeddingID(id);
    }
    se->registerSubsystem(qa);
  }
}

void Distortions_QA()
{
   int verbosity = std::max(Enable::QA_VERBOSITY, Enable::TRACKING_VERBOSITY);

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();
  
  std::cout<<"This is my Distortions_QA()"<<std::endl;
  auto qa = new QAG4SimulationDistortions();
  qa->Verbosity(verbosity);
  qa->set_trackmap_name("SvtxSiliconMMTrackMap");
  qa->disableAverageCorr();
  se->registerSubsystem(qa);
}

#endif
