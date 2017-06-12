///////////////////////////////////////////////////////////////////
// GenatinoRecording.cpp
///////////////////////////////////////////////////////////////////

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Writers/IWriterT.hpp"
#include "ACTFW/MaterialMapping/GeantinoRecording.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"

int
main()
{
  // job steering: the number of events
  size_t               nEvents      = 1000;
  size_t               nTracks      = 100;
  int                  randomSeed1  = 536235167;
  int                  randomSeed2  = 729237523;
  
  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  //
  // set up the geometry service 
  FWDD4hep::GeometryService::Config gsConfig("GeometryService",
                                              Acts::Logging::INFO);
  gsConfig.xmlFileName
      = "file:Examples/DD4hepExample/DD4hepDetector/compact/FCCTracker.xml";
  auto geometrySvc = std::make_shared<FWDD4hep::GeometryService>(gsConfig);
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = geometrySvc->trackingGeometry();

  // DD4Hep to Geant4 conversion
  //
  FWDD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4", Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc     = std::make_shared<FWDD4hepG4::DD4hepToG4Svc>(dgConfig);
  
  // --------------------------------------------------------------------------------
  // Geant4 JOB:
  // --------------------------------------------------------------------------------
  // set up the writer for 
  FWRoot::RootMaterialTrackWriter::Config g4WriterConfig("MaterialTrackWriter",
                                                            Acts::Logging::INFO);
  g4WriterConfig.fileName = "GeantMaterialTracks3.root";
  g4WriterConfig.treeName = "GeantMaterialTracks";
  auto g4TrackRecWriter
      = std::make_shared<FWRoot::RootMaterialTrackWriter>(g4WriterConfig);  
  
  // set up the algorithm writing out the material map
  FWA::GeantinoRecording::Config g4rConfig;
  g4rConfig.materialTrackWriter = g4TrackRecWriter;
  g4rConfig.geant4Service       = dd4hepToG4Svc;
  g4rConfig.tracksPerEvent      = nTracks;
  g4rConfig.seed1               = randomSeed1;
  g4rConfig.seed2               = randomSeed2;
  // create the geant4 algorithm
  auto g4rAlgorithm    = std::make_shared<FWA::GeantinoRecording>(
      g4rConfig, Acts::getDefaultLogger("GeantinoRecording", Acts::Logging::INFO));
  
  // Geant4 job - these can be many Geant4 jobs, indeed
  //
  // create the config object for the sequencer
  FW::Sequencer::Config g4SeqConfig;
  // now create the sequencer
  FW::Sequencer g4Sequencer(g4SeqConfig);
  // the writer is a service as it needs initialize, finalize
  g4Sequencer.addServices({g4TrackRecWriter});
  g4Sequencer.appendEventAlgorithms({g4rAlgorithm});
  // initialize loop
  g4Sequencer.initializeEventLoop();
  // run the loop
  g4Sequencer.processEventLoop(nEvents);
  // finalize loop
  g4Sequencer.finalizeEventLoop();
  // --------------------------------------------------------------------------------
  
}