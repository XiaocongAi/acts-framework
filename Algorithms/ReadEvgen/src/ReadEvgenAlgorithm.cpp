#include <iostream>
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/ReadEvgen/ReadEvgenAlgorithm.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FWA::ReadEvgenAlgorithm::ReadEvgenAlgorithm(
    const Config&                       cfg,
    std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(cfg), m_logger(std::move(logger))
{
}

FW::ProcessCode
FWA::ReadEvgenAlgorithm::skip(size_t nEvents)
{
  // there is a hard scatter evgen reader
  std::vector<Acts::ParticleProperties> skipParticles;  
  if (m_cfg.hardscatterParticleReader && 
      m_cfg.hardscatterParticleReader->read(skipParticles, nEvents)
      == FW::ProcessCode::ABORT){
    // error and abort
    ACTS_ERROR("Could not skip " << nEvents << ". Aborting.");
    return FW::ProcessCode::ABORT;    
  }
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FWA::ReadEvgenAlgorithm::read(const FW::AlgorithmContext context) const
{

  // Retrieve relevant information from the execution context
  size_t eventNumber = context.eventContext->eventNumber;
  auto   eventStore  = context.eventContext->eventStore;

  ACTS_DEBUG("Reading in genertated event info for event no. " << eventNumber);

  // Create a random number generator
  FW::RandomNumbersSvc::Generator rngPileup
      = m_cfg.pileupRandomNumbers->spawnGenerator(context);

  FW::RandomNumbersSvc::Generator rngVertexT
      = m_cfg.pileupVertexDistT->spawnGenerator(context);

  FW::RandomNumbersSvc::Generator rngVertexZ
      = m_cfg.pileupVertexDistZ->spawnGenerator(context);

  // prepare the output vector
  std::vector<Acts::ParticleProperties>* eventParticles
      = new std::vector<Acts::ParticleProperties>;

  // get the hard scatter if you have it
  std::vector<Acts::ParticleProperties> hardscatterParticles = {};
  if (m_cfg.hardscatterParticleReader && 
      m_cfg.hardscatterParticleReader->read(hardscatterParticles) 
      == FW::ProcessCode::ABORT){
      ACTS_ERROR("Could not read hard scatter event. Aborting.");
      return FW::ProcessCode::ABORT;
  }
  ACTS_VERBOSE("- [HS X] number of hard scatter particles   : "
               << (hardscatterParticles.size() > 0 ? 1 : 0));

  // generate the number of pileup events
  size_t nPileUpEvents = m_cfg.pileupRandomNumbers
      ? size_t(rngPileup.drawPoisson())
      : 0;

  ACTS_VERBOSE("- [PU X] number of in-time pileup events : " << nPileUpEvents);

  // reserve a lot
  eventParticles->reserve((nPileUpEvents)*hardscatterParticles.size() * 2);

  //
  // reserve quite a lot of space
  double vertexX = rngVertexT.drawGauss();
  double vertexY = rngVertexT.drawGauss();
  double vertexZ = rngVertexZ.drawGauss();

  Acts::Vector3D vertex(vertexX, vertexY, vertexZ);

  // fill in the particles
  barcode_type pCounter = 0;
  for (auto& hsParticle : hardscatterParticles) {
    // shift the particle by the vertex
    hsParticle.shift(vertex);
      hsParticle.assign(m_cfg.barcodeSvc->generate(0,pCounter++));
    // now push-back
    eventParticles->push_back(hsParticle);
  }

  // loop over the pile-up vertices
  for (size_t ipue = 0; ipue < nPileUpEvents; ++ipue) {
    // reserve quite a lot of space
    double         puVertexX = rngVertexT.drawGauss();
    double         puVertexY = rngVertexT.drawGauss();
    double         puVertexZ = rngVertexZ.drawGauss();
    Acts::Vector3D puVertex(puVertexX, puVertexY, puVertexZ);
    // get the vertices per pileup event
    std::vector<Acts::ParticleProperties> pileupPartiles = {};
    if (m_cfg.pileupParticleReader && 
          m_cfg.pileupParticleReader->read(pileupPartiles) 
          == FW::ProcessCode::ABORT){
          ACTS_ERROR("Could not read pile up event " << ipue << ". Aborting.");
          return FW::ProcessCode::ABORT;
    }
    pCounter = 0;
    ACTS_VERBOSE("- [PU " << ipue << "] number of pile-up particles : "
                          << pileupPartiles.size()
                          << " - with z vertex position: "
                          << puVertexZ);
    // loop over pileupParicles
    for (auto& puParticle : pileupPartiles) {
      // shift to the pile-up vertex
      puParticle.shift(puVertex);
      puParticle.assign(m_cfg.barcodeSvc->generate(ipue+1,pCounter++));
      // now store the particle
      eventParticles->push_back(puParticle);
    }
  }

  // write to file if you have
  if (m_cfg.particleWriter
      && m_cfg.particleWriter->write(*eventParticles)
          == FW::ProcessCode::ABORT) {
    ACTS_WARNING(
        "Could not write colleciton of particles to file. Aborting.");
    return FW::ProcessCode::ABORT;
  }

  // write to the EventStore
  if (eventStore
      && eventStore->writeT(eventParticles, m_cfg.evgenParticlesCollection)
          == FW::ProcessCode::ABORT) {
    ACTS_WARNING(
        "Could not write colleciton of process vertices to event store.");
    return FW::ProcessCode::ABORT;
  }

  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FWA::ReadEvgenAlgorithm::initialize(std::shared_ptr<FW::WhiteBoard> jStore)
{
  m_cfg.jBoard = jStore;
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FWA::ReadEvgenAlgorithm::finalize()
{
  return FW::ProcessCode::SUCCESS;
}
