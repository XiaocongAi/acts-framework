#include "RandomNumbersAlgorithm.hpp"

#include <iostream>

#include "ACTFW/Random/RandomNumbersSvc.hpp"

FWE::RandomNumbersAlgorithm::RandomNumbersAlgorithm(
    const FWE::RandomNumbersAlgorithm::Config& cfg,
    std::unique_ptr<Acts::Logger>              logger)
  : FW::Algorithm(cfg, std::move(logger)), m_cfg(cfg)
{
}

FWE::RandomNumbersAlgorithm::~RandomNumbersAlgorithm()
{
}

/** Framework finalize mehtod */
FW::ProcessCode
FWE::RandomNumbersAlgorithm::initialize(std::shared_ptr<FW::WhiteBoard> jStore)
{
  // call the algorithm initialize for setting the stores
  if (FW::Algorithm::initialize(jStore) != FW::ProcessCode::SUCCESS) {
    ACTS_FATAL("Algorithm::initialize() did not succeed!");
    return FW::ProcessCode::SUCCESS;
  }
  ACTS_VERBOSE("initialize successful.");
  return FW::ProcessCode::SUCCESS;
}

/** Framework execode method */
FW::ProcessCode
FWE::RandomNumbersAlgorithm::execute(const FW::AlgorithmContext context) const
{
  // Create a random number generator
  FW::RandomNumbersSvc::Generator rng =
    m_cfg.randomNumbers->spawnGenerator(context);

  for (size_t idraw = 0; idraw < m_cfg.drawsPerEvent; ++idraw) {
    double gauss   = rng.draw(FW::Distribution::gauss);
    double uniform = rng.draw(FW::Distribution::uniform);
    double landau  = rng.draw(FW::Distribution::landau);
    double gamma   = rng.draw(FW::Distribution::gamma);

    ACTS_VERBOSE("Gauss   : " << gauss);
    ACTS_VERBOSE("Uniform : " << uniform);
    ACTS_VERBOSE("Landau  : " << landau);
    ACTS_VERBOSE("Gamma   : " << gamma);
  }
  return FW::ProcessCode::SUCCESS;
}

/** Framework finalize mehtod */
FW::ProcessCode
FWE::RandomNumbersAlgorithm::finalize()
{
  ACTS_VERBOSE("initialize successful.");
  return FW::ProcessCode::SUCCESS;
}
