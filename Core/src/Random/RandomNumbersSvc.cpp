//
//  RandomNumbersSvc.cpp
//  ACTFW
//
//  Created by Andreas Salzburger on 17/05/16.
//
//

#include "ACTFW/Random/RandomNumbersSvc.hpp"

FW::RandomNumbersSvc::RandomNumbersSvc(
    const Config&                       cfg,
    std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(cfg), m_logger(std::move(logger))
{
}

std::string
FW::RandomNumbersSvc::name() const
{
  return "RandomNumbersSvc";
}

FW::ProcessCode
FW::RandomNumbersSvc::initialize()
{
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::RandomNumbersSvc::finalize()
{
  return FW::ProcessCode::SUCCESS;
}

FW::RandomEngine
FW::RandomNumbersSvc::spawnGenerator(const AlgorithmContext& context) const
{
  const auto         eventContext = context.eventContext;
  const unsigned int generatorID
      = context.algorithmNumber * eventContext->jobContext->eventCount
      + eventContext->eventNumber;
  return RandomEngine(m_cfg.seed + generatorID);
}
