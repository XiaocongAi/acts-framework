// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/CalibrationDecorator.hpp"

FW::Contextual::CalibrationDecorator::CalibrationDecorator(
    const FW::Contextual::CalibrationDecorator::Config& cfg,
    std::unique_ptr<const Acts::Logger>                 logger)
  : m_cfg(cfg), m_logger(std::move(logger))
{
  if (!m_cfg.randomNumberSvc) {
    throw std::invalid_argument("Missing random numbers tool");
  }
}

FW::ProcessCode
FW::Contextual::CalibrationDecorator::decorate(AlgorithmContext& context)
{
  // We need to lock the Decorator
  std::lock_guard<std::mutex> calibrationLock(m_calibrationMutex);

  // In which iov batch are we?
  unsigned int iov = context.eventNumber / m_cfg.iovSize;

  // Detect if we have a new calibration range
  auto calibStore_it = m_calibrationStore.find(iov);
  if (calibStore_it == m_calibrationStore.end()) {
    ACTS_INFO("New IOV detected at event " << context.eventNumber
                                           << ", emulate new calibration.");

    // Create an algorithm local random number generator
    RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

    // Create a calibration map
    std::map<const Acts::Surface*, CalibrationData> calibMap;
    for (auto& [surface, calibDataGen] : m_cfg.calibrationDataGenerator) {
      if (surface->associatedDetectorElement()) {
        auto calibData = calibDataGen(rng);
        // push the calibration data for this surface
        calibMap.emplace(surface, calibData);
      }
    }
    // push the new iov calibration map to the store
    m_calibrationStore.emplace(iov, std::move(calibMap));
  }

  // This creates a calib Context from the calibration map
  FW::Contextual::CalibrationContextType calibContext{
      m_calibrationStore.at(iov)};
  context.calibContext
      = std::make_any<FW::Contextual::CalibrationContextType>(calibContext);

  return ProcessCode::SUCCESS;
}
