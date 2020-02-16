// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/GenericDetector/GenericDetector.hpp"

#include <boost/program_options.hpp>
#include "ACTFW/Framework/CalibrationDecorator.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/GenericDetector/BuildGenericDetector.hpp"
#include "ACTFW/GenericDetector/GenericDetectorElement.hpp"
#include "ACTFW/GenericDetector/GenericDetectorOptions.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

void
GenericDetector::addOptions(
    boost::program_options::options_description& opt) const
{
  FW::Options::addGenericGeometryOptions(opt);
  opt.add_options()(
      "calib-seed",
      boost::program_options::value<size_t>()->default_value(1324354657),
      "Seed for the decorator random numbers.")(
      "calib-iovsize",
      boost::program_options::value<size_t>()->default_value(100),
      "Size of a valid IOV.")(
      "calib-flushsize",
      boost::program_options::value<size_t>()->default_value(200),
      "Span until garbage collection is active.")(
      "miscalib-sigma",
      boost::program_options::value<double>()->default_value(10.),
      "Sigma of the measurement miscalibration in [um]")(
      "calib-loglevel",
      boost::program_options::value<size_t>()->default_value(3),
      "Output log level of the calibration decorator.");
}

auto
GenericDetector::finalize(
    const boost::program_options::variables_map&    vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators>
{
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;

  auto buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  Acts::Logging::Level surfaceLogLevel
      = Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel
      = Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel
      = Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  bool buildProto
      = (vm["mat-input-type"].template as<std::string>() == "proto");

  /// Return the generic detector
  TrackingGeometryPtr gGeometry
      = FW::Generic::buildDetector<DetectorElement>(nominalContext,
                                                    detectorStore,
                                                    buildLevel,
                                                    std::move(mdecorator),
                                                    buildProto,
                                                    surfaceLogLevel,
                                                    layerLogLevel,
                                                    volumeLogLevel);
  Acts::Logging::Level decoratorLogLevel
      = Acts::Logging::Level(vm["calib-loglevel"].template as<size_t>());

  // Let's create a reandom number service
  FW::RandomNumbers::Config randomNumberConfig;
  randomNumberConfig.seed = vm["calib-seed"].template as<size_t>();
  auto randomNumberSvc
      = std::make_shared<FW::RandomNumbers>(randomNumberConfig);

  // The miscalibration generators
  std::map<const Acts::Surface*, FW::Contextual::CalibDataGenFunction>
      misCalibGenerator;
  for (const auto& layer : detectorStore) {
    for (const auto& det : layer) {
      // TODO: consider the dimension of measurements
      FW::Contextual::MisCalibrationGenerator miscalibGen{
          {{Acts::ParDef::eLOC_0, vm["miscalib-sigma"].template as<double>()},
           {Acts::ParDef::eLOC_1, vm["miscalib-sigma"].template as<double>()}}};
      misCalibGenerator.emplace(&det->surface(), std::move(miscalibGen));
    }
  }

  // Calibration decorator service
  FW::Contextual::CalibrationDecorator::Config calibDecoConfig;
  calibDecoConfig.calibrationDataGenerator = std::move(misCalibGenerator);
  calibDecoConfig.iovSize         = vm["calib-iovsize"].template as<size_t>();
  calibDecoConfig.randomNumberSvc = randomNumberSvc;

  // Make the context decorators
  ContextDecorators gContextDeocrators
      = {std::make_shared<FW::Contextual::CalibrationDecorator>(
          calibDecoConfig,
          Acts::getDefaultLogger("CalibrationDecorator", decoratorLogLevel))};

  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(gGeometry), std::move(gContextDeocrators));
}
