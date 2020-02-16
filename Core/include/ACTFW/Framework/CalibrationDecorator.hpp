// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <mutex>
#include <vector>
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {

namespace Contextual {
  using CalibrationData      = std::vector<std::pair<Acts::ParID_t, double>>;
  using CalibDataGenFunction = std::function<CalibrationData(RandomEngine&)>;
  ///
  /// @class MisCalibrationGenerator
  ///
  struct MisCalibrationGenerator
  {
    // The sigma of miscalibration for local measurements
    std::vector<std::pair<Acts::ParID_t, double>> parSigma;

    // The function to generate miscalibration for local measurements
    std::vector<std::pair<Acts::ParID_t, double>>
    operator()(RandomEngine& rng) const
    {
      std::vector<std::pair<Acts::ParID_t, double>> miscalib;
      for (const auto& [parID, sigma] : parSigma) {
        miscalib.push_back(
            std::pair(parID, std::normal_distribution<double>(0., sigma)(rng)));
      }
      return std::move(miscalib);
    }
  };

  ///
  /// @class CalibrationContextType
  ///
  struct CalibrationContextType
  {
    // The calibration map
    std::map<const Acts::Surface*, CalibrationData> calibrationMap;
  };

  /// @brief
  ///
  ///
  class CalibrationDecorator : public IContextDecorator
  {
  public:
    /// @brief nested configuration struct
    struct Config
    {
      /// Miscalibration generator
      std::map<const Acts::Surface*, CalibDataGenFunction>
          calibrationDataGenerator;

      /// Calibration frequency - every X events
      unsigned int iovSize = 100;

      /// Random numbers tool
      std::shared_ptr<RandomNumbers> randomNumberSvc = nullptr;
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param logger The logging framework
    CalibrationDecorator(const Config&                       cfg,
                         std::unique_ptr<const Acts::Logger> logger
                         = Acts::getDefaultLogger("CalibrationDecorator",
                                                  Acts::Logging::INFO));

    /// Virtual destructor
    virtual ~CalibrationDecorator() = default;

    /// @brief decorates (adds, modifies) the AlgorithmContext
    /// with calibration data per event
    ///
    /// @param context the bare (or at least non-const) Event context
    ProcessCode
    decorate(AlgorithmContext& context) final override;

    /// @brief decorator name() for screen output
    const std::string&
    name() const final override
    {
      return m_name;
    }

  private:
    Config                              m_cfg;     ///< the configuration class
    std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
    std::string                         m_name = "CalibrationDecorator";

    ///< protect multiple calibration to be loaded at once
    std::mutex m_calibrationMutex;

    /// Calibration store
    //@TODO: store the calibration data in detector element
    std::map<size_t, std::map<const Acts::Surface*, CalibrationData>>
        m_calibrationStore;

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_logger;
    }
  };
}  // namespace Contextual

}  // namespace FW
