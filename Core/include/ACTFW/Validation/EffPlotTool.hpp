// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <string>
#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/Utilities/Helpers.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class EffPlotTool
{
  using Identifier  = Data::SimSourceLink;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;

public:
  /// @brief The nested configuration struct
  struct Config
  {
    std::map<std::string, PlotHelpers::Binning> varBinning
        = {{"Eta", PlotHelpers::Binning("#eta", 50, -4, 4)},
           {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
           {"Pt", PlotHelpers::Binning("pT [GeV/c]", 100, 0, 1000)}};
  };

  /// @brief Nested Cache struct
  struct EffPlotCache
  {

    TEfficiency* trackeff_vs_eta;  ///< Tracking efficiency vs eta
    TEfficiency* trackeff_vs_phi;  ///< Tracking efficiency vs phi
    TEfficiency* trackeff_vs_pT;   ///< Tracking efficiency vs pT
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  EffPlotTool(const Config&        cfg,
              Acts::Logging::Level level = Acts::Logging::INFO);

  /// @brief book the efficiency plots
  /// @param effPlotCache the cache for efficiency plots
  void
  book(EffPlotCache& effPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void
  fill(EffPlotCache&            effPlotCache,
       const Data::SimParticle& truthParticle,
       bool                     status) const;

  /// @brief write the efficiency plots to file
  /// @param effPlotCache cache object for efficiency plots
  void
  write(const EffPlotCache& effPlotCache) const;

  /// @brief delete the efficiency plots
  /// @param effPlotCache cache object for efficiency plots
  void
  clear(EffPlotCache& effPlotCache) const;

private:
  Config                              m_cfg;     ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW
