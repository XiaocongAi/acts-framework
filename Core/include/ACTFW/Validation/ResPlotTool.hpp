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
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace FW {

// Tools to make hists to show residual (smoothed_parameter - truth_parameter)
// and pull (smoothed_parameter - truth_parameter)/smoothed_paramter_error
// distributions of track parameters.
// For the moment, the entry is state-wise.
class ResPlotTool
{
  using Identifier  = Data::SimSourceLink;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;
  using ParVector_t = typename Acts::BoundParameters::ParVector_t;

public:
  /// @brief Nested configuration struct
  struct Config
  {
    /// parameter sets to do plots
    std::vector<std::string> paramNames
        = {"LOC_0", "LOC_1", "PHI", "THETA", "QOP", "TIME"};
    /// Binning info for variables
    std::map<std::string, PlotHelpers::Binning> varBinning
        = {{"Eta", PlotHelpers::Binning("#eta", 50, -4, 4)},
           {"R", PlotHelpers::Binning("r [mm]", 100, 0, 1100)},
           {"Z", PlotHelpers::Binning("z [mm]", 100, -3500, 3500)},
           {"Residual", PlotHelpers::Binning("residual", 100, -0.05, 0.05)},
           {"Pull", PlotHelpers::Binning("pull", 100, -5, 5)}};
  };

  /// @brief Nested Cache struct
  struct ResPlotCache
  {

    std::map<std::string, TH1F*> res;  ///< Residual distribution

    std::map<std::string, TH2F*> res_vs_eta;  ///< Residual vs eta scatter plot
    std::map<std::string, TH1F*>
        resmean_vs_eta;  ///< Residual mean vs eta distribution
    std::map<std::string, TH1F*>
        reswidth_vs_eta;  ///< Residual width vs eta distribution

    std::map<std::string, TH2F*>
        res_vs_r;  ///< Residual vs global r scatter plot
    std::map<std::string, TH1F*>
        resmean_vs_r;  ///< Residual mean vs global r distribution
    std::map<std::string, TH1F*>
        reswidth_vs_r;  ///< Residual width vs global r distribution

    std::map<std::string, TH2F*>
        res_vs_z;  ///< Residual vs global z scatter plot
    std::map<std::string, TH1F*>
        resmean_vs_z;  ///< Residual mean vs global z distribution
    std::map<std::string, TH1F*>
        reswidth_vs_z;  ///< Residual width vs global z distribution

    std::map<std::string, TH1F*> pull;  ///< Pull distribution

    std::map<std::string, TH2F*> pull_vs_eta;  ///< Pull vs eta scatter plot
    std::map<std::string, TH1F*>
        pullmean_vs_eta;  ///< Pull mean vs eta distribution
    std::map<std::string, TH1F*>
        pullwidth_vs_eta;  ///< Pull width vs eta distribution

    std::map<std::string, TH2F*> pull_vs_r;  ///< Pull vs global r scatter plot
    std::map<std::string, TH1F*>
        pullmean_vs_r;  ///< Pull mean vs global r distribution
    std::map<std::string, TH1F*>
        pullwidth_vs_r;  ///< Pull width vs global r distribution

    std::map<std::string, TH2F*> pull_vs_z;  ///< Pull vs global z scatter plot
    std::map<std::string, TH1F*>
        pullmean_vs_z;  ///< Pull mean vs global z distribution
    std::map<std::string, TH1F*>
        pullwidth_vs_z;  ///< Pull width vs global z distribution
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  ResPlotTool(const Config&        cfg,
              Acts::Logging::Level level = Acts::Logging::INFO);

  /// @brief book the histograms
  /// @param resPlotCache the cache for residual/pull histograms
  void
  book(ResPlotCache& resPlotCache) const;

  /// @brief fill the histograms
  /// @param resPlotCache the cache for residual/pull histograms
  /// @param track a vector of trackState for the moment
  /// @param truth a vector of truthHit
  void
  fill(ResPlotCache&                                               resPlotCache,
       const Acts::GeometryContext&                                gctx,
       const std::pair<size_t, Acts::MultiTrajectory<Identifier>>& trajectory)
      const;

  /// @brief extract the details of the residual/pull plots and fill details
  /// into separate histograms
  /// @param resPlotCache the cache object for residual/pull histograms
  void
  refinement(ResPlotCache& resPlotCache) const;

  /// @brief write the histograms to output file
  /// @param resPlotCache the cache object for residual/pull histograms
  void
  write(const ResPlotCache& resPlotCache) const;

  /// @brief delele the histograms
  /// @param resPlotCache the cache object for residual/pull histograms
  void
  clear(ResPlotCache& resPlotCache) const;

private:
  Config                              m_cfg;     ///< The config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW
