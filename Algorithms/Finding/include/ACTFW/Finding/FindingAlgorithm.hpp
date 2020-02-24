// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <vector>

// this include is a work-around. its missing in the KalmanFitter header
// @TODO: include map in KalmanFitter
#include <map>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/TrackFinder/CKFSourceLinkSelector.hpp>
#include <Acts/TrackFinder/CombinatorialKalmanFilter.hpp>

#include "ACTFW/EventData/CKFTrack.hpp"
#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"

namespace FW {

class FindingAlgorithm final : public BareAlgorithm
{
public:
  using FinderResult
      = Acts::Result<Acts::CombinatorialKalmanFilterResult<Data::SimSourceLink>>;
  /// Fit function that takes input measurements, initial track state and fitter
  /// options and returns some fit-specific result.
  using FinderFunction = std::function<FinderResult(
      const std::vector<Data::SimSourceLink>&,
      const TrackParameters&,
      const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>&)>;

  /// Create the fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static FinderFunction
  makeFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant                        magneticField,
      Acts::Logging::Level                          lvl);

  struct Config
  {
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased fitter function.
    FinderFunction find;

    Acts::CKFSourceLinkSelector::Config slsCfg;
  };

  /// Constructor of the Finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  FindingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  FW::ProcessCode
  execute(const FW::AlgorithmContext& ctx) const final override;

private:
  Config m_cfg;
};

}  // namespace FW
