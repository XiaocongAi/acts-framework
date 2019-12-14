// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Fitting/FittingAlgorithm.hpp"

#include <stdexcept>
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::FittingAlgorithm::FittingAlgorithm(Config cfg, Acts::Logging::Level level)
  : FW::BareAlgorithm("FittingAlgorithm", level), m_cfg(std::move(cfg))
{
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

FW::ProcessCode
FW::FittingAlgorithm::execute(const FW::AlgorithmContext& ctx) const
{

  // Read input data
  const auto sourceLinks
      = ctx.eventStore.get<SimSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto protoTracks
      = ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  const auto initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters");
    return ProcessCode::ABORT;
  }

  // Prepare the output data with MultiTrajectory
  TrajectoryContainer trajectories;
  trajectories.reserve(protoTracks.size());

  // Perform the fit for each input track
  std::vector<Data::SimSourceLink> trackSourceLinks;
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack    = protoTracks[itrack];
    const auto& initialParams = initialParameters[itrack];

    // We can have empty tracks which must give empty fit results
    if (protoTrack.empty()) {
      trajectories.push_back({});
      continue;
    }
    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links via their indices from the container
    for (auto hitIndex : protoTrack) {
      auto sourceLink = sourceLinks.nth(hitIndex);
      if (sourceLink == sourceLinks.end()) {
        ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
      trackSourceLinks.push_back(*sourceLink);
    }

    // Set the target surface
    const Acts::Surface* rSurface = &initialParams.referenceSurface();

    // Set the KalmanFitter options
    Acts::KalmanFitterOptions kfOptions(
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext, rSurface);

    ACTS_DEBUG("Invoke fitter");
    auto result = m_cfg.fit(trackSourceLinks, initialParams, kfOptions);
    if (result.ok()) {
      // Get the fit output object
      const auto& fitOutput = result.value();
      trajectories.push_back(
          std::make_pair(fitOutput.trackTip, fitOutput.fittedStates));

      // fitted parameters on the reference surface are only used for debug
      // output for now. only the fit trajectory is actually stored.
      // TODO store single track parameters in separate container.
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.get();
        ACTS_DEBUG("Fitted paramemeters for track " << itrack);
        ACTS_DEBUG("  position: " << params.position().transpose());
        ACTS_DEBUG("  momentum: " << params.momentum().transpose());
      }
    } else {
      ACTS_WARNING("Fit failed for track " << itrack << " with error"
                                           << result.error());
      //      trajectories.push_back({});
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
