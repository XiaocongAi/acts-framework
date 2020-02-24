// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdexcept>

#include "ACTFW/EventData/CKFTrack.hpp"
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Finding/FindingAlgorithm.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::FindingAlgorithm::FindingAlgorithm(Config cfg, Acts::Logging::Level level)
  : FW::BareAlgorithm("FindingAlgorithm", level), m_cfg(std::move(cfg))
{
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
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
FW::FindingAlgorithm::execute(const FW::AlgorithmContext& ctx) const
{

  // Read input data
  const auto sourceLinks
      = ctx.eventStore.get<SimSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Prepare the output data with MultiTrajectory
  CKFTrajectoryContainer trajectories;
  trajectories.reserve(initialParameters.size());

  // Perform the fit for each input track
  std::vector<Data::SimSourceLink> trackSourceLinks;
  for (std::size_t itrack = 0; itrack < initialParameters.size(); ++itrack) {
    const auto& initialParams = initialParameters[itrack];

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(initialParameters.size());

    for (const auto& sl : sourceLinks) { trackSourceLinks.push_back(sl); }

    // Set the target surface
    const Acts::Surface* rSurface = &initialParams.referenceSurface();

    // Set the KalmanFitter options
    using SourceLinkSelectorType = typename Acts::CKFSourceLinkSelector;
    using SourceLinkSelectorConfigType =
        typename SourceLinkSelectorType::Config;

    Acts::CombinatorialKalmanFilterOptions<SourceLinkSelectorType> tfOptions(
        ctx.geoContext,
        ctx.magFieldContext,
        ctx.calibContext,
        m_cfg.slsCfg,
        rSurface);

    ACTS_DEBUG("Invoke FindingAlgorithm");
    auto result = m_cfg.find(trackSourceLinks, initialParams, tfOptions);
    if (result.ok()) {
      // Get the actual TrackFinderResult object
      const auto& fitOutput = result.value();

      if (!fitOutput.fittedParameters.empty()) {
        // const auto& params = fitOutput.fittedParameters.get();
        // ACTS_VERBOSE("Fitted parameters for track " << itrack);
        // ACTS_VERBOSE("  position: " << params.position().transpose());
        // ACTS_VERBOSE("  momentum: " << params.momentum().transpose());

        // Construct a track using trajectory and
        // track parameter
        trajectories.emplace_back(fitOutput.trackTips,
                                  fitOutput.fittedStates,
                                  fitOutput.fittedParameters);
      } else {
        ACTS_DEBUG("No fitted parameters for track " << itrack);
        // Construct a track using trajectory

        trajectories.emplace_back(fitOutput.trackTips, fitOutput.fittedStates);
      }
    } else {
      ACTS_WARNING("Fit failed for track " << itrack << " with error"
                                           << result.error());
      // Fit failed, but still create a empty truth fit track
      trajectories.push_back(CKFTrack());
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
