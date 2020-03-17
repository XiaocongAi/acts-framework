// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Fitting/FittingAlgorithm.hpp"

#include <stdexcept>
#include <tbb/tbb.h>
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

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
  // automatically determine the number of concurrent threads to use
  if (m_cfg.numFittingThreads < 0) {
    m_cfg.numFittingThreads = tbb::task_scheduler_init::default_num_threads();
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

  // Synchronize the access to the fitting results (trajectories)
  tbb::queuing_mutex trajectoriesMutex;

  auto addTrajectory = [&](TruthFitTrack&& track) {
    // Synchronize writes
    tbb::queuing_mutex::scoped_lock lock(trajectoriesMutex);
    trajectories.push_back(std::move(track));
  };

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3D{0., 0., 0.});

  // Setup a task arena for the the parallel loop (because this loop is
  // imbricated in the Sequencer parallel loop) to ensure: (a) better execution
  // time and (b) lower memory footprint
  tbb::task_arena arena(m_cfg.numFittingThreads);
  ACTS_INFO("Starting tracks loop with " << m_cfg.numFittingThreads
                                         << " threads");

  arena.execute([&]() {
    // Perform the fit for each input track
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, protoTracks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
          for (auto itrack = r.begin(); itrack != r.end(); ++itrack) {

            // The list of hits and the initial start parameters
            const auto& protoTrack    = protoTracks[itrack];
            const auto& initialParams = initialParameters[itrack];

            // We can have empty tracks which must give empty fit results
            if (protoTrack.empty()) {
              addTrajectory(TruthFitTrack());
              ACTS_WARNING("Empty track " << itrack << " found.");
              continue;
            }

            // Clear & reserve the right size
            std::vector<Data::SimSourceLink> trackSourceLinks;
            trackSourceLinks.clear();
            trackSourceLinks.reserve(protoTrack.size());

            // Fill the source links via their indices from the container
            for (auto hitIndex : protoTrack) {
              auto sourceLink = sourceLinks.nth(hitIndex);
              if (sourceLink == sourceLinks.end()) {
                ACTS_FATAL("Proto track " << itrack
                                          << " contains invalid hit index"
                                          << hitIndex);
                return ProcessCode::ABORT;
              }
              trackSourceLinks.push_back(*sourceLink);
            }

            // Set the KalmanFitter options
            Acts::KalmanFitterOptions kfOptions(ctx.geoContext,
                                                ctx.magFieldContext,
                                                ctx.calibContext,
                                                &(*pSurface));
            ACTS_DEBUG("Invoke fitter");
            auto result = m_cfg.fit(trackSourceLinks, initialParams, kfOptions);
            if (result.ok()) {
              // Get the fit output object
              const auto& fitOutput = result.value();
              if (fitOutput.fittedParameters) {
                const auto& params = fitOutput.fittedParameters.value();
                ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
                ACTS_VERBOSE("  position: " << params.position().transpose());
                ACTS_VERBOSE("  momentum: " << params.momentum().transpose());
                // Construct a truth fit track using trajectory and
                // track parameter
                addTrajectory(TruthFitTrack(fitOutput.trackTip,
                                            std::move(fitOutput.fittedStates),
                                            std::move(params)));
              } else {
                ACTS_DEBUG("No fitted paramemeters for track " << itrack);
                // Construct a truth fit track using trajectory
                addTrajectory(TruthFitTrack(fitOutput.trackTip,
                                            std::move(fitOutput.fittedStates)));
              }
            } else {
              ACTS_WARNING("Fit failed for track " << itrack << " with error"
                                                   << result.error());
              // Fit failed, but still create a empty truth fit track
              addTrajectory(TruthFitTrack());
            }
          }  // end for
          return FW::ProcessCode::SUCCESS;
        }  // end parallel_for
    );
  });  // end task arena

  if (m_cfg.sortedTracks) {
    // Make sure that the trajectories are in the same order as provided by the
    // sequential impl (ie sorted by particle id)
    std::sort(trajectories.begin(),
              trajectories.end(),
              [](const TruthFitTrack& t1, const TruthFitTrack& t2) -> bool {
                return t1.identifyMajorityParticle().front().particleId
                    < t2.identifyMajorityParticle().front().particleId;
              });
  }
  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
