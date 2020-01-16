// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/Validation/ProtoTrackClassification.hpp"

#include <boost/none.hpp>
#include <boost/optional.hpp>

namespace FW {

/// @brief struct for track finding/fitting result with
/// CombinatorialKalmanFilter(CKF)
/// It could consists of a multitrajectory, a vector of entry indices for
/// individual trajectories and a map of fitted parameters for each trajectory
///
/// @Todo Use a track proxy or helper to retrieve the detailed info, such as
/// number of measurments, holes, truth info etc.
struct CKFTrack
{
public:
  // Default constructor
  CKFTrack() = default;

  /// Constructor from found multitrajectory
  ///
  /// @param tTips The multiTrajectory entry points
  /// @param trajectory The found multiTrajectory
  CKFTrack(const std::vector<size_t>&                        tTips,
           const Acts::MultiTrajectory<Data::SimSourceLink>& trajectory)
    : m_trackTips(tTips), m_trajectory(trajectory)
  {
  }

  /// Constructor from fitted track parameter
  ///
  /// @param parameter The fitted track parameter
  CKFTrack(const std::map<size_t, Acts::BoundParameters>& parameters)
    : m_trackParameters(parameters)
  {
  }

  /// Constructor from found multitrajectory and fitted track parameter
  ///
  /// @param tTips The multiTrajectory entry points
  /// @param trajectory The found multiTrajectory
  /// @param parameter The fitted track parameter
  CKFTrack(const std::vector<size_t>&                        tTips,
           const Acts::MultiTrajectory<Data::SimSourceLink>& trajectory,
           const std::map<size_t, Acts::BoundParameters>&    parameters)
    : m_trackTips(tTips)
    , m_trajectory(trajectory)
    , m_trackParameters(parameters)
  {
  }

  /// Indicator for having found multitrajectory or not
  bool
  hasTrajectory() const
  {
    return m_trajectory ? true : false;
  }

  /// Indicator for having fitted track parameter or not for trajectory with
  /// given entry index
  bool
  hasTrackParameters(const size_t& entryIndex) const
  {
    if (not m_trackParameters.empty()) {
      auto it = m_trackParameters.find(entryIndex);
      if (it != m_trackParameters.end()) { return true; }
    }
    return false;
  }

  /// Get multitrajectory along with the entry points
  const std::pair<std::vector<size_t>,
                  Acts::MultiTrajectory<Data::SimSourceLink>>
  trajectory() const
  {
    if (m_trajectory and not m_trackTips.empty()) {
      return std::make_pair(m_trackTips, *m_trajectory);
    } else {
      throw std::runtime_error("No multi-trajectory!");
    };
  }

  /// Get fitted track parameter for trajectory with given entry index
  const Acts::BoundParameters&
  trackParameters(const size_t& entryIndex) const
  {
    if (not m_trackParameters.empty()) {
      auto it = m_trackParameters.find(entryIndex);
      if (it != m_trackParameters.end()) {
        return it->second;
      } else {
        throw std::runtime_error(
            "No fitted track parameter for trajectory with entry index = "
            + entryIndex);
      }
    } else {
      throw std::runtime_error(
          "No fitted track parameters for this multi-trajectory!");
    }
  }

  /// Get number of trajectories in the multitrajectory
  size_t
  numTrajectories() const
  {
    return m_trackTips.size();
  }

  /// Get number of track states that have measurements for trajectory with
  /// given entry index
  size_t
  numMeasurements(const size_t& entryIndex) const
  {
    size_t nMeasurements = 0;
    if (m_trajectory) {
      // Find if the multitrajectory has an entry index as given
      auto it = std::find_if(
          m_trackTips.begin(), m_trackTips.end(), [=](const size_t& trackTip) {
            return trackTip == entryIndex;
          });
      if (it != m_trackTips.end()) {
        (*m_trajectory).visitBackwards(entryIndex, [&](const auto& state) {
          if (state.hasUncalibrated()) { nMeasurements++; }
        });
      }
    }
    return nMeasurements;
  }

  /// Find the majority particle for trajectory with given entry index
  std::vector<ParticleHitCount>
  identifyMajorityParticle(const size_t& entryIndex) const
  {
    std::vector<ParticleHitCount> particleHitCount;

    if (m_trajectory) {
      (*m_trajectory).visitBackwards(entryIndex, [&](const auto& state) {
        // No truth info with non-measurement state
        if (not state.hasUncalibrated()) { return true; }
        // Find the truth particle associated with this state
        const auto& particle = state.uncalibrated().truthHit().particle;

        // Get the barcode
        auto particleId = particle.barcode();

        // Find if the particle already exists
        auto it = std::find_if(particleHitCount.begin(),
                               particleHitCount.end(),
                               [=](const ParticleHitCount& phc) {
                                 return phc.particleId == particleId;
                               });

        // Either increase count if we saw the particle before or add it
        if (it != particleHitCount.end()) {
          it->hitCount += 1;
        } else {
          particleHitCount.push_back({particleId, 1u});
        }
        return true;
      });
    }

    if (not particleHitCount.empty()) {
      // sort by hit count, i.e. majority particle first
      std::sort(particleHitCount.begin(),
                particleHitCount.end(),
                [](const ParticleHitCount& lhs, const ParticleHitCount& rhs) {
                  return lhs.hitCount < rhs.hitCount;
                });
    }

    return std::move(particleHitCount);
  }

private:
  // The optional found multitrajectory
  boost::optional<Acts::MultiTrajectory<Data::SimSourceLink>> m_trajectory{
      boost::none};

  // This is the indices of the 'tip' of the tracks stored in multitrajectory.
  std::vector<size_t> m_trackTips = {};

  // The parameters at the provided surface
  std::map<size_t, Acts::BoundParameters> m_trackParameters = {};
};

}  // namespace FW
