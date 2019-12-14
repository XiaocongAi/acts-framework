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

namespace FW {

namespace TruthFitTrackStatus {
  using Status = std::bitset<8>;

  constexpr Status Trajectory{1 << 0};
  constexpr Status Measurements{1 << 1};
  constexpr Status TrackParameter{1 << 2};
}  // namespace TruthFitTrackStatus

/// @brief struct for truth fitting result
///
/// @note A truth particle must be provided regardless of fitting status
/// for investigation of track efficiency etc.
///
struct TruthFitTrack
{
public:
  /// Deleted default constructor
  TruthFitTrack() = delete;

  /// Constructor from only truth particle
  ///
  /// @param tParticle The truth particle associated to this trajectory
  TruthFitTrack(const Data::SimParticle& tParticle) : m_truthParticle(tParticle)
  {
  }

  /// Constructor from truth particle and fitted states
  ///
  /// @param tParticle The truth particle associated to this trajectory
  /// @param tTip The fitted multiTrajectory entry point
  /// @param trajectory The fitted multiTrajectory
  TruthFitTrack(const Data::SimParticle&                          tParticle,
                size_t                                            tTip,
                const Acts::MultiTrajectory<Data::SimSourceLink>& trajectory)
    : m_truthParticle(tParticle), m_trackTip(tTip), m_trajectory(trajectory)
  {
  }

  /// Constructor from truth particle, fitted states and fitted track parameter
  ///
  /// @param tParticle The truth particle associated to this trajectory
  /// @param tTip The fitted multiTrajectory entry point
  /// @param trajectory The fitted multiTrajectory
  /// @param parameter The fitted track parameter
  TruthFitTrack(const Data::SimParticle&                          tParticle,
                size_t                                            tTip,
                const Acts::MultiTrajectory<Data::SimSourceLink>& trajectory,
                const Acts::BoundParameters&                      parameter)
    : m_truthParticle(tParticle)
    , m_trackTip(tTip)
    , m_trajectory(trajectory)
    , m_trackParameter(parameter)
  {
  }

  /// Get truth particle
  const Data::SimParticle&
  truthParticle() const
  {
    return m_truthParticle;
  }

  /// Get trajectory along with the entry point
  const std::pair<size_t, Acts::MultiTrajectory<Data::SimSourceLink>>
  trajectory() const
  {
    if (m_trajectory) {
      return std::make_pair(m_trackTip, *m_trajectory);
    } else {
      throw std::runtime_error("No fitted states on this trajectory!");
    };
  }

  /// Get fitted track parameter
  const Acts::BoundParameters
  parameters() const
  {
    if (m_trackParameter) {
      return *m_trackParameter;
    } else {
      throw std::runtime_error(
          "No fitted track parameter for this trajectory!");
    }
  }

  /// Get number of track states
  size_t
  numStates() const
  {
    size_t nStates = 0;
    if (m_trajectory) {
      (*m_trajectory).visitBackwards(m_trackTip, [&](const auto& state) {
        nStates++;
      });
    }
    return nStates;
  }

  /// Get number of track states that have measurements
  size_t
  numMeasurements() const
  {
    size_t nMeasurements = 0;
    if (m_trajectory) {
      (*m_trajectory).visitBackwards(m_trackTip, [&](const auto& state) {
        if (state.hasUncalibrated()) { nMeasurements++; }
      });
    }
    return nMeasurements;
  }

  /// Get track status
  TruthFitTrackStatus::Status
  status() const
  {
    TruthFitTrackStatus::Status status{0};
    if (m_trajectory) { status |= TruthFitTrackStatus::Trajectory; }
    if (m_trackParameter) { status |= TruthFitTrackStatus::TrackParameter; }
    if (numMeasurements() > 0) { status |= TruthFitTrackStatus::Measurements; }
    return std::move(status);
  }

private:
  // Truth particle associated to the trajectory
  const Data::SimParticle m_truthParticle;

  // The optional fitted multitrajectory
  boost::optional<Acts::MultiTrajectory<Data::SimSourceLink>> m_trajectory{
      boost::none};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  size_t m_trackTip = SIZE_MAX;

  // The optional Parameters at the provided surface
  boost::optional<Acts::BoundParameters> m_trackParameter{boost::none};
};

}  // namespace FW
