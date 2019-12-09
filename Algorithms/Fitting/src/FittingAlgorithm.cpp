// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Fitting/FittingAlgorithm.hpp"

#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

FW::FittingAlgorithm::FittingAlgorithm(Config cfg, Acts::Logging::Level level)
  : FW::BareAlgorithm("FittingAlgorithm", level), m_cfg(std::move(cfg))
{
  if (m_cfg.simulatedHitCollection.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  } else if (m_cfg.simulatedEventCollection.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  } else if (m_cfg.trackCollection.empty()) {
    throw std::invalid_argument("Missing output track collection");
  }
}

FW::ProcessCode
FW::FittingAlgorithm::execute(const FW::AlgorithmContext& context) const
{
  // Create a random number generator
  FW::RandomEngine generator = m_cfg.randomNumberSvc->spawnGenerator(context);
  std::normal_distribution<double> gauss(0, 1);

  // Read truth particles from input collection
  const auto& simulatedEvent
      = context.eventStore.get<std::vector<Data::SimVertex>>(
          m_cfg.simulatedEventCollection);
  ACTS_DEBUG("Read collection '" << m_cfg.simulatedEventCollection << "' with "
                                 << simulatedEvent.size() << " vertices");

  // Read truth hits from input collection
  const auto& simHits
      = context.eventStore.get<SimHits>(m_cfg.simulatedHitCollection);
  ACTS_DEBUG("Retrieved hit data '" << m_cfg.simulatedHitCollection
                                    << "' from event store.");

  // Prepare the output data
  std::vector<std::pair<size_t, Acts::MultiTrajectory<Data::SimSourceLink>>>
      fittedTracks;

  // Prepare the measurements for KalmanFitter
  ACTS_DEBUG("Prepare the measurements and then tracks");
  std::map<barcode_type, std::vector<Data::SimSourceLink>> sourceLinkMap;
  for (const Data::SimHit& hit : simHits) {
    // get the barcode of the particle associated to the hit
    barcode_type barcode = hit.particle.barcode();
    // get the surface of hit
    const auto& hitSurface = *hit.surface;

    // transform global into local position
    Acts::Vector2D local(0, 0);
    Acts::Vector3D mom(1, 1, 1);
    hitSurface.globalToLocal(context.geoContext, hit.position, mom, local);

    // smear the truth hit with a gaussian and set the covariance
    double                  resX = m_cfg.measurementSigma[0] * Acts::units::_um;
    double                  resY = m_cfg.measurementSigma[1] * Acts::units::_um;
    Acts::ActsSymMatrixD<2> cov2D;
    cov2D << resX * resX, 0., 0., resY * resY;

    double dx = resX * gauss(generator);
    double dy = resY * gauss(generator);

    // move a ,LOC_0, LOC_1 measurement
    // make source link which includes smeared hit

    Acts::BoundMatrix cov{};
    cov.topLeftCorner<2, 2>() = cov2D;

    Acts::BoundVector loc{};
    loc.head<2>() << local[Acts::ParDef::eLOC_0] + dx,
        local[Acts::ParDef::eLOC_1] + dy;

    Data::SimSourceLink sourceLink{&hit, 2, loc, cov};

    // push the truth hits for this particle
    sourceLinkMap[barcode].push_back(sourceLink);
  }

  ACTS_DEBUG("There are " << sourceLinkMap.size() << " tracks for this event ");

  // Get the truth particle
  ACTS_DEBUG("Get truth particle.");
  std::map<barcode_type, Data::SimParticle> particles;
  for (auto& vertex : simulatedEvent) {
    for (auto& particle : vertex.outgoing) {
      particles.insert(std::make_pair(particle.barcode(), particle));
    }
  }

  // Start to perform fit to the prepared tracks
  int itrack = 0;
  for (auto& [barcode, sourceLinks] : sourceLinkMap) {
    if (!barcode) continue;
    itrack++;
    auto particle = particles.find(barcode)->second;

    ACTS_DEBUG("Start processing itrack = "
               << itrack << " with nStates = " << sourceLinks.size()
               << " and truth particle id = " << barcode);

    // get the truth particle info
    Acts::Vector3D pos = particle.position();
    Acts::Vector3D mom = particle.momentum();
    double         q   = particle.q();
    ACTS_DEBUG("truth position = " << pos[0] << " : " << pos[1] << " : "
                                   << pos[2]);
    ACTS_DEBUG("truth momentum = " << mom[0] << " : " << mom[1] << " : "
                                   << mom[2]);

    // smear the truth particle momentum and position
    Acts::Vector3D dir = mom.normalized();

    double fphi   = Acts::VectorHelpers::phi(mom);
    double ftheta = Acts::VectorHelpers::theta(mom);
    double fp     = mom.norm();
    double fqOp   = q / fp;

    // set the smearing error
    double loc0Res  = m_cfg.parameterSigma[0] * Acts::units::_um;
    double loc1Res  = m_cfg.parameterSigma[1] * Acts::units::_um;
    double phiRes   = m_cfg.parameterSigma[2];
    double thetaRes = m_cfg.parameterSigma[3];
    double qOpRes
        = -q / (fp * fp) * m_cfg.parameterSigma[4] * Acts::units::_GeV;

    // prepare the covariance
    Acts::BoundSymMatrix cov;
    cov.setZero();
    cov.diagonal() << pow(loc0Res, 2), pow(loc1Res, 2), pow(phiRes, 2),
        pow(thetaRes, 2), pow(qOpRes, 2), 0;

    // prepare the initial momentum
    double         rPhi   = fphi + phiRes * gauss(generator);
    double         rTheta = ftheta + thetaRes * gauss(generator);
    double         rP     = q / (fqOp + qOpRes * gauss(generator));
    Acts::Vector3D rMom(rP * sin(rTheta) * cos(rPhi),
                        rP * sin(rTheta) * sin(rPhi),
                        rP * cos(rTheta));

    // prepare the initial position
    // the rotation of the starting surface
    Acts::Vector3D T = rMom.normalized();
    Acts::Vector3D U = std::abs(T.dot(Acts::Vector3D::UnitZ())) < 0.99
        ? Acts::Vector3D::UnitZ().cross(T).normalized()
        : Acts::Vector3D::UnitX().cross(T).normalized();
    Acts::Vector3D V = T.cross(U);
    Acts::Vector3D displaced
        = U * loc0Res * gauss(generator) + V * loc1Res * gauss(generator);
    Acts::Vector3D rPos = pos + displaced;

    // then create the start parameters using the prepared momentum and position
    StartParameters rStart(cov, rPos, rMom, q, 0);

    // set the target surface
    const Acts::Surface* rSurface = &rStart.referenceSurface();

    // set the KalmanFitter options
    Acts::KalmanFitterOptions kfOptions(context.geoContext,
                                        context.magFieldContext,
                                        context.calibContext,
                                        rSurface);

    // perform the fit with KalmanFitter
    // auto fittedResult = m_cfg.kFitter.fit(sourceLinks, rStart, kfOptions);
    ACTS_DEBUG("Invoke fitter");
    auto fittedResult = m_cfg.fitFunction(sourceLinks, rStart, kfOptions);
    ACTS_DEBUG("Finish the fitting.");

    // get the fitted parameters
    if (fittedResult.fittedParameters) {
      ACTS_DEBUG("Get the fitted parameters.");
      auto fittedParameters = fittedResult.fittedParameters.get();
      auto fittedPos        = fittedParameters.position();
      auto fittedMom        = fittedParameters.momentum();
      ACTS_DEBUG("Fitted Position at target = " << fittedPos[0] << " : "
                                                << fittedPos[1] << " : "
                                                << fittedPos[2]);
      ACTS_DEBUG("fitted Momentum at target = " << fittedMom[0] << " : "
                                                << fittedMom[1] << " : "
                                                << fittedMom[2]);
    } else
      ACTS_WARNING("No fittedParameter!");

    // get the fitted states
    if (fittedResult.processedStates != 0) {
      ACTS_DEBUG("Get the fitted states.");
      fittedTracks.push_back(
          std::make_pair(fittedResult.trackTip, fittedResult.fittedStates));
    }

    // Make sure the fitting is deterministic
    // auto fittedAgainTrack      = m_cfg.kFitter.fit(measurements, rStart,
    // rSurface);
    // ACTS_DEBUG("Finish fitting again");
    // auto fittedAgainParameters = fittedAgainTrack.fittedParameters.get();

    ACTS_DEBUG("Finish processing the track with particle id = " << barcode);
  }

  // Write the fitted tracks to the EventStore
  context.eventStore.add(m_cfg.trackCollection, std::move(fittedTracks));

  return FW::ProcessCode::SUCCESS;
}
