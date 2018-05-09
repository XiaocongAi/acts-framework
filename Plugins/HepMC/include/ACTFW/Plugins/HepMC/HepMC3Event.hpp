// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "HepMC/FourVector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepPID/ParticleIDMethods.hh"

#include "ACTFW/Utilities/SimEvent.hpp"
#include "ACTFW/Utilities/SimParticle.hpp"
#include "ACTFW/Utilities/SimVertex.hpp"

namespace FW {

/// @struct SimulatedEvent<HepMC::GenEvent>
///
/// This struct is an explicit implementation of FW::SimulatedEvent for the
/// translation of HepMC::GenEvent objects into Acts.
///
template <>
struct SimulatedEvent<HepMC::GenEvent>
{
public:
  ///
  /// Setter
  ///

  /// @brief Sets new units for momentums
  /// @note The allowed units are MeV and Gev
  /// @param event event in HepMC data type
  /// @param momentumUnit new unit of momentum
  static void
  momentumUnit(std::shared_ptr<HepMC::GenEvent> event,
               const double                     momentumUnit);

  /// @brief Sets new units for lengths
  /// @note The allowed units are mm and cm
  /// @param event event in HepMC data type
  /// @param lengthUnit new unit of length
  static void
  lengthUnit(std::shared_ptr<HepMC::GenEvent> event, const double lengthUnit);

  /// @brief Shifts the positioning of an event in space and time
  /// @param event event in HepMC data type
  /// @param deltaPos relative spatial shift that will be applied
  /// @param deltaTime relative time shift that will be applied
  static void
  shiftPositionBy(std::shared_ptr<HepMC::GenEvent> event,
                  const Acts::Vector3D&            deltaPos,
                  const double                     deltaTime);

  /// @brief Shifts the positioning of an event to a paint in space and time
  /// @param event event in HepMC data type
  /// @param pos new position of the event
  /// @param time new time of the event
  static void
  shiftPositionTo(std::shared_ptr<HepMC::GenEvent> event,
                  const Acts::Vector3D&            pos,
                  const double                     time);

  /// @brief Shifts the positioning of an event to a paint in space
  /// @param event event in HepMC data type
  /// @param pos new position of the event
  static void
  shiftPositionTo(std::shared_ptr<HepMC::GenEvent> event,
                  const Acts::Vector3D&            pos);

  /// @brief Shifts the positioning of an event to a paint in time
  /// @param event event in HepMC data type
  /// @param time new time of the event
  static void
  shiftPositionTo(std::shared_ptr<HepMC::GenEvent> event, const double time);

  ///
  /// Adder
  ///

  /// @brief Adds a new particle
  /// @param event event in HepMC data type
  /// @param particle new particle that will be added
  static void
  addParticle(std::shared_ptr<HepMC::GenEvent>          event,
              std::shared_ptr<Acts::ParticleProperties> particle);

  /// @brief Adds a new vertex
  /// @param event event in HepMC data type
  /// @param vertex new vertex that will be added
  /// @note The statuses are not represented in Acts and therefore set to 0
  static void
  addVertex(std::shared_ptr<HepMC::GenEvent>           event,
            const std::shared_ptr<Acts::ProcessVertex> vertex);
  ///
  /// Remover
  ///

  /// @brief Removes a particle from the record
  /// @param event event in HepMC data type
  /// @param particle particle that will be removed
  static void
  removeParticle(std::shared_ptr<HepMC::GenEvent>                 event,
                 const std::shared_ptr<Acts::ParticleProperties>& particle);

  /// @brief Removes a vertex from the record
  /// @note The identification of the vertex is potentially unstable (c.f.
  /// HepMC3Event::compareVertices())
  /// @param event event in HepMC data type
  /// @param vertex vertex that will be removed
  static void
  removeVertex(std::shared_ptr<HepMC::GenEvent>            event,
               const std::shared_ptr<Acts::ProcessVertex>& vertex);

  ///
  /// Getter
  ///

  /// @brief Getter of the unit of momentum used
  /// @param event event in HepMC data type
  /// @return unit of momentum
  static double
  momentumUnit(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Getter of the unit of length used
  /// @param event event in HepMC data type
  /// @return unit of length
  static double
  lengthUnit(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Getter of the position of the event
  /// @param event event in HepMC data type
  /// @return vector to the location of the event
  static Acts::Vector3D
  eventPos(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Getter of the time of the event
  /// @param event event in HepMC data type
  /// @return time of the event
  static double
  eventTime(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Get list of particles
  /// @param event event in HepMC data type
  /// @return List of particles
  static std::vector<std::unique_ptr<Acts::ParticleProperties>>
  particles(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Get list of vertices
  /// @param event event in HepMC data type
  /// @return List of vertices
  static std::vector<std::unique_ptr<Acts::ProcessVertex>>
  vertices(const std::shared_ptr<HepMC::GenEvent> event);

  /// @brief Get beam particles
  /// @param event event in HepMC data type
  /// @return List of beam particles
  static std::vector<std::unique_ptr<Acts::ParticleProperties>>
  beams(const std::shared_ptr<HepMC::GenEvent> event);

private:
  /// @brief Converts an Acts::ParticleProperties into HepMC::GenParticle
  /// @note The conversion ignores HepMC status codes
  /// @param actsParticle Acts particle that will be converted
  /// @return converted particle
  static HepMC::GenParticlePtr
  actsParticleToGen(std::shared_ptr<Acts::ParticleProperties> actsParticle);

  /// @brief Converts an Acts vertex to a HepMC::GenVertexPtr
  /// @note The conversion ignores HepMC status codes
  /// @param actsVertex Acts vertex that will be converted
  /// @return Converted Acts vertex to HepMC::GenVertexPtr
  static HepMC::GenVertexPtr
  createGenVertex(const std::shared_ptr<Acts::ProcessVertex>& actsVertex);

  /// @brief Compares an Acts vertex with a HepMC::GenVertex
  /// @note An Acts vertex does not store a barcode. Therefore the content of
  /// both vertices is compared. The position, time and number of incoming and
  /// outgoing particles will be compared. Since a second vertex could exist in
  /// the record with identical informations (although unlikely), this
  /// comparison could lead to false positive results. On the other hand, a
  /// numerical deviation of the parameters could lead to a false negative.
  /// @param actsVertex Acts vertex
  /// @param genVertex HepMC::GenVertex
  /// @return boolean result if both vertices are identical
  static bool
  compareVertices(const std::shared_ptr<Acts::ProcessVertex>& actsVertex,
                  const HepMC::GenVertexPtr&                  genVertex);
};
}  // FW
