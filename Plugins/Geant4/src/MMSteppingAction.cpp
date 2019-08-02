// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include <stdexcept>
#include "Acts/Utilities/Units.hpp"
#include "G4Material.hh"
#include "G4Step.hh"

FW::Geant4::MMSteppingAction* FW::Geant4::MMSteppingAction::fgInstance
    = nullptr;

FW::Geant4::MMSteppingAction*
FW::Geant4::MMSteppingAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::Geant4::MMSteppingAction::MMSteppingAction()
  : G4UserSteppingAction(), m_steps()
// m_volMgr(MaterialRunAction::Instance()->getGeant4VolumeManager())
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::Geant4::MMSteppingAction::~MMSteppingAction()
{
  fgInstance = nullptr;
}

void
FW::Geant4::MMSteppingAction::UserSteppingAction(const G4Step* step)
{
  // get the material
  G4Material* material = step->GetPreStepPoint()->GetMaterial();

  if (material && material->GetName() != "Vacuum"
      && material->GetName() != "Air") {
    // go through the elements of the material & weigh it with its fraction
    const G4ElementVector* elements  = material->GetElementVector();
    const G4double*        fraction  = material->GetFractionVector();
    size_t                 nElements = material->GetNumberOfElements();
    double                 A         = 0.;
    double                 Z         = 0.;
    double                 X0        = material->GetRadlen();
    double                 L0        = material->GetNuclearInterLength();
    double rho        = material->GetDensity() * CLHEP::mm3 / CLHEP::gram;
    double steplength = step->GetStepLength() / CLHEP::mm;
    if (nElements == 1) {
      A = material->GetA() * CLHEP::mole / CLHEP::gram;
      Z = material->GetZ();
    } else {
      for (size_t i = 0; i < nElements; i++) {
        A += elements->at(i)->GetA() * fraction[i];
        Z += elements->at(i)->GetZ() * fraction[i];
      }
      if (A != 0.) A /= nElements;
      if (A != 0.) Z /= nElements;
    }

    // create the RecordedMaterialProperties
    const auto& prePos  = step->GetPreStepPoint()->GetPosition();
    const auto& postPos = step->GetPostStepPoint()->GetPosition();

    Acts::Vector3D presPos(prePos.x(), prePos.y(), prePos.z());
    Acts::Vector3D possPos(postPos.x(), postPos.y(), postPos.z());

    Acts::MaterialInteraction mInteraction;
    mInteraction.position       = 0.5 * (presPos + possPos);
    mInteraction.direction      = (possPos - presPos).normalized();
    mInteraction.pathCorrection = step->GetStepLength();
    mInteraction.materialProperties
        = Acts::MaterialProperties(X0, L0, A, Z, rho, steplength);
    m_steps.push_back(mInteraction);
  }
}

void
FW::Geant4::MMSteppingAction::Reset()
{
  m_steps.clear();
}
