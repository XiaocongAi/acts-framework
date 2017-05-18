#include "ACTFW/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

FWG4::MMPrimaryGeneratorAction* FWG4::MMPrimaryGeneratorAction::fgInstance = 0;

FWG4::MMPrimaryGeneratorAction::MMPrimaryGeneratorAction(
    const G4String& particleName,
    G4double        energy,
    G4int           randomSeed1,
    G4int           randomSeed2)
  : G4VUserPrimaryGeneratorAction()
  , fParticleGun(0)
{
  // configure the run
  fgInstance         = this;
  G4int nofParticles = 1;
  fParticleGun       = new G4ParticleGun(nofParticles);
  // default particle kinematic
  G4ParticleTable*      particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(energy);
  G4UnitDefinition::PrintUnitsTable();
  // set the random seeds 
  CLHEP::HepRandom::getTheEngine()->setSeed(randomSeed1,randomSeed2);
}

FWG4::MMPrimaryGeneratorAction::~MMPrimaryGeneratorAction()
{
  fgInstance = 0;
  delete fParticleGun;
}

FWG4::MMPrimaryGeneratorAction*
FWG4::MMPrimaryGeneratorAction::Instance()
{
  // Static acces function via G4RunManager

  return fgInstance;
}

void
FWG4::MMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  G4double phi   = -M_PI + G4UniformRand() * 2. * M_PI;
  G4double theta = G4UniformRand() * M_PI;
  // build a direction
  m_direction
    = G4ThreeVector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  m_position = G4ThreeVector(0., 0., 0.); /// @todo make configurable G4RandGauss::shoot(0., 150.));
  // set to the particle gun and 
  fParticleGun->SetParticleMomentumDirection(m_direction);
  fParticleGun->SetParticlePosition(m_position);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
