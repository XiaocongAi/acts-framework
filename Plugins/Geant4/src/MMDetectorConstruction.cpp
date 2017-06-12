#include "ACTFW/Plugins/Geant4/MMDetectorConstruction.hpp"
// VGM demo
//#include "Geant4GM/volumes/Factory.h"
//#include "RootGM/volumes/Factory.h"
#include "TGeoManager.h"
// end VGM demo
#include "G4GDMLParser.hh"

FWG4::MMDetectorConstruction::MMDetectorConstruction():
G4VUserDetectorConstruction(),
m_tgeoNode(nullptr),
m_gdmlFile(nullptr)
{
}


FWG4::MMDetectorConstruction::~MMDetectorConstruction()
{
}

G4VPhysicalVolume* FWG4::MMDetectorConstruction::Construct()
{
    if (m_tgeoNode) {
        // Import geometry from Root to VGM
   /*     RootGM::Factory rtFactory;
        rtFactory.SetDebug(1);
        rtFactory.Import(m_tgeoNode);
    
        // Export VGM geometry to Geant4
        Geant4GM::Factory g4Factory;
        g4Factory.SetDebug(1);
        rtFactory.Export(&g4Factory);
        G4VPhysicalVolume* world = g4Factory.World();
        return world;*/
        return nullptr;
    }
    else if(m_gdmlFile) {
        G4GDMLParser parser;
        parser.Read(*m_gdmlFile);
        return parser.GetWorldVolume();
    }
    else return nullptr; // and error Message
}

void FWG4::MMDetectorConstruction::setTGeoGeometry(TGeoNode* tgeoNode)
{
    m_tgeoNode = tgeoNode;
}

void FWG4::MMDetectorConstruction::setGdmlInput(std::string gdmlFile)
{
    m_gdmlFile = new std::string(gdmlFile);
}