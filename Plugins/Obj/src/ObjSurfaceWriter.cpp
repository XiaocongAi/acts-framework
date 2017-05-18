#include <iostream>
#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"
#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/GeometryID.hpp"
#include "ACTS/Layers/Layer.hpp"

FWObj::ObjSurfaceWriter::ObjSurfaceWriter(
    const FWObj::ObjSurfaceWriter::Config& cfg)
  : FW::IWriterT<Acts::Surface>()
  , m_cfg(cfg)
{}

FWObj::ObjSurfaceWriter::~ObjSurfaceWriter()
{
}

FW::ProcessCode
FWObj::ObjSurfaceWriter::initialize()
{
  // write out the the file 
  if (!(m_cfg.outputStream)) return FW::ProcessCode::SUCCESS;
  (*(m_cfg.outputStream)) << m_cfg.filePrefix << '\n';
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FWObj::ObjSurfaceWriter::finalize()
{
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FWObj::ObjSurfaceWriter::write(const Acts::Surface& surface)
{

  std::lock_guard<std::mutex> lock(m_write_mutex);
  
  // check
  if (!(m_cfg.outputStream)) return FW::ProcessCode::SUCCESS;
  ACTS_DEBUG(">>Obj: Writer for Surface object called.");

  auto  scalor = m_cfg.outputScalor;
  // let's get the bounds & the transform
  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  auto sTransform = surface.transform();
  
  // dynamic_cast to PlanarBounds
  const Acts::PlanarBounds* planarBounds = 
    dynamic_cast<const Acts::PlanarBounds*>(&surfaceBounds);
  // only continue if the cast worked
  if (planarBounds && m_cfg.outputSensitive){
    ACTS_VERBOSE(">>Obj: Writing out a PlaneSurface");
    // set the precision - just to be sure
    (*(m_cfg.outputStream)) << '\n';
    (*(m_cfg.outputStream)) << std::setprecision(m_cfg.outputPrecision);
    // get the vertices
    auto planarVertices = planarBounds->vertices();
    // loop over the vertices
    std::vector<const Acts::Vector3D> vertices;
    vertices.reserve(planarVertices.size());
    for (auto pv : planarVertices){
      // get the point in 3D
      Acts::Vector3D v3D(sTransform*Acts::Vector3D(pv.x(), pv.y(), 0.));
      vertices.push_back(v3D);   
    }
    // get the thickness and vertical faces
    double thickness = 0.;
    std::vector<unsigned int> vfaces = {};
    if (surface.associatedDetectorElement()){
      // get the thickness form the detector element
      thickness = surface.associatedDetectorElement()->thickness();
      vfaces = { 1, 1, 1, 1 };
    }
    // output to file
    FWObjHelper::writePlanarFace(*(m_cfg.outputStream),
                                 m_vtnCounter,
                                 scalor,
                                 vertices,
                                 thickness,
                                 vfaces);
    (*(m_cfg.outputStream)) << '\n';
  }

 // check if you have layer and check what your have
 //dynamic cast to CylinderBounds work the same 
 const Acts::CylinderBounds* cylinderBounds = 
   dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
 if (cylinderBounds && m_cfg.outputLayerSurface){
   ACTS_VERBOSE(">>Obj: Writing out a CylinderSurface with r = " << cylinderBounds->r());
   // name the object
   auto layerID = surface.geoID().value(Acts::GeometryID::layer_mask);
   (*(m_cfg.outputStream)) << " o Cylinder_" << std::to_string(layerID) << '\n';
   // output to the file
   FWObjHelper::writeTube(*(m_cfg.outputStream),
                          m_vtnCounter,
                          scalor,
                          m_cfg.outputPhiSegemnts,
                          sTransform,
                          cylinderBounds->r(),
                          cylinderBounds->halflengthZ(),
                          m_cfg.outputThickness);
   (*(m_cfg.outputStream)) << '\n';
 }
 
 ////dynamic cast to RadialBounds or disc bounds work the same 
 const Acts::RadialBounds* radialBounds = 
   dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
 if (radialBounds && m_cfg.outputLayerSurface){
   ACTS_VERBOSE(">>Obj: Writing out a DiskSurface at z = " << sTransform.translation().z());
   // name the object
   auto layerID = surface.geoID().value(Acts::GeometryID::layer_mask);
   (*(m_cfg.outputStream)) << "o Disk_" << std::to_string(layerID) << '\n';
   // we use the tube writer in the other direction
   double rMin = radialBounds->rMin();
   double rMax = radialBounds->rMax();
   double thickness = rMax-rMin;
   // output to the file
   FWObjHelper::writeTube(*(m_cfg.outputStream),
                          m_vtnCounter,
                          scalor,
                          m_cfg.outputPhiSegemnts,
                          sTransform,
                          0.5*(rMin+rMax),
                          m_cfg.outputThickness,
                          thickness);
   (*(m_cfg.outputStream)) << '\n';
 }

  // return success 
  return FW::ProcessCode::SUCCESS;
}


