// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"
#include "ODDModuleHelper.hpp"

using namespace std;
using namespace dd4hep;
static void
completeStaveStructure(Detector&   oddd,
                       xml_comp_t& x_stave,
                       Assembly&   staveAssembly,
                       double      staveHlength,
                       double      ylength)
{

  unsigned int nModules = x_stave.nmodules();

  // Place carbon foam structure
  if (x_stave.hasChild(_U(subtraction)) and x_stave.hasChild(_U(tube))) {
    xml_comp_t x_sub  = x_stave.child(_U(subtraction));
    xml_comp_t x_trd  = x_sub.child(_U(trd));
    xml_comp_t x_tubs = x_sub.child(_U(tubs));
    xml_comp_t x_pipe = x_stave.child(_U(tube));

    // Create the two shapes first
    Trapezoid foamShape(x_trd.x1(),
                        x_trd.x2(),
                        staveHlength + x_trd.dz(),
                        staveHlength + x_trd.dz(),
                        x_trd.thickness());

    Tube foamCutout(x_tubs.rmin(), x_tubs.rmax(), staveHlength + x_tubs.dz());

    // Create the substraction
    Volume foamVolume("CarbonFoam",
                      SubtractionSolid(foamShape,
                                       foamCutout,
                                       Transform3D(RotationX(0.5 * M_PI))),
                      oddd.material(x_sub.materialStr()));
    foamVolume.setVisAttributes(oddd, x_sub.visStr());
    // Place the volume in the stave
    PlacedVolume placedFoam = staveAssembly.placeVolume(
        foamVolume,
        Position(x_sub.x_offset(), x_sub.y_offset(), x_sub.z_offset()));

    Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(), staveHlength + x_pipe.dz());
    // Create the substraction
    Volume pipeVolume(
        "CoolingPipe", coolingPipe, oddd.material(x_pipe.materialStr()));
    pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

    // Place the pipe in the stave
    PlacedVolume placedPipe = staveAssembly.placeVolume(
        pipeVolume,
        Transform3D(
            RotationX(0.5 * M_PI),
            Position(x_pipe.x_offset(), x_pipe.y_offset(), x_pipe.z_offset())));

    xml_comp_t x_cable = x_stave.child(_U(eltube));
    // Place the support cables for the modules
    for (unsigned int modCable = 0; modCable < 0.5 * nModules; ++modCable) {

      double cableLength = staveHlength - modCable * ylength;

      for (int side = -1; side < 2; side += 2) {
        Tube cable(x_cable.rmin(), x_cable.rmax(), 0.5 * cableLength);
        // Create the cable volume
        Volume cableVolume(
            "Cable", cable, oddd.material(x_cable.materialStr()));
        cableVolume.setVisAttributes(oddd, x_cable.visStr());

        // Place the pipe in the stave
        PlacedVolume placedCable = staveAssembly.placeVolume(
            cableVolume,
            Transform3D(
                RotationX(0.5 * M_PI),
                Position(
                    x_cable.x_offset() + 2.05 * modCable * x_cable.rmax(),
                    side * (staveHlength - 0.5 * cableLength + x_cable.dz()),
                    x_cable.z_offset())));
      }
    }
  }
}

static Ref_t
create_element(Detector& oddd, xml_h xml, SensitiveDetector sens)
{

  xml_det_t x_det   = xml;
  string    detName = x_det.nameStr();

  // Make DetElement
  DetElement barrelDetector(detName, x_det.id());

  // Add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  volConfig.isBarrel             = true;
  Acts::ActsExtension* detvolume = new Acts::ActsExtension(volConfig);
  barrelDetector.addExtension<Acts::IActsExtension>(detvolume);

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string                 barrelShapeName = x_det_dim.nameStr();

  // The Shape and Volume
  Tube   barrelShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume barrelVolume(detName, barrelShape, oddd.air());
  barrelVolume.setVisAttributes(oddd, x_det.visStr());

  // Create the stave volume and DetElement tree
  xml_comp_t x_stave = x_det.child(_U(stave));
  Assembly   staveAssembly("stave");
  // Visualization
  staveAssembly.setVisAttributes(oddd, x_stave.visStr());
  // DetElement tree
  DetElement staveElementTemplate("StaveElementTemplate", 0);

  // Create the module components
  xml_comp_t x_module = x_det.child(_U(module));
  double     ylength  = 0.;
  auto       module = assembleRectangularModule(oddd, sens, x_module, ylength);

  // Place the modules into the stave
  double       gap          = x_stave.gap();
  unsigned int nModules     = x_stave.nmodules();
  double       ystep        = ylength + gap;
  double       ymin         = (nModules * 0.5 - 0.5) * ylength;
  double       staveHlength = ymin + 0.5 * ylength;

  // Loop over the modules and place them in the stave
  for (unsigned int moduleNum = 0; moduleNum < nModules; ++moduleNum) {
    // Place them along local y
    PlacedVolume placedModule = staveAssembly.placeVolume(
        module.first, Position(0., -ymin + moduleNum * ystep, 0.));
    placedModule.addPhysVolID("module", moduleNum);

    string moduleName = _toString((int)moduleNum, "module%d");
    // Clone the detector element
    auto moduleElement = module.second.clone(moduleName, moduleNum);
    moduleElement.setPlacement(placedModule);
    // Assign it as child to the stave template
    staveElementTemplate.add(moduleElement);
  }

  // Complete the staveStructure
  completeStaveStructure(oddd, x_stave, staveAssembly, staveHlength, ylength);

  // Loop over the layers to build staves
  size_t layerNum = 0;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layerNum) {
    xml_comp_t x_layer = lay;

    // The Layer envelope volume
    string layerName = detName + std::to_string(layerNum);
    Volume layerVolume(
        layerName,
        Tube(x_layer.rmin(), x_layer.rmax(), staveHlength + x_layer.outer_z()),
        oddd.air());
    // Visualization
    layerVolume.setVisAttributes(oddd, x_layer.visStr());

    // The DetElement tree, keep it flat
    DetElement layerElement(barrelDetector, layerName, layerNum);

    // Place the staves in the layer
    unsigned int nStaves = x_layer.nphi();
    double       phiStep = 2. * M_PI / nStaves;
    double       phiTilt = x_layer.phi_tilt();
    double       phi0    = x_layer.phi0();
    double       r       = x_layer.r();

    // Loop over the staves and place them
    for (unsigned int staveNum = 0; staveNum < nStaves; ++staveNum) {
      string staveName = _toString((int)staveNum, "stave%d");
      // position of the stave
      double phi = phi0 + staveNum * phiStep;
      double x   = r * cos(phi);
      double y   = r * sin(phi);
      // Now place the stave
      PlacedVolume placedStave = layerVolume.placeVolume(
          staveAssembly,
          Transform3D(RotationY(0.5 * M_PI) * RotationZ(0.5 * M_PI)
                          * RotationY(phi + phiTilt),
                      Position(x, y, 0.)));
      placedStave.addPhysVolID("stave", staveNum);

      // Clone the stave element from the template
      DetElement staveElement = staveElementTemplate.clone(staveName, staveNum);
      staveElement.setPlacement(placedStave);
      // Add to the layer element
      layerElement.add(staveElement);
    }

    // Place the support cylinder
    if (x_layer.hasChild(_U(support))) {
      xml_comp_t x_support = x_layer.child(_U(support));
      // Create the volume of the support structure
      Volume supportVolume(
          "SupportCylinder",
          Tube(x_support.rmin(), x_support.rmax(), x_support.dz()),
          oddd.material(x_support.materialStr()));
      supportVolume.setVisAttributes(oddd, x_support.visStr());
      // Place the support structure
      PlacedVolume placedSupport = layerVolume.placeVolume(supportVolume);
    }

    // Cleanup the templates
    // dd4hep::detail::destroyHandle(staveElementTemplate);
    // dd4hep::detail::destroyHandle(module.second);

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension::Config layerConfig;
    layerConfig.isLayer                 = true;
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension(layerConfig);
    layerElement.addExtension<Acts::IActsExtension>(layerExtension);

    PlacedVolume placedLayer = barrelVolume.placeVolume(layerVolume);
    placedLayer.addPhysVolID("layer", layerNum);

    // Assign layer DetElement to layer volume
    layerElement.setPlacement(placedLayer);
  }

  // Place Volume
  Volume       motherVolume = oddd.pickMotherVolume(barrelDetector);
  Position     translation(0., 0., x_det_dim.z());
  PlacedVolume placedBarrel
      = motherVolume.placeVolume(barrelVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedBarrel.addPhysVolID("system", barrelDetector.id());
  barrelDetector.setPlacement(placedBarrel);

  // And return the detector element
  return barrelDetector;
}

DECLARE_DETELEMENT(ODDPixelBarrel, create_element)
