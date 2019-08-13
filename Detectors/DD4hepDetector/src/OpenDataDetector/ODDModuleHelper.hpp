// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

namespace Acts {
class DigitizationModule;
}

static std::pair<Assembly, DetElement>
assembleRectangularModule(Detector&          oddd,
                          SensitiveDetector& sens,
                          const xml_comp_t&  x_module,
                          double&            ylength)
{

  // The Module envelope volume
  Assembly moduleAssembly("module");
  // Visualization
  moduleAssembly.setVisAttributes(oddd, x_module.visStr());

  // The module detector element
  DetElement moduleElement("ModuleElementTemplate", 0);

  // Place the components inside the module
  unsigned int compNum   = 0;
  unsigned int sensorNum = 0;

  for (xml_coll_t comp(x_module, _U(module_component)); comp;
       ++comp, ++compNum) {
    xml_comp_t x_comp = comp;

    // Component volume
    string componentName = _toString((int)compNum, "component%d");
    Box    boxShape(0.5 * x_comp.dx(), 0.5 * x_comp.dy(), 0.5 * x_comp.dz());
    // Standard component volume without cutout
    Volume componentVolume(
        componentName, boxShape, oddd.material(x_comp.materialStr()));

    // overwrite if you have a subtraction
    if (x_comp.hasChild(_U(subtraction))) {
      std::cout << "Subtraction detected !" << std::endl;
      xml_comp_t x_sub = x_comp.child(_U(subtraction));
      Tube       tubeCutout(x_sub.rmin(), x_sub.rmax(), x_comp.dy());

      // Create the substraction
      componentVolume
          = Volume(componentName,
                   SubtractionSolid(boxShape,
                                    tubeCutout,
                                    Transform3D(RotationX(0.5 * M_PI),
                                                Position(x_sub.x_offset(),
                                                         x_sub.y_offset(),
                                                         x_sub.z_offset()))),
                   oddd.material(x_comp.materialStr()));

      // place a fitting pipe if available
      if (x_comp.hasChild(_U(tube))) {
        xml_comp_t x_pipe = x_comp.child(_U(tube));
        Tube       coolingPipe(x_pipe.rmin(), x_pipe.rmax(), 0.5 * x_comp.dy());
        // Create the substraction
        Volume pipeVolume(
            "CoolingPipe", coolingPipe, oddd.material(x_pipe.materialStr()));
        pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

        PlacedVolume pacedPipe = componentVolume.placeVolume(
            pipeVolume,
            Transform3D(RotationX(0.5 * M_PI),
                        Position(x_pipe.x_offset(),
                                 x_pipe.y_offset(),
                                 x_pipe.z_offset())));
      }
    }
    componentVolume.setVisAttributes(oddd, x_comp.visStr());

    // Calculate the module dimension
    double cylength
        = 2. * abs(std::copysign(0.5 * x_comp.dy(), x_comp.y_offset())
                   + x_comp.y_offset());
    ylength = cylength > ylength ? cylength : ylength;

    // Visualization
    componentVolume.setVisAttributes(oddd, x_comp.visStr());
    // Place Module Box Volumes in layer
    double       stereoAlpha     = x_comp.alpha();
    PlacedVolume placedComponent = moduleAssembly.placeVolume(
        componentVolume,
        Transform3D(
            RotationZ(stereoAlpha),
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset())));

    // Deal with the sensitive sensor
    if (x_comp.isSensitive()) {

      componentVolume.setSensitiveDetector(sens);
      placedComponent.addPhysVolID("sensor", sensorNum++);

      // Create the sensor element and place it
      string     sensorName = _toString((int)sensorNum, "sensor%d");
      DetElement sensorElement(moduleElement, sensorName, sensorNum);
      sensorElement.setPlacement(placedComponent);

      // Add the sensor extension
      // std::shared_ptr<const Acts::DigitizationModule> digiModule = nullptr;
      // Acts::ActsExtension* sensorExtension = new
      // Acts::ActsExtension(digiModule);
      // sensorElement.addExtension<Acts::IActsExtension>(sensorExtension);
    }
  }

  // return the module assembly
  return std::pair<Assembly, DetElement>(moduleAssembly, moduleElement);
}
