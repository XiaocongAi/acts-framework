// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

static Ref_t
create_element(Detector& oddd, xml_h xml, SensitiveDetector sens)
{

  xml_det_t x_det   = xml;
  string    detName = x_det.nameStr();

  // Make DetElement
  DetElement cylinderElement(detName, x_det.id());

  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  bool                        isBeamPipe = x_det.hasChild(_U(beampipe));
  volConfig.isPassiveCylinder            = true;
  if (isBeamPipe) {
    volConfig.isBeampipe = true;
  }
  Acts::ActsExtension* detvolume = new Acts::ActsExtension(volConfig);
  cylinderElement.addExtension<Acts::IActsExtension>(detvolume);

  // Make Volume
  xml_comp_t x_det_tubs = x_det.child(_U(tubs));
  string     shapeName  = x_det_tubs.nameStr();
  Tube       tubeShape(
      shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(), x_det_tubs.dz());
  Volume tube_vol(detName, tubeShape, oddd.material(x_det_tubs.materialStr()));
  tube_vol.setVisAttributes(oddd, x_det.visStr());

  // Place it in the mother
  Volume       mother_vol = oddd.pickMotherVolume(cylinderElement);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol);
  placedTube.addPhysVolID(detName, cylinderElement.id());
  cylinderElement.setPlacement(placedTube);

  // And return the element for further parsing
  return cylinderElement;
}

DECLARE_DETELEMENT(ODDCylinder, create_element)
