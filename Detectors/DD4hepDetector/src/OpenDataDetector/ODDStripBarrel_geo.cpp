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
create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens)
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

  Tube   barrelShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume barrelVolume(detName, barrelShape, lcdd.air());
  barrelVolume.setVisAttributes(lcdd, x_det.visStr());

  // Place Volume
  Volume       motherVolume = lcdd.pickMotherVolume(barrelDetector);
  Position     translation(0., 0., x_det_dim.z());
  PlacedVolume placedBarrel
      = motherVolume.placeVolume(barrelVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedBarrel.addPhysVolID("system", barrelDetector.id());
  barrelDetector.setPlacement(placedBarrel);

  // And return it
  return barrelDetector;
}

DECLARE_DETELEMENT(ODDStripBarrel, create_element)
