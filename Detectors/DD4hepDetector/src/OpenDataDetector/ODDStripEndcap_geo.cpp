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
  DetElement endcapDetector(detName, x_det.id());

  // Add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  volConfig.isEndcap             = true;
  Acts::ActsExtension* detvolume = new Acts::ActsExtension(volConfig);
  endcapDetector.addExtension<Acts::IActsExtension>(detvolume);

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string                 endcapShapeName = x_det_dim.nameStr();

  Tube   endcapShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume endcapVolume(detName, endcapShape, lcdd.air());
  endcapVolume.setVisAttributes(lcdd, x_det.visStr());

  // Place Volume
  Volume       motherVolume = lcdd.pickMotherVolume(endcapDetector);
  Position     translation(0., 0., x_det_dim.z());
  PlacedVolume placedEndcap
      = motherVolume.placeVolume(endcapVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedEndcap.addPhysVolID("system", endcapDetector.id());
  endcapDetector.setPlacement(placedEndcap);

  // And return it
  return endcapDetector;
}

DECLARE_DETELEMENT(ODDStripEndcap, create_element)
