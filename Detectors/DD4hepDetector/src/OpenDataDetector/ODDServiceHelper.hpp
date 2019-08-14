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

template <typename volume_t>
void
buildBarrelRouting(Detector&                  oddd,
                   volume_t&                  barrelVolume,
                   const xml_comp_t&          x_routing,
                   const std::vector<double>& layerR)
{

  // Grab the cables & route them outwards
  unsigned int nphi = x_routing.nphi();

  double phiStep = 2 * M_PI / nphi;
  double phi0    = x_routing.phi0();
  double rmin    = x_routing.rmin();
  double rmax    = x_routing.rmax();
  double n       = x_routing.number();

  for (int side = -1; side < 2; side += 2) {
    // Loop over the layer routings
    for (unsigned int ib = 1; ib < layerR.size(); ++ib) {
      for (unsigned int iphi = 0; iphi < nphi; ++iphi) {
        // Calculate the phi
        double phi = phi0 + iphi * phiStep;

        // The layer position
        double gap     = x_routing.gap();
        double clength = layerR[ib] - layerR[ib - 1] - 2. * gap;
        double rpos    = 0.5 * (layerR[ib] + layerR[ib - 1]);
        double xpos    = rpos * cos(phi);
        double ypos    = rpos * sin(phi);
        double zpos    = side * x_routing.z_offset();

        Assembly cableboxAssembly("CableBox");
        if (x_routing.hasChild(_U(box))) {
          // The box plate for the cables
          xml_comp_t x_box = x_routing.child(_U(box));
          Box        box(x_box.dz(), n * ib * rmax, 0.5 * clength);
          Volume     boxVolume(
              "CableBand", box, oddd.material(x_routing.materialStr()));
          boxVolume.setVisAttributes(oddd, x_box.visStr());

          PlacedVolume pacedBox = cableboxAssembly.placeVolume(
              boxVolume, Position(side * (rmax + x_box.dz()), 0., 0.));
        }

        Tube   cable(rmin, rmax, 0.5 * clength);
        Volume cableVolume(
            "Cable", cable, oddd.material(x_routing.materialStr()));
        cableVolume.setVisAttributes(oddd, x_routing.visStr());

        for (unsigned int icable = 0; icable < n * ib; ++icable) {

          // Place the pipe in the stave
          PlacedVolume placedCable = cableboxAssembly.placeVolume(
              cableVolume, Position(0., (-n * ib + 1 + 2 * icable) * rmax, 0.));
        }
        // Place the pipe in the stave
        PlacedVolume placedCableBox = barrelVolume.placeVolume(
            cableboxAssembly,
            Transform3D(RotationZ(phi) * RotationY(0.5 * M_PI),
                        Position(xpos, ypos, zpos)));
      }
    }
  }
}

template <typename volume_t>
void
buildEndcapRouting(Detector&                  oddd,
                   volume_t&                  barrelVolume,
                   const xml_comp_t&          x_routing,
                   const std::vector<double>& endcapZ)
{

  // Grab the cables & route them outwards
  unsigned int nphi = x_routing.nphi();

  double phiStep = 2 * M_PI / nphi;
  double phi0    = x_routing.phi0();
  double rmin    = x_routing.rmin();
  double rmax    = x_routing.rmax();
  double r       = x_routing.r();
  double n       = x_routing.number();

  // Loop over the layer routings
  for (unsigned int iec = 1; iec < endcapZ.size(); ++iec) {
    for (unsigned int iphi = 0; iphi < nphi; ++iphi) {
      // Calculate the phi
      double phi = phi0 + iphi * phiStep;

      // The layer position
      double gap     = x_routing.gap();
      double clength = std::abs(endcapZ[iec] - endcapZ[iec - 1]) - 2. * gap;
      double xpos    = r * cos(phi);
      double ypos    = r * sin(phi);
      double zpos    = 0.5 * (endcapZ[iec] + endcapZ[iec - 1]);

      Assembly cableboxAssembly("CableBox");
      if (x_routing.hasChild(_U(box))) {
        // The box plate for the cables
        xml_comp_t x_box = x_routing.child(_U(box));
        Box        box(x_box.dz(), n * iec * rmax, 0.5 * clength);
        Volume     boxVolume(
            "CableBand", box, oddd.material(x_routing.materialStr()));
        boxVolume.setVisAttributes(oddd, x_box.visStr());

        PlacedVolume pacedBox = cableboxAssembly.placeVolume(
            boxVolume, Position(rmax + x_box.dz(), 0., 0.));
      }

      Tube   cable(rmin, rmax, 0.5 * clength);
      Volume cableVolume(
          "Cable", cable, oddd.material(x_routing.materialStr()));
      cableVolume.setVisAttributes(oddd, x_routing.visStr());

      for (unsigned int icable = 0; icable < n * iec; ++icable) {

        // Place the pipe in the stave
        PlacedVolume placedCable = cableboxAssembly.placeVolume(
            cableVolume, Position(0., (-n * iec + 1 + 2 * icable) * rmax, 0.));
      }
      // Place the pipe in the stave
      PlacedVolume placedCableBox = barrelVolume.placeVolume(
          cableboxAssembly,
          Transform3D(RotationZ(+phi), Position(xpos, ypos, zpos)));
    }
  }
}
