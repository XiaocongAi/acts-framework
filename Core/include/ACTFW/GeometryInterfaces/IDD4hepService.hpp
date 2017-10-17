// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
/// IDD4hepService.hpp
///////////////////////////////////////////////////////////////////

#ifndef GEOMETRYINTERFACES_IDD4HEPSERVICE_H
#define GEOMETRYINTERFACES_IDD4HEPSERVICE_H

#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace dd4hep {
class DetElement;
class Detector;
}

namespace FW {

/// @class IDD4hepService
///
/// Interface class for a service returning the DD4hep geometry.

class IDD4hepService : public IService
{

public:
  /// Virtual destructor
  virtual ~IDD4hepService() = default;
  /// Access to the DD4hep geometry
  /// @return The world DD4hep DetElement
  virtual dd4hep::DetElement
  dd4hepGeometry()
      = 0;
  /// Access to the interface of the DD4hep geometry
  virtual dd4hep::Detector*
  lcdd()
      = 0;
};
}
#endif  // GEOMETRYINTERFACES_IDD4HEPSERVICE_H
