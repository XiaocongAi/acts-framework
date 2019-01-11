// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/SurfaceMaterial.hpp"

/// @brief a void material reader
struct VoidMaterialReader
{

  ///@brief operator to be called to retrieve the material
  ///@return a SurfaceMaterialMap
  Acts::SurfaceMaterialMap
  operator()() const
  {
    return Acts::SurfaceMaterialMap();
  }
};
