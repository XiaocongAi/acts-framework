// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Json/JsonMaterialWriter.hpp"
#include <fstream>
#include <ios>
#include <iostream>
#include <stdexcept>
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Utilities/GeometryID.hpp"

FW::Json::JsonMaterialWriter::JsonMaterialWriter(
    const FW::Json::JsonGeometryConverter::Config& cfg,
    const std::string&                             fileName)
  : m_cfg(cfg), m_fileName(fileName)
{
  // Validate the configuration
  if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing service name");
  }
}

FW::Json::JsonMaterialWriter::~JsonMaterialWriter()
{
}

void
FW::Json::JsonMaterialWriter::write(
    const Acts::DetectorMaterialMaps& detMaterial)
{
  // Evoke the converter
  FW::Json::JsonGeometryConverter jmConverter(m_cfg);
  auto jout = jmConverter.materialMapsToJson(detMaterial);
  // And write the file
  std::ofstream ofj(m_fileName);
  ofj << std::setw(4) << jout << std::endl;
}

void
FW::Json::JsonMaterialWriter::write(const Acts::TrackingGeometry& tGeometry)
{
  // Evoke the converter
  FW::Json::JsonGeometryConverter jmConverter(m_cfg);
  auto jout = jmConverter.trackingGeometryToJson(tGeometry);
  // And write the file
  std::ofstream ofj(m_fileName);
  ofj << std::setw(4) << jout << std::endl;
}