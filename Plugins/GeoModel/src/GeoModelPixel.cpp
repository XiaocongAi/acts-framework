// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/GeoModel/GeoModelPixel.hpp"

#include <iostream>
#include "Acts/Utilities/Definitions.hpp"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoTube.h"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "GeoModelKernel/GeoShapeShift.h"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"

// Units
#include "GeoModelKernel/Units.h"
#define SYSTEM_OF_UNITS GeoModelKernelUnits

namespace {
double conv = Acts::units::_mm / SYSTEM_OF_UNITS::mm;
}

std::unique_ptr<Acts::SurfaceArray>
FW::GeoModelPixel::surfaceArray(GeoVPhysVol const* lay) const
{
	// TODO: debug exit!
	return nullptr;
	
	std::cout << "\tThis is my lay: " << lay->getLogVol()->getName() << "\t" << lay->getLogVol()->getShape()->type() << std::endl;
	auto tube = dynamic_cast<GeoTube const*>(lay->getLogVol()->getShape());
	std::cout << "\t" << tube->getRMin() << "\t" << tube->getRMax() << "\t" << tube->getZHalfLength() << std::endl;
	
	// Walk through the content of the layer
	for(unsigned int i = 0; i < lay->getNChildVols(); i++)
	{
		GeoVPhysVol const* child = &(*(lay->getChildVol(i)));
		GeoShape const* shape = child->getLogVol()->getShape();
				
		// Only consider the "Ladder"
		if(child->getLogVol()->getName() == "Ladder")
		{
			std::cout << "ccc: " << i << "\n" << lay->getXToChildVol(i).matrix() << std::endl;
			
			if(shape->type() == "Box")
			{
				auto box = dynamic_cast<GeoBox const*>(shape);
				std::cout << "Box: " << box->getXHalfLength() << "\t" << box->getYHalfLength() << "\t" << box->getZHalfLength() << std::endl;
			}
			
			// TODO: The following lines look for the ModuleBrl
			// Search for the modules in the ladder
			for(unsigned int j = 0; j < child->getNChildVols(); j++)
			{
				if(child->getChildVol(j)->getLogVol()->getName() != "ModuleBrl")
				{
					continue;
				}
				GeoVPhysVol const* moduleBrl = &(*(child->getChildVol(j)));
				
				/**
				//~ std::cout << "bc: " << j << "\t" << moduleBrl->getLogVol()->getName() << "\t" << moduleBrl->getLogVol()->getShape()->type() << "\t" << moduleBrl->getNChildVols() << std::endl;
				//~ std::cout << "trafo: " << child->getXToChildVol(j).matrix() << std::endl; // TODO: conversion of units

				if(moduleBrl->getLogVol()->getShape()->type() == "Box") // TODO: Shift also exists
				{
					auto box = dynamic_cast<GeoBox const*>(moduleBrl->getLogVol()->getShape());
					
					auto rBounds = std::make_shared<Acts::RectangleBounds>(box->getYHalfLength(), box->getZHalfLength()); // TODO: conversion of units
					
					//~ Acts::PlaneSurface plane(std::make_shared<Acts::Transform3D>(noseToTheWind->getX()), rBounds);
				}
				else
				if(moduleBrl->getLogVol()->getShape()->type() == "Shift")
				{
					auto shift = dynamic_cast<GeoShapeShift const*>(moduleBrl->getLogVol()->getShape());	
					auto box = dynamic_cast<GeoBox const*>(shift->getOp());
					
					auto rBounds = std::make_shared<Acts::RectangleBounds>(box->getYHalfLength(), box->getZHalfLength()); // TODO: conversion of units
					
					//~ std::cout << "shifttrafo: " << shift->getOp()->type() << std::endl << shift->getX().matrix() << std::endl;
				}
				**/
			}
		}
	}
	return nullptr;
}

std::vector<std::shared_ptr<const Acts::Layer>>
FW::GeoModelPixel::buildLayers(GeoVPhysVol const* vol) const
{
	// The resulting layers
	std::vector<std::shared_ptr<const Acts::Layer>> layers;
	
	// Search through all elements of the volume
	unsigned int ncChildren = vol->getNChildVols();
	for(unsigned int j = 0; j < ncChildren; j++)
	{
		auto layer = dynamic_cast<GeoVPhysVol const*>(&(*(vol->getChildVol(j))));
		// Test if the layer is one that we search for
		if(m_layerKeys.find(layer->getLogVol()->getName()) != m_layerKeys.end())
		{
			// Build the transformation
			std::shared_ptr<const Acts::Transform3D> transformation = nullptr;
			if(dynamic_cast<GeoFullPhysVol const*>(layer))
			{
				transformation = std::make_shared<Acts::Transform3D>(dynamic_cast<GeoVFullPhysVol const*>(layer)->getDefAbsoluteTransform());
			}
			
			// Build the bounds & thickness
			std::shared_ptr<const Acts::CylinderBounds> bounds = nullptr;
			double thickness = 0.;
			if(dynamic_cast<GeoTube const*>(layer->getLogVol()->getShape()))
			{
				GeoTube const* tube = dynamic_cast<GeoTube const*>(layer->getLogVol()->getShape());
				bounds = std::make_shared<const Acts::CylinderBounds>(tube->getRMax() * conv, tube->getZHalfLength() * conv);
				thickness = (tube->getRMax() - tube->getRMin()) * conv;
			}
		
			// Build the surface array
			std::unique_ptr<Acts::SurfaceArray> surArray = 	surfaceArray(&(*(vol->getChildVol(j))));
			
			// Create and store the layer
			layers.push_back(Acts::CylinderLayer::create(transformation, bounds, std::move(surArray), thickness, nullptr, Acts::LayerType::active));
		}
	}
	return layers;
}

std::unique_ptr<const Acts::LayerArray>
FW::GeoModelPixel::buildLayerArray(GeoVPhysVol const* vol) const
{
	// Build the layers
	std::vector<std::shared_ptr<const Acts::Layer>> layers = buildLayers(vol);
	
	// Get the minimum/maximum of the volume
	double min, max;		  
	if(dynamic_cast<GeoTube const*>(vol->getLogVol()->getShape()))
	{
		GeoTube const* tube = dynamic_cast<GeoTube const*>(vol->getLogVol()->getShape());
		min = tube->getRMin() * conv;
		max = tube->getRMax() * conv;
	}
	
	// Build the layer array
	Acts::LayerArrayCreator layArrayCreator;
	return layArrayCreator.layerArray(layers, min, max, Acts::BinningType::arbitrary, Acts::BinningValue::binR);
}

std::shared_ptr<Acts::TrackingVolume>
FW::GeoModelPixel::buildPixelBarrel(GeoVPhysVol const* pdBarrel) const
{
	// Find the transformation
	std::shared_ptr<const Acts::Transform3D> transformation = nullptr;
	if(dynamic_cast<GeoFullPhysVol const*>(pdBarrel))
	{
		transformation = std::make_shared<Acts::Transform3D>(dynamic_cast<GeoVFullPhysVol const*>(pdBarrel)->getDefAbsoluteTransform());
	}

	// Look the bounds up
	std::shared_ptr<const Acts::VolumeBounds> bounds = nullptr;
	if(dynamic_cast<GeoTube const*>(pdBarrel->getLogVol()->getShape()))
	{
		auto tube = dynamic_cast<GeoTube const*>(pdBarrel->getLogVol()->getShape());
		bounds = std::make_shared<const Acts::CylinderVolumeBounds>(tube->getRMin(), tube->getRMax(), tube->getZHalfLength());
	}
	
	// Get the layers
	std::unique_ptr<const Acts::LayerArray> layArray = buildLayerArray(pdBarrel);

	// Build the volume
	return Acts::TrackingVolume::create(transformation, bounds, nullptr, std::move(layArray), {}, {}, {}, "ACTlaS::Pixel::Barrel");
}

std::shared_ptr<Acts::TrackingVolume>
FW::GeoModelPixel::buildPixelEndCap(GeoVPhysVol const* pdEndCap) const
{
	return nullptr;
}

std::shared_ptr<Acts::TrackingVolume>
FW::GeoModelPixel::buildPixel(GeoVPhysVol const* pd) const
{

	std::vector<std::shared_ptr<Acts::TrackingVolume>> volumes;

  // Walk over all children of the pixel volume
  unsigned int nChildren = pd->getNChildVols();
  for (unsigned int i = 0; i < nChildren; i++) {
	  GeoVPhysVol const* child = &(*(pd->getChildVol(i)));
	  // Test if it is a layer which has sensitive surfaces
	  if(m_volumeKeys.find(child->getLogVol()->getName()) != m_volumeKeys.end())
	  {
		  if(child->getLogVol()->getName() == "Barrel")
		  {
			volumes.push_back(buildPixelBarrel(child));
		  }
		  else
			if(child->getLogVol()->getName() == "EndCap")
			{
				volumes.push_back(buildPixelEndCap(child));
			}
	  }
  }
  return nullptr;
}