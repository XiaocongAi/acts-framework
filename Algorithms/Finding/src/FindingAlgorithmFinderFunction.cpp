// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Finding/FindingAlgorithm.hpp"

#include <random>
#include <stdexcept>

#include <Acts/Fitter/GainMatrixSmoother.hpp>
#include <Acts/Fitter/GainMatrixUpdater.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>

#include "ACTFW/Plugins/BField/ScalableBField.hpp"

namespace {
template <typename Finder>
struct FinderFunctionImpl
{
  Finder finder;

  FinderFunctionImpl(Finder&& f) : finder(std::move(f)) {}

  FW::FindingAlgorithm::FinderResult
  operator()(const std::vector<FW::Data::SimSourceLink>& sourceLinks,
             const FW::TrackParameters&                  initialParameters,
             const Acts::TrackFinderOptions<Acts::CKFSourceLinkSelector>&
                 options) const
  {
    return finder.findTracks(sourceLinks, initialParameters, options);
  };
};
}  // namespace

FW::FindingAlgorithm::FinderFunction
FW::FindingAlgorithm::makeFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    Options::BFieldVariant                        magneticField,
    Acts::Logging::Level                          lvl)
{
  using Updater  = Acts::GainMatrixUpdater<Acts::BoundParameters>;
  using Smoother = Acts::GainMatrixSmoother<Acts::BoundParameters>;

  // unpack the magnetic field variant and instantiate the corresponding finder.
  return std::visit(
      [trackingGeometry, lvl](auto&& inputField) -> FinderFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper       = Acts::EigenStepper<MagneticField>;
        using Navigator     = Acts::Navigator;
        using Propagator    = Acts::Propagator<Stepper, Navigator>;
        using SLS           = Acts::CKFSourceLinkSelector;
        using Finder = Acts::TrackFinder<Propagator, Updater, Smoother, SLS>;

        // construct all components for the finder
        MagneticField field(std::move(inputField));
        Stepper       stepper(std::move(field));
        Navigator     navigator(trackingGeometry);
        navigator.resolvePassive   = false;
        navigator.resolveMaterial  = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Finder     finder(std::move(propagator),
                      Acts::getDefaultLogger("TrackFinder", lvl));

        // build the finder functions. owns the finder object.
        return FinderFunctionImpl<Finder>(std::move(finder));
      },
      std::move(magneticField));
}
