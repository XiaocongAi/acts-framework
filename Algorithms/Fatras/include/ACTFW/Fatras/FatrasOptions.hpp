// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>

#include <boost/program_options.hpp>

#include "ACTFW/Fatras/FatrasAlgorithm.hpp"
#include "ACTFW/Utilities/OptionsFwd.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Physics/StandardPhysicsLists.hpp"

namespace FW {
namespace Options {

  /// Add Fatras options.
  ///
  /// @param desc The options description to add options to
  void
  addFatrasOptions(Description& desc);

  /// Read Fatras options to create the algorithm config.
  ///
  /// @tparam simulator_t type of the simulation kernel
  /// @param vars         the variables to read from
  /// @param simulator    the simulation kernel
  template <typename simulator_t>
  typename FatrasAlgorithm<simulator_t>::Config
  readFatrasConfig(const Variables& variables, simulator_t&& simulator)
  {
    using namespace Acts::UnitLiterals;
    using Config = typename FatrasAlgorithm<simulator_t>::Config;
    using PMin   = ActsFatras::Min<ActsFatras::Casts::P>;

    Config cfg(std::forward<simulator_t>(simulator));

    // simulation particle cuts both for input and physics lists output
    const auto pmin = variables["fatras-pmin-gev"].as<double>() * 1_GeV;
    cfg.simulator.selectCharged.template get<PMin>().valMin = pmin;
    cfg.simulator.selectNeutral.template get<PMin>().valMin = pmin;
    cfg.simulator.charged.physics
        = ActsFatras::makeChargedElectroMagneticPhysicsList(pmin);

    // all physics process are enabled by default
    if (not variables["fatras-em-scattering"].as<bool>()) {
      cfg.simulator.charged.physics
          .template disable<ActsFatras::detail::StandardScattering>();
    }
    if (not variables["fatras-em-ionisation"].as<bool>()) {
      cfg.simulator.charged.physics
          .template disable<ActsFatras::detail::StandardBetheBloch>();
    }
    if (not variables["fatras-em-radiation"].as<bool>()) {
      cfg.simulator.charged.physics
          .template disable<ActsFatras::detail::StandardBetheHeitler>();
    }

    // select hit surfaces for charged particles
    const auto hits = variables["fatras-hits"].as<std::string>();
    if (hits == "sensitive") {
      simulator.charged.selectHitSurface.sensitive == true;
      simulator.charged.selectHitSurface.material == false;
      simulator.charged.selectHitSurface.passive == false;
    } else if (hits == "material") {
      simulator.charged.selectHitSurface.sensitive == false;
      simulator.charged.selectHitSurface.material == true;
      simulator.charged.selectHitSurface.passive == false;
    } else if (hits == "all") {
      simulator.charged.selectHitSurface.sensitive == false;
      simulator.charged.selectHitSurface.material == false;
      // this includes sensitive and material surfaces
      simulator.charged.selectHitSurface.passive == true;
    } else {
      // none or unknown type -> record nothing
      simulator.charged.selectHitSurface.sensitive == false;
      simulator.charged.selectHitSurface.material == false;
      simulator.charged.selectHitSurface.passive == false;
    }

    return cfg;
  }

}  // namespace Options
}  // namespace FW
