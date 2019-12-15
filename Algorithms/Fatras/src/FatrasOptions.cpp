// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Fatras/FatrasOptions.hpp"

#include <string>

void
FW::Options::addFatrasOptions(FW::Options::Description& desc)
{
  using boost::program_options::value;

  desc.add_options()("fatras-sim-particles",
                     value<std::string>()->default_value("fatras-particles"),
                     "The collection of simulated particles.")(
      "fatras-sim-hits",
      value<std::string>()->default_value("fatras-hits"),
      "The collection of simulated hits")(
      "fatras-em-ionisation",
      value<bool>()->default_value(true),
      "Switch on ionisiation loss of charged particles")(
      "fatras-em-radiation",
      value<bool>()->default_value(false),
      "Switch on radiation for charged particles")(
      "fatras-em-scattering",
      value<bool>()->default_value(true),
      "Switch on multiple scattering")("fatras-em-conversions",
                                       value<bool>()->default_value(false),
                                       "Switch on gamma conversions")(
      "fatras-had-interaction",
      value<bool>()->default_value(false),
      "Switch on hadronic interaction")("fatras-debug-output",
                                        value<bool>()->default_value(false),
                                        "Switch on debug output on/off");
}
