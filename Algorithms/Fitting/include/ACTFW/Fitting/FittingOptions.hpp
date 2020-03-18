// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include "ACTFW/Utilities/Options.hpp"
#include "FittingAlgorithm.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  template <typename aopt_t>
  void
  addFittingOptions(aopt_t& opt)
  {
    opt.add_options()("num-fitting-threads",
                      po::value<int>()->default_value(1),
                      "The number of threads used for fitting; default is 1; "
                      "negative value for automatic.")(
        "sorted-tracks",
        po::value<bool>()->default_value(0),
        "Decide if to sort tracks. Deafult is false.");
  }

  FW::FittingAlgorithm::Config
  readFittingConfig(const boost::program_options::variables_map& vm)
  {
    FittingAlgorithm::Config cfg;
    if (not vm["num-fitting-threads"].empty())
      cfg.numFittingThreads = vm["num-fitting-threads"].as<int>();
    cfg.sortedTracks = vm["sorted-tracks"].as<bool>();
    return cfg;
  }

}  // namespace Options
}  // namespace FW
