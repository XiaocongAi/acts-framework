// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <iostream>
#include "FittingAlgorithm.hpp"

using namespace boost::program_options;

namespace FW {

namespace Options {

  void
  addFittingOptions(boost::program_options::options_description& opt)
  {
    opt.add_options()("num-fitting-threads",
                      value<int>()->default_value(1),
                      "The number of threads used for fitting; default is 1; "
                      "negative value for automatic.")(
        "sorted-tracks",
        value<bool>()->default_value(false),
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
