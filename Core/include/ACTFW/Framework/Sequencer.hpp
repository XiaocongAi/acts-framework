// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>
#include <tbb/task_scheduler_init.h>

#include <Acts/Utilities/Logger.hpp>

#include "ACTFW/Framework/IAlgorithm.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/IWriter.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// @class  Sequencer
///
/// This is the backbone of the mini framework, it initializes all algorithms,
/// calls execute per event and deals with the event store.
///
class Sequencer
{
public:
  struct Config
  {
    /// job store logging level
    Acts::Logging::Level jobStoreLogLevel = Acts::Logging::INFO;
    /// event store logging level
    Acts::Logging::Level eventStoreLogLevel = Acts::Logging::INFO;
  };

  /// Constructor
  ///
  /// @param cfg is the configuration object
  Sequencer(const Config&                       cfg,
            std::unique_ptr<const Acts::Logger> logger
            = Acts::getDefaultLogger("Sequencer", Acts::Logging::INFO));

  /// Add a service to the set of services.
  ///
  /// @throws std::invalid_argument if the service is NULL.
  void
  addService(std::shared_ptr<IService> service);
  /// Add a context decorator to the set of context decorators.
  ///
  /// @throws std::invalid_argument if the decorator is NULL.
  void
  addContextDecorator(std::shared_ptr<IContextDecorator> decorator);
  /// Add a reader to the set of readers.
  ///
  /// @throws std::invalid_argument if the reader is NULL.
  void
  addReader(std::shared_ptr<IReader> reader);
  /// Append an algorithm to the sequence of algorithms.
  ///
  /// @throws std::invalid_argument if the algorithm is NULL.
  void
  addAlgorithm(std::shared_ptr<IAlgorithm> algorithm);
  /// Add a writer to the set of writers.
  ///
  /// @throws std::invalid_argument if the writer is NULL.
  void
  addWriter(std::shared_ptr<IWriter> writer);

  /// Run the event loop over the given number of events.
  ///
  /// @param events (optional) Number of events to process
  /// @note This parameter is optional when input is read from a file. In this
  /// scenario, leaving it unset will process events until the end of the file,
  /// and setting it will put an upper bound on the number of events to be
  /// processed.
  /// @param skip Number of events to skip before processing
  ///
  /// This will run all configured algorithms for each event, potentially in
  /// parallel, then invoke the endRun hook of writers and services.
  ProcessCode
  run(boost::optional<size_t> events, size_t skip = 0);

private:
  std::vector<std::shared_ptr<IService>>          m_services;
  std::vector<std::shared_ptr<IContextDecorator>> m_decorators;
  std::vector<std::shared_ptr<IReader>>           m_readers;
  std::vector<std::shared_ptr<IAlgorithm>>        m_algorithms;
  std::vector<std::shared_ptr<IWriter>>           m_writers;
  Config                                          m_cfg;
  std::unique_ptr<const Acts::Logger>             m_logger;
  tbb::task_scheduler_init                        m_tbb_init;

  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW
